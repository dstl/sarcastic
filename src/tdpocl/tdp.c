/***************************************************************************
 *
 *       Module:    tdp.c
 *      Program:    tdpocl
 *   Created by:    Darren Muff on 12/03/2013.
 *                  Copyright (c) 2013 [Dstl]. All rights reserved.
 *
 *   Description:
 *   Main routine for Time Domain Processon.
 *    This module is responsible for collecting user input, determining
 *    memory requirements and then farming out the work to key tdp_core
 *    processing modules.
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  18/09/2014
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 *
 ***************************************************************************/

#include "tdp.h"
#include "tdpcore.h"
#include <unistd.h>
#include "TxPowerPerRay.h"
#define DEBUGLVL 10
#define RCSCORRECT
//#define     INLINEKERNELCODE
#if defined INLINEKERNELCODE
#include "tdpKernel.h"
#else
// Path to kernel source code
static const char *kernelCode     = KERNELDIR"/tdpKernel.cl";
#endif

// To better optimise GPU cards, the following structure
// allows you to set custom sizes for the device localWorkSize.
// A negative value will trigger a calculation to try and work it out
// for you.
//
const int    ndevdata = 5 ;
const struct oclLWS devdata[] = {
    { "Tesla C2075",                                     32, 1 },
    { "Intel(R) Xeon(R) CPU           E5620  @ 2.40GHz", 0, 0  },
    { "Intel(R) Core(TM) i7-3820QM CPU @ 2.70GHz",       8, 1  },
    { "GeForce GT 650M",                                 16, 8 },
    { "Tesla K10.G2.8GB",                                32,2  }

} ;

int create_surface(SPVector sc_pos, SPVector sensor_pos, char * surfaceFilename, int nx, int ny, double surface_sx, double surface_sy, SPStatus * status);
KSpaceSupportImage * load_kspace_trimming(FILE * fp, int process_reference, int nx, int ny, int surface_ny, int surface_oy, int npulse, int nsamp);
void *  do_work(void * info);

int main(int argc, char *argv[])
{
    SPImageLineSaveInfo save_info;
    SPStatus status;
    SPVector corner[4];
    SPImage image;
    SPImage cphd;
    SPCmplx tmp;
    SPCmplx pcorr;
    
    double phs;
    double EIRP = DEFAULT_EIRP ;
    
    char * cphdFilename;
    char * outputFilename;
    char * t_outname = NULL;
    char * old_t_outname = NULL;
    char * surfaceFilename;
    
    CPHDHeader hdr;
    int startPulse;
    int nPulses;
    int x;
    int y;
    int surface_nx;
    int surface_ny;
    int surface_xstart ;
    int surface_ystart ;
    int chunkOffset ;
    int scan =0;
    int iPulse;
    int topotrim=0;
    int USEGPU=0;
    TDPCoreStruct tdpCore_s;
    
    banner();
        
    /* This sets the debug level used (some functions in mtlib print out debug info if this is set higher)
     in this program (and all the threads it creates) */
    status.debug = DEBUGLVL;
    fftwf_init_threads();
    im_init_lib(&status, "tdpocl", argc, argv);
    
    // Initialise the fftw library to use threads
    //
    fftwf_plan_with_nthreads((int)sysconf(_SC_NPROCESSORS_ONLN));
    
    // Process User Input
    //
    getUserInput(&USEGPU, &cphdFilename, &surfaceFilename, &surface_nx, &surface_ny, &surface_xstart, &surface_ystart,
                 &nPulses, &startPulse, &EIRP, &outputFilename, &status) ;
    CHECK_STATUS_NON_PTR(status) ;
    
    
    // Initialise OpenCL and load the relevent information into the
    // platform structure. OCL tasks will use this later
    //
    cl_int err;
    cl_uint ndevs;
    OCLPlatform platform;
       
    platform.clSelectedPlatformID = NULL;
    CL_CHECK(oclGetPlatformID (&platform, &status));
    
    oclPrintPlatform(platform);
    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform.clSelectedPlatformID,
        0
    };
    platform.props = props ;
   
    // Find number of devices on this platform
    //
    if (USEGPU){
        err = oclGetNamedGPUDevices(&platform, "NVIDIA", "", &platform.device_ids, &ndevs, &status);
        if(err == CL_SUCCESS && ndevs != 0 ){
            char cbuf[1024];
            printf("GPU DEVICES                 : %d\n",ndevs);
            for(int d=0; d<ndevs; d++){
                printf("DEVICE                      : %d\n",d);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_NAME, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE NAME               : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_VENDOR, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE VENDOR             : %s\n",cbuf);
                
            }
        }else{
            printf("No GPU devices available with compute capability: %d.%d\n", GPUCAPABILITY_MAJOR, GPUCAPABILITY_MINOR);
            exit(1);
        }
    }else{
        err = oclGetCPUDevices(&platform, &platform.device_ids, &ndevs, &status);
        char cbuf[1024];
        if(err == CL_SUCCESS){
            printf("CPU DEVICES                 : %d\n",ndevs);
            for(int d=0; d<ndevs; d++){
                printf("DEVICE                      : %d\n",d);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_NAME, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE NAME               : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_VENDOR, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE VENDOR             : %s\n",cbuf);
                
            }
        }
    }
    platform.nDevs = ndevs ;
    
    // Now initialise the task structure that will perform the
    // image formation
    //
    oclCreateTask(&platform, ndevs, "tdpocl", kernelCode, &(tdpCore_s.task_s), &status);

    // Find out the maximum amount of memory that all OpenCL devices have
    // For this task
    //
    cl_ulong memSize,devMemSize;
    devMemSize = 100e9;    // big number 100 GB
    for (int dev=0; dev<tdpCore_s.task_s.devsToUse; dev++){
        err = clGetDeviceInfo(tdpCore_s.task_s.device_ids[dev], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &memSize, NULL);
        if(err != CL_SUCCESS){
            printf("Error [%d] : couldnt obtain device info\n",err);
        }
        devMemSize = (memSize < devMemSize) ? memSize : devMemSize;
    }
    long kernelsize = 1024*1024*368;  // 368 MBytes - determined through experimentation. Will change if you change the kernel
    devMemSize -= kernelsize ;

    printf("DEVICES USED                : %d\n",tdpCore_s.task_s.devsToUse);
   
    // Now that we know the number of devices to use allocate arrays to the processing
    // structures
    //
    SPImage *images, *surfaces;
    images             = (SPImage *)sp_malloc(sizeof(SPImage)*tdpCore_s.task_s.devsToUse);
    surfaces           = (SPImage *)sp_malloc(sizeof(SPImage)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.images   = images ;
    tdpCore_s.surfaces = surfaces ;
    
    // The kernel performance is heavily dependent upon the localworksize values specified.
    // These are best found by experimentation. Use the device name to see if LWS sizes have
    // been determined and if so then use them. If not then set the LWS values to be negative
    // which will flag the core routines to try and calculate the best values to use
    //
    
    tdpCore_s.localWorkSize[0] = 0 ;
    tdpCore_s.localWorkSize[1] = 0 ;
    for (int i=0; i<ndevdata; i++){
        char cbuf[1024];
        clGetDeviceInfo(platform.device_ids[0], CL_DEVICE_NAME, sizeof(cbuf), &cbuf, NULL);
        if(! strcmp(cbuf, devdata[i].devName) ){
            tdpCore_s.localWorkSize[0] = devdata[i].x ;
            tdpCore_s.localWorkSize[1] = devdata[i].y ;
        }
    }
    printf("LocalWorkSize set to : %ld, %ld\n",tdpCore_s.localWorkSize[0],tdpCore_s.localWorkSize[1]);
    
    // Reduce the surface size in x a little bit if required so that we can equally
    // divide the processing run amongst the devices.
    //
    while ( (surface_nx % tdpCore_s.task_s.devsToUse ) != 0 ) surface_nx-- ;
    
    // initialise cphd which is the phase history store for a block of pulses
    // Also initialise the task phd store in tdpCore_s.phd. This the block
    // of range compressed pulses to be processed by the ocl kernel.
    // Read in teh cphd header so that we can use the parameters later on.
    // Also create the final image. This will be written a part at a time during the loop below.
    //
    im_init(&tdpCore_s.phd, &status);
    readCPHDHeader(cphdFilename, &hdr, &status);
    im_init(&cphd, &status);
//    load_cphd(&cphd, &hdr, 200, 201, &status);
//    im_destroy(&cphd, &status);
    im_create(&image, ITYPE_CMPL_FLOAT, surface_nx, surface_ny, 1.0, 1.0, &status);
    CHECK_STATUS_NON_PTR(status);
    
    // t_outname is a temporary name used to store intermediate files.
    //
    t_outname = sp_calloc(strlen(outputFilename)+40, sizeof(char));
        
    
    //     We need to efficiently break the problem into manageable chunks.
    //     The philosophy is to form an image from only a block of the pulses. The block of pulses steps down
    //     through all the pulses and each image is coherently added to the next. In this way the resolution
    //     is progressively refined.
    //
    //     For each block of pulses we also need to break down the image surface into manageable chunks.
    //     For this we spilt the surface in the X-direction assigning an equal portion to each device.
    //     The amount of y-samples per device is determined by the available memory on each device.
    //
    //     The final complication is that the range compression for the block of pulses is performed on the
    //     host machine.
    //
  
    int nPulsesPerBlock ;
    int nPulseBlocks ;
    int nSampsPerPulse ;
    int nSurfaceLinesPerChunk, nSurfaceLinesThisChunk ;
    int nChunks ;
    int remainder ;
    cl_ulong memForPulse , memForAll , memAvailForSurfaceLine , devMemForSurfaceLine, memForRComp ;
    nSampsPerPulse = hdr.nsamp ;
    memForAll      = sizeof(float) * (NPOINTS * OVERSAMP + 1) ; // Interpolation Kernel
    memForPulse    = sizeof(SPCmplx) * nSampsPerPulse * RCOMPOVERSAMPLE  // Samples in pulse
                   + sizeof(SPVector)                   // Rx Vector
                   + sizeof(SPVector)                   // Tx Vector
                   + sizeof(double)                     // SRP ranges
                   + sizeof(double);                    // FxStep
   
    size_t cpumem = getMemorySize() ;
    memForRComp   = sizeof(SPCmplx) * nSampsPerPulse * RCOMPOVERSAMPLE ;  // a wild stab at the amount of space used by fftw to range comp on CPU
    long progsize = 4294967296; // guess prog size to be 4 gig
    
    // calculate the memory required for a single surface line. This is used later to determine the number of lines
    // that can be processed in a single chunk.
    // mem = size of surface locations on input + sizeof image formed for output / number of devices
    //
    devMemForSurfaceLine = ((sizeof(SPVector) *  surface_nx) + (sizeof(SPCmplx) * surface_nx)) / tdpCore_s.task_s.devsToUse ;
  
    // t is the number of pulses that can be stuffed onto the device (CPU or GPU) based upon device memory available
    // and the size of a single pulse (including required other data such as Rx,Tx and SRP positions
    //
    long t = (long)floor((double)((devMemSize - devMemForSurfaceLine) - memForAll) / (double)memForPulse );
   
    // s is the number of pulses that can be stuffed onto the device based upon the pulse's size and
    // the amount of memory required for range compression on the CPU
    //
    long s = (long)floor((double)((cpumem - progsize) / memForRComp));
    
    // If s is less than t (usually from long pulse lengths) then we use s for the pulsesPerBlock to ensure that
    // the fftf does not keep swapping in and out of memory
    // if t is less than s then we use this value.
    //
    t = (s < t) ? s : t ;
   
    nPulseBlocks    = (int)ceil((double)nPulses / (double)t) ;
    nPulsesPerBlock = nPulses / nPulseBlocks ;
    // Calculate the number of surface lines to process in a chunk
    //
    memAvailForSurfaceLine = devMemSize - memForAll - (nPulsesPerBlock * memForPulse) ;
    nSurfaceLinesPerChunk = (int)floor( memAvailForSurfaceLine / devMemForSurfaceLine ) ;
    if(nSurfaceLinesPerChunk > surface_ny)nSurfaceLinesPerChunk=surface_ny ;
    nChunks = (int)floor((double)surface_ny / (double)nSurfaceLinesPerChunk) ;
    remainder = surface_ny - (nChunks * nSurfaceLinesPerChunk) ;
    if (remainder > 0) {
        nChunks++;
    }
    // number of chunks is number of device cycles - ie number of whole chunks + remainder chunk
    //
    printf("-----------------------------------------------------------------------------------------\n");
    printf("                      Device Memory Summary\n");
    printf("Memory per Device                : %lld MByte\n",(unsigned long long)devMemSize/1024/1024);
    printf("CPU memory for pulse compression : %lld MByte\n",(unsigned long long)cpumem/1024/1024);
    printf("Number of Devices                : %d\n",tdpCore_s.task_s.devsToUse);
    printf("Number of threads                : %d\n",(int)sysconf(_SC_NPROCESSORS_ONLN)) ;
    printf("Number of pulse blocks           : %d\n",nPulseBlocks);
    printf("Number of Pulses per block       : %d\n",nPulsesPerBlock);
    printf("X surface elements               : %d\n",surface_nx);
    printf("X surface elements per device    : %d\n",surface_nx / tdpCore_s.task_s.devsToUse);
    printf("Y surface chunks                 : %d\n",nChunks);
    printf("Y surface rows per chunk         : %d\n",nSurfaceLinesPerChunk) ;
    printf("Last chunk rows                  : %d\n",remainder) ;
    printf("-----------------------------------------------------------------------------------------\n");

    // The following items are required throughout the program so initialise them here
    //
    tdpCore_s.RXs        = (SPVector *)sp_malloc(sizeof(SPVector)*nPulsesPerBlock) ;
    tdpCore_s.TXs        = (SPVector *)sp_malloc(sizeof(SPVector)*nPulsesPerBlock) ;
    tdpCore_s.FXStep     =    (float *)sp_malloc(sizeof(float)*nPulsesPerBlock) ;
    tdpCore_s.SRPRanges  =   (double *)sp_malloc(sizeof(double)*nPulsesPerBlock);
    
    // Allocate memory for items used by all devices
    //
    tdpCore_s.dPHD       = (cl_mem *)sp_malloc(sizeof(cl_mem)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.dTXs       = (cl_mem *)sp_malloc(sizeof(cl_mem)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.dRXs       = (cl_mem *)sp_malloc(sizeof(cl_mem)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.dSRPRanges = (cl_mem *)sp_malloc(sizeof(cl_mem)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.dFXStep    = (cl_mem *)sp_malloc(sizeof(cl_mem)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.dSurf      = (cl_mem *)sp_malloc(sizeof(cl_mem)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.dIm        = (cl_mem *)sp_malloc(sizeof(cl_mem)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.diKernel   = (cl_mem *)sp_malloc(sizeof(cl_mem)*tdpCore_s.task_s.devsToUse);
    tdpCore_s.fcent      = hdr.freq_centre ;
    tdpCore_s.phseSgn    = hdr.phaseSgn ;
    
    // build a sinc kernel and convert it to floats for quick GPU processing
    //
    double * ikernel_dbl = generate_sinc_kernel(OVERSAMP, NPOINTS);
    if (ikernel_dbl == NULL) {
        printf("Failed to create ikernel!\n");
        exit(614);
    }
    tdpCore_s.ikernel = (float *)sp_calloc(OVERSAMP * NPOINTS + 1, sizeof(float));
    if (tdpCore_s.ikernel == NULL) {
        printf("Failed to create ikernel!\n");
        exit(825);
    }
    for(int i = 0; i < (OVERSAMP * NPOINTS+1); i++) {
        tdpCore_s.ikernel[i] = ikernel_dbl[i];
    }
    free(ikernel_dbl);
    
    // Pre-calculate the SRP ranges for each pulse. This saves us calculation time in the image
    // formation kernel
    //
    double *SRPRanges =   (double *)sp_malloc(sizeof(double)  *nPulses);
    SPVector *TXs     = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses);
    SPVector *RXs     = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses);
    SPVector *SRP     = (SPVector *)sp_malloc(sizeof(SPVector)*nPulses);
    
    for (int i=0; i<nPulses; i++){
        TXs[i] = hdr.pulses[i+startPulse].sat_ps_tx ;
        RXs[i] = hdr.pulses[i+startPulse].sat_ps_rx ;
        SRP[i] = hdr.pulses[i+startPulse].srp ;
    }
    
    // Now initialise the task structure that will perform the
    // image formation
    //
    oclSRPRanges(platform, nPulses, TXs, RXs, SRP, &SRPRanges, &status);
    
    free(TXs);
    free(RXs);
    free(SRP);
    
    // Initialisation done - start the clock
    //
    struct timeval tv ;
    gettimeofday(&tv, NULL) ;
    double start = (double)tv.tv_sec+((double)tv.tv_usec / 1.0e6);
    
    // We go through the raw data in blocks of pulses. (nPulsesPerBlock is the number of pulses in a block.)
    // For each block of pulses, we then create the number of threads to work on that block.
    // Each thread then works (pixels by pixel) on a chunk of the image which is image_width / devsToUses wide, and
    // nSurfaceLinesPerChunk high. When all lines in a chunk are completed the next chunk is processed and finally the
    //  remainder lines are processed. Then the next block of pulses is tackled. */
    //
    for(int blk = 0; blk < nPulseBlocks; blk++) {
        printf("Processing pulses %6d - %6d out of %6d [%d%%]\n", blk*nPulsesPerBlock, (blk+1)*nPulsesPerBlock-1, nPulses,100*blk*nPulsesPerBlock/nPulses);
        
        // Load the raw data and range compress it
        // This uses the fftw3 library and multiple threads as set up in line 89
        //
        load_cphd(&cphd, &hdr, blk*nPulsesPerBlock+startPulse, blk*nPulsesPerBlock+nPulsesPerBlock+startPulse, &status);
        im_zpad(&cphd, (int)(cphd.nx * RCOMPOVERSAMPLE), cphd.ny, &tdpCore_s.phd, &status);
        im_fftw(&tdpCore_s.phd, FFT_X_ONLY+REV, &status);
        im_circshift(&tdpCore_s.phd, tdpCore_s.phd.nx/2, 0, &status);

        // Do the phase correct thats required because the start frequency doesn't stay the same pulse to pulse
        //
        for(y = 0; y < tdpCore_s.phd.ny; y++) {
            iPulse = blk*nPulsesPerBlock+startPulse+y ;
            phs = (((hdr.pulses[iPulse].fx0 - hdr.freq_centre) / hdr.pulses[iPulse].fx_step_size)) * 2.0 * M_PI / (double)tdpCore_s.phd.nx;
            
            for(x = 0; x < tdpCore_s.phd.nx; x++) {
                pcorr.r = hdr.pulses[iPulse].amp_sf0 * cos(phs * (x - tdpCore_s.phd.nx/2));
                pcorr.i = hdr.pulses[iPulse].amp_sf0 * sin(phs * (x - tdpCore_s.phd.nx/2));
                
                CMPLX_MULT(tdpCore_s.phd.data.cmpl_f[x + y * tdpCore_s.phd.nx],  pcorr, tmp);
                tdpCore_s.phd.data.cmpl_f[x + y * tdpCore_s.phd.nx] = tmp;
            }
            
            // also load up tx,rx and srp positions
            //
            tdpCore_s.RXs[y]       = hdr.pulses[iPulse].sat_ps_rx ;
            tdpCore_s.TXs[y]       = hdr.pulses[iPulse].sat_ps_tx ;
            tdpCore_s.FXStep[y]    = hdr.pulses[iPulse].fx_step_size ;
            tdpCore_s.SRPRanges[y] = SRPRanges[y] ;
        }
        
        if (topotrim) {
            /* KSpace trim mode needs the data in the frequency domain and in the same pixel position it
             came in down from the sensor (hence why the shift is im.nx/2 rather than tdpCore_s.phd->nx/2 */
            im_fftw(&tdpCore_s.phd, FFT_X_ONLY+FWD+NOSCALE, &status);
            im_circshift(&tdpCore_s.phd, tdpCore_s.phd.nx/2, 0, &status);
        }
        
        if (blk != 0) {
            if (old_t_outname) free(old_t_outname);
            old_t_outname = strdup(t_outname);
        }
        snprintf(t_outname, strlen(outputFilename)+40, "%s_pulse%d.dat", outputFilename, blk*nPulsesPerBlock+startPulse);
        im_save_line_init(&image, t_outname, &save_info, &status);
        
        // Loop through the surface chunks for this pulse block
        //
        for (int iChunk=0; iChunk < nChunks; iChunk++){
                       
            if(iChunk == 0){
                tdpCore_s.run = TDPINIT ;
            }
            if(iChunk == nChunks-1){
                tdpCore_s.run += TDPCLOSE ;
            }
            chunkOffset = iChunk * nSurfaceLinesPerChunk ;

            if(iChunk == nChunks-1 && remainder !=0 ){
                nSurfaceLinesThisChunk = remainder ;
            }else{
                nSurfaceLinesThisChunk = nSurfaceLinesPerChunk ;
            }
            
            // Split the surface and the image in the X-direction across the devices
            //
           
            for (int dev = 0; dev<tdpCore_s.task_s.devsToUse; dev++) {
                im_create(&(surfaces[dev]), ITYPE_VECTOR    , surface_nx/tdpCore_s.task_s.devsToUse, nSurfaceLinesThisChunk, 1.0, 1.0, &status);
                im_load_subset(&(surfaces[dev]), surfaceFilename, dev * surface_nx/tdpCore_s.task_s.devsToUse, chunkOffset, &status);
                
                im_create(&(images[dev])  , ITYPE_CMPL_FLOAT, surface_nx/tdpCore_s.task_s.devsToUse, nSurfaceLinesThisChunk, 1.0, 1.0, &status);
                if (blk != 0) {
                    im_load_subset(&(images[dev]), old_t_outname, dev*surface_nx/tdpCore_s.task_s.devsToUse, chunkOffset, &status);
                } else {
                    im_zero_data(&(images[dev]), &status);
                }
                CHECK_STATUS_NON_PTR(status);

            }
            if (chunkOffset == 0) {
                corner[0] = surfaces[0].data.vect[0];
                corner[1] = surfaces[tdpCore_s.task_s.devsToUse-1].data.vect[surfaces[tdpCore_s.task_s.devsToUse-1].nx - 1];
            }
            if (chunkOffset + surfaces[0].ny == surface_ny) {
                corner[3] = surfaces[0].data.vect[(surfaces[0].ny - 1) * surfaces[0].nx];
                corner[2] = surfaces[tdpCore_s.task_s.devsToUse-1].data.vect[(surfaces[tdpCore_s.task_s.devsToUse-1].ny * surfaces[tdpCore_s.task_s.devsToUse-1].nx) - 1];
            }
            
            // At this point there are various options to form the image chunks. To make things cleaner we will
            // have each variant in a seperate submodule
            // Each chunk is split across devsToUse OpenCL devices
            //
            if(topotrim) {
                // Call topotrim variant. This allows each pixel to have its own unique pulses and samples
                //
                /*if (kspacebnds) {
                    free(kspacebnds->data);
                    free(kspacebnds);
                }
                kspacebnds = load_kspace_trimming(topotrim_fp, process_reference, surface_nx, surface_ny, surface_ny, chunkOffset, hdr.num_azi, hdr.nsamp);
                rewind(topotrim_fp) ;*/
            }else if(scan){
                
            }else{
                // perform image formation. Take phase history data from tdpCore_s.phd and use it to
                // focus pixels for each location defined by tdpCore_s.surface. Each pixel is then stored in tdpCore_s.im which
                // has to be the same size as surface. Values in tdpCore_s.im when tdpcore is called are coherently added to
                // new values.
                
                tdpcore( &tdpCore_s) ;
            }
            
            // Combine the images from the devices back into a single image
            //
            SPImage im;
            im_create (&im,      ITYPE_CMPL_FLOAT,  surface_nx, nSurfaceLinesThisChunk, 1.0, 1.0, &status);
            for (int dev = 0; dev<tdpCore_s.task_s.devsToUse; dev++) {
#ifdef RCSCORRECT
                printf("  Performing RCS Correction...");
                double TxRange, RxRange, rcs, m, p, PtG, gainRx, maxrcs = -1e6;
                SPVector txtopix_v,rxtopix_v;
                PtG = TxPowerPerRay(&hdr, &gainRx);
                for(int y=0; y<surfaces[dev].ny; y++){
                    for(int x=0; x<surfaces[dev].nx; x++){
                        VECT_SUB(surfaces[dev].data.vect[y*surfaces[dev].nx+x], hdr.pulses[startPulse+nPulses/2].sat_ps_tx, txtopix_v);
                        VECT_SUB(surfaces[dev].data.vect[y*surfaces[dev].nx+x], hdr.pulses[startPulse+nPulses/2].sat_ps_rx, rxtopix_v);
                        TxRange  = VECT_MAG(txtopix_v);
                        RxRange  = VECT_MAG(rxtopix_v);
                        m        = CMPLX_MAG(images[dev].data.cmpl_f[y*surfaces[dev].nx+x]) ;  // E in volts
                        p        = CMPLX_PHASE(images[dev].data.cmpl_f[y*surfaces[dev].nx+x]); // phase
                        rcs      = RCOMPOVERSAMPLE * m * m * (4 * SIPC_pi * TxRange * TxRange) * (4 * SIPC_pi * RxRange * RxRange) / (PtG * PtG) ;
                        maxrcs   = (rcs > maxrcs ) ? rcs : maxrcs ;
                        images[dev].data.cmpl_f[y*surfaces[dev].nx+x].r = rcs * cos(p);
                        images[dev].data.cmpl_f[y*surfaces[dev].nx+x].i = rcs * sin(p);
                    }
                }
                printf("Done. \n  Max RCS in scene is %8.2f m^2 (%5.0f dBm^2)\n",maxrcs, 10*log10(maxrcs));
#endif
                im_insert(&(images[dev]), dev * images[dev].nx, 0, &im, &status);
            }
            
            // Save out the intermediate image file
            //
            im_save_line(&im, &save_info, &status);
            im_destroy(&im,   &status);
            for (int dev = 0; dev<tdpCore_s.task_s.devsToUse; dev++) {
                im_destroy(&(surfaces[dev]), &status);
                im_destroy(&(images[dev])  , &status);
            }
        }
        
        im_save_line_close_with_geo_ecf(&save_info, corner[0], corner[1], corner[2], corner[3], &status);
        destroy_cphd_pulses(&hdr, blk*nPulsesPerBlock+startPulse, blk*nPulsesPerBlock + nPulsesPerBlock +startPulse, &status);
        im_destroy(&cphd, &status) ;
        im_destroy(&tdpCore_s.phd, &status) ;
    }
    
   
    for(int dev=0; dev<tdpCore_s.task_s.devsToUse; dev++){
        clReleaseContext(tdpCore_s.task_s.contexts[dev]);
        clReleaseProgram(tdpCore_s.task_s.programs[dev]);
        clReleaseKernel(tdpCore_s.task_s.kernels[dev]);
    }
       
    free ( tdpCore_s.ikernel  );
    free ( old_t_outname      );
    free ( tdpCore_s.RXs      );
    free ( tdpCore_s.TXs      );
    free ( tdpCore_s.FXStep   );
    free ( tdpCore_s.dPHD     );
    free ( tdpCore_s.dTXs     );
    free ( tdpCore_s.dRXs     );
    free ( tdpCore_s.dFXStep  );
    free ( tdpCore_s.dSurf    );
    free ( tdpCore_s.dIm      );
    free ( tdpCore_s.diKernel );
    free ( images             );
    free ( surfaces           );
    
    rename(t_outname, outputFilename);
    char * tmp1 = sp_calloc(strlen(t_outname)+40, 1);
    char * tmp2 = sp_calloc(strlen(outputFilename)+40, 1);
    sprintf(tmp1, "%s.hdr", t_outname);
    sprintf(tmp2, "%s.hdr", outputFilename);
    rename(tmp1, tmp2);
    free(tmp1);
    free(tmp2);
    
    gettimeofday(&tv, NULL);
    double end = (double)tv.tv_sec+((double)tv.tv_usec / 1.0e6);
    printf("Thinking completed in %f secs\n",end-start);
    printf("  pixels                      : %d\n",surface_nx*surface_ny);
    printf("  pixels per second           : %3.2f Pixels/Sec\n",(double)(surface_nx*surface_ny)/(end-start));
    printf("  Seconds per pixel           : %3.2e Sec/Pixel\n",(end-start)/((double)(surface_nx*surface_ny)));
    printf("  Seconds per pixel per pulse : %3.2e Sec/pixel/pulse\n", ((end-start)/((double)(surface_nx*surface_ny)))/nPulses);
    printf("  Pulse pixels per second     : %3.2f MPPpS\n", (double)(surface_nx*surface_ny)*nPulses/(end-start)/1e6);
    
    
    free(outputFilename);
    free(t_outname);
    
    im_destroy(&image, &status);
    im_destroy(&tdpCore_s.phd, &status);
    
    im_close_lib(&status);
    
    return (0);
}

