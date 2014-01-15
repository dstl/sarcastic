/***************************************************************************
 *
 *       Module:    threadCore.c
 *      Program:    GOSS
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      This function runs on a single thread and controls a single OpenCL device
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  15/01/2014
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
#include "GOSS.h"
#include <fftw3.h>

#define FASTPULSEBUILD
#define NPOINTS (32)
#define OVERSAMP (512)

static char *kernelCodePath = "/Users/darren/Development/GOSS/GOSS/GOSSKernelCode.cl" ;
void packSinc(SPCmplx point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel);
static void ham1dx(double * data, int nx);
double * sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing);

void * devPulseBlock ( void * threadArg ) {
    
    struct threadData *td;
    td = (struct threadData *) threadArg;
    
    cl_device_id devId ;
    cl_context context ;
    cl_program program ;
    cl_kernel kernel ;
    cl_command_queue commandQ ;
    cl_int err ;
    cl_mem dTriangles, dTextures, dKdTree, dtriListData, dtriListPtrs, drnp ;
    rangeAndPower * rnp ;
    SPVector aimdir ;

    
    int nAzBeam             = td->nAzBeam ;
    int nElBeam             = td->nElBeam ;
    SPVector RxPos;
    SPVector TxPos;
    double dAz              = td->dAz ;
    double dEl              = td->dEl ;
    double raySolidAng      = td->raySolidAng ;
    double TxPowPerRay      = td->TxPowPerRay ;
    AABB SceneBoundingBox   = td->SceneBoundingBox ;
    double Aeff             = td->Aeff ;
    int bounceToShow        = td->bounceToShow ;
    int debug               = td->debug;
    int debugX              = td->debugX;
    int debugY              = td->debugY;
    int nrnpItems;
    double derampRange ;
#ifdef FASTPULSEBUILD
    SPCmplx targ, pcorr, tmp;
    SPImage pulseLine ;
    double rangeLabel, phasecorr,resolution,sampSpacing, *ikernel ;
    double bandwidth = td->chirpRate * td->pulseDuration ;
    fftwf_init_threads();
    printf("Using FAST pulse building routines...\n");
#else
    double freqSampsToCentreStart, ADCSampsToTargStart ;
    double A = 4.0*SIPC_pi*td->chirpRate / (SIPC_c * td->ADRate);
    double B = 4.0*SIPC_pi*td->oneOverLambda ;
#endif

    int tid ;
    tid   = td->devIndex ;
    devId = td->platform.device_ids[tid] ;
    context = CL_CHECK_ERR(clCreateContext(td->platform.props, 1, &devId, NULL, NULL, &_err));
    program = CL_CHECK_ERR(clCreateProgramWithSource(context, 1, (const char **) &kernelCodePath, NULL, &_err));
    
    char compilerOptions[255];
    sprintf(compilerOptions, "-Werror -D MAXBOUNCES=%d",MAXBOUNCES) ;
    if (debug) {
        char dbgxy[128];
        sprintf(dbgxy, " -D DEBUG=%d -D DBGX=%d -D DBGY=%d",debug,debugX,debugY);
        strcat(compilerOptions, dbgxy);
    }
    
    err = clBuildProgram(program, 0, NULL, compilerOptions, NULL, NULL);
    if (err != CL_SUCCESS){
        size_t len;
        char buffer[32768];
        printf("[thread:%d], Error: Failed to build program executable!\n",tid);
        clGetProgramBuildInfo(program, devId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("err: %d. Buffer:\n",err);
        printf("%s\n", buffer);
        exit(-3);
    }
    
    kernel   = CL_CHECK_ERR(clCreateKernel(program, "rayTraceBeam", &_err));
    commandQ = CL_CHECK_ERR(clCreateCommandQueue(context, devId, 0, &_err));
    
    
    // Calculate global and local work sizes
    //
    int tx,ty ;
    size_t localWorkSize [2] ;
    size_t globalWorkSize[2] ;
    globalWorkSize[0] = nAzBeam ;
    globalWorkSize[1] = nElBeam ;
    best2DWorkSize(kernel, devId, globalWorkSize[0], globalWorkSize[1], &tx, &ty, &td->status) ;
    localWorkSize[0]  = tx ;
    localWorkSize[1]  = ty ;
    
    // Create device buffers
    //
    dTriangles   = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Triangle)*td->nTriangles, NULL, &_err));
    dTextures    = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Texture)*td->nTextures, NULL, &_err));
    dKdTree      = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(KdData)*td->nTreeNodes, NULL, &_err));
    dtriListData = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, td->triListDataSize, NULL, &_err));
    dtriListPtrs = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int)*td->nLeaves, NULL, &_err));
    drnp = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(rangeAndPower)*nAzBeam*nElBeam*MAXBOUNCES, NULL, &_err));

    // Load the device buffers
    //
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dTriangles,   CL_TRUE, 0, sizeof(Triangle)*td->nTriangles, td->Triangles,    0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dTextures,    CL_TRUE, 0, sizeof(Texture)*td->nTextures,td->Textures,        0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dKdTree,      CL_TRUE, 0, sizeof(KdData)*td->nTreeNodes, td->KdTree,         0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtriListData, CL_TRUE, 0, td->triListDataSize, td->triListData,              0, NULL, NULL);
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtriListPtrs, CL_TRUE, 0, sizeof(int)*td->nLeaves, td->triPtrs,              0, NULL, NULL)));
    
    Timer threadTimer ;
    SPStatus status;
    startTimer(&threadTimer, &status) ;
    int hrs,min,sec;
    time_t current,complete;
    struct tm *p;
    int reportN = 500 ;
    char ct[1000];
    double pCentDone ;
    double sexToGo ;
    int nx;
    nx = (int)td->phd->nx ;

#ifdef FASTPULSEBUILD

    // build a sinc kernel and convert it to floats for quick GPU processing
    //
    
    resolution = SIPC_c / (2 * bandwidth ) ;
    sampSpacing = SIPC_c /  (2.0 * (nx/(td->pulseDuration * td->ADRate)) * bandwidth) ;
    ikernel = sinc_kernel(OVERSAMP, NPOINTS, resolution, sampSpacing);
    if (ikernel == NULL) {
        printf("Failed to create ikernel!\n");
        exit(614);
    }
    
#endif
    
    // **** loop  start here
    //
    for (int pulse=0; pulse<td->nPulses; pulse++){
        int pulseIndex = (tid * td->nPulses) + pulse ;
        if ( td->devIndex == 0 ) {
            if( pulse % reportN == 0 && td->nPulses != 1){
                pCentDone = 100.0*pulse/td->nPulses ;
                printf("Processing pulses %6d - %6d out of %6d [%2d%%]",  pulse*td->nThreads, (pulse+reportN)*td->nThreads, td->nPulses*td->nThreads,(int)pCentDone);
                if(pulse != 0 ){
                    current  = time(NULL);
                    sexToGo  = estimatedTimeRemaining(&threadTimer, pCentDone, &status);
                    hrs      = (int)floor(sexToGo/60/60) ;
                    min      = (int)floor(sexToGo / 60) - (hrs*60);
                    sec      = (int)sexToGo % 60;
                    complete = current + sexToGo ;
                    p        = localtime(&complete) ;
                    strftime(ct, 1000, "%a %b %d %H:%M", p);
                    printf("  ETC %s (in %2.0dh:%2.0dm:%2.0ds) \n",ct,hrs,min,sec);
                }else{
                    printf("  Calculating ETC...\n");
                }
            }
        }
        // Set correct parameters for beam to ray trace
        //
        TxPos = td->TxPositions[pulseIndex] ;
        RxPos = td->RxPositions[pulseIndex] ;
        
        // Set up the kernel arguments
        //
        
        CL_CHECK(clSetKernelArg(kernel, 0,   sizeof(int), &nAzBeam));
        CL_CHECK(clSetKernelArg(kernel, 1,   sizeof(int), &nElBeam));
        CL_CHECK(clSetKernelArg(kernel, 2,   sizeof(SPVector), &RxPos));
        CL_CHECK(clSetKernelArg(kernel, 3,   sizeof(SPVector), &TxPos));
        CL_CHECK(clSetKernelArg(kernel, 4,   sizeof(double), &dAz));
        CL_CHECK(clSetKernelArg(kernel, 5,   sizeof(double), &dEl));
        CL_CHECK(clSetKernelArg(kernel, 6,   sizeof(double), &raySolidAng));
        CL_CHECK(clSetKernelArg(kernel, 7,   sizeof(double), &TxPowPerRay));
        CL_CHECK(clSetKernelArg(kernel, 8,   sizeof(AABB), &SceneBoundingBox));
        CL_CHECK(clSetKernelArg(kernel, 9,   sizeof(double), &Aeff));
        CL_CHECK(clSetKernelArg(kernel, 10,  sizeof(int), &bounceToShow));
        CL_CHECK(clSetKernelArg(kernel, 11,  sizeof(cl_mem), &dTriangles));
        CL_CHECK(clSetKernelArg(kernel, 12,  sizeof(cl_mem), &dTextures));
        CL_CHECK(clSetKernelArg(kernel, 13,  sizeof(cl_mem), &dKdTree));
        CL_CHECK(clSetKernelArg(kernel, 14,  sizeof(cl_mem), &dtriListData));
        CL_CHECK(clSetKernelArg(kernel, 15,  sizeof(cl_mem), &dtriListPtrs));
        CL_CHECK(clSetKernelArg(kernel, 16,  sizeof(cl_mem), &drnp));
        CL_CHECK(clSetKernelArg(kernel, 17,  sizeof(int), &pulseIndex));


        // allocate buffer for answers and make sure its zeroed (using calloc)
        //
        rnp = (rangeAndPower *)calloc(nAzBeam*nElBeam*MAXBOUNCES, sizeof(rangeAndPower));
        if (rnp == NULL){
            printf("Error : failed to malloc rangeAndPower array in thread %d (size %ld)\n"
                   ,tid, sizeof(rangeAndPower)*nAzBeam*nElBeam*MAXBOUNCES);
            exit(-6);
        }
        
        CL_CHECK(clEnqueueWriteBuffer(commandQ,drnp,CL_TRUE,0,sizeof(rangeAndPower)*nAzBeam*nElBeam*MAXBOUNCES,rnp,0,NULL,NULL));
        CL_CHECK(clEnqueueNDRangeKernel(commandQ, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL));

        // Read results from device
        //
        CL_CHECK(clEnqueueReadBuffer(commandQ,drnp, CL_TRUE,0,sizeof(rangeAndPower)*nAzBeam*nElBeam*MAXBOUNCES,rnp,0,NULL,NULL));
        
        int cnt=0;
        for (int i=0; i<nAzBeam*nElBeam*MAXBOUNCES; i++){
            if (rnp[i].power != 0 && rnp[i].range !=0) cnt++ ;
        }
        if(cnt > 0){
            nrnpItems = cnt;
        }else{
            printf("ERROR: No intersections on device %d\n",tid);
            exit(-1);
        }
        
// DEBUG
//        nrnpItems = 1;
// DEBUG
        rnpData_t * rnpData = (rnpData_t *)malloc(sizeof(rnpData_t)*nrnpItems) ;
        if(rnpData == NULL){
            printf("Error : malloc fail for rnpData on device %d. request data for %d items\n",tid,nrnpItems) ;
            exit(-15);
        }
        
        VECT_MINUS(td->RxPositions[pulseIndex], aimdir) ;
        derampRange = VECT_MAG(aimdir);
        
        cnt = 0;
        for (int i=0; i<nAzBeam*nElBeam*MAXBOUNCES; i++){
            if (rnp[i].power != 0 && rnp[i].range !=0){
                rnpData[cnt].power = rnp[i].power ;
                rnpData[cnt].range = rnp[i].range ;
                rnpData[cnt].rdiff = -1*(derampRange - rnp[i].range);
                cnt++ ;
            }
        }
        
        ////////// DEBUG CODE to inject a single point target
        
//        float xpos = -5 ;
//        float ypos = 0 ;
//        SPVector tp ; // Target Point in metres from scene centre in E,N,Alt system
//        SPVector rvect1, rvect2 ;
//        SPVector R,sensor_pos, sc_pos, R_,KDP_,JDP_,IDP_;
//        SPVector rng_dir, azi_dir ;
//        sensor_pos = td->TxPositions[td->nPulses/2] ;
//        VECT_CREATE(0, 0, 0, sc_pos);
//        VECT_SUB(sensor_pos, sc_pos, R);
//        VECT_UNIT(R, R_);
//        VECT_CREATE(0, 0, 1, KDP_);     // == zbar
//        VECT_PROJ(R_, KDP_, JDP_);
//        VECT_CROSS(JDP_, KDP_, IDP_);
//        VECT_SCMULT(JDP_, -1*ypos, rng_dir);
//        VECT_SCMULT(IDP_,  xpos, azi_dir);
//        VECT_ADD(azi_dir, rng_dir, tp) ;
//        VECT_SUB(tp, TxPos, rvect1);
//        VECT_SUB(tp, RxPos, rvect2) ;
//        double r1,r2 ;
//        r1 = VECT_MAG(rvect1);
//        r2 = VECT_MAG(rvect2);
//        int phase_sign = -1.0;
//        rnpData[0].power = 1.0 ;
//        rnpData[0].range = (r1 + r2) / 2.0 ;
//        rnpData[0].rdiff = phase_sign * (derampRange - rnpData[0].range) ;
        
        ////////// DEBUG CODE to inject a single point target
        
#ifdef FASTPULSEBUILD
        
        
        im_create(&pulseLine, ITYPE_CMPL_FLOAT, nx, 1, 1.0, 1.0, &status);
        
        for (int i=0; i<nrnpItems; i++){
            
            double phse  = (-4.0 * SIPC_pi * rnpData[i].rdiff * td->oneOverLambda) ;
            targ.r = rnpData[i].power*cos(phse) ;
            targ.i = rnpData[i].power*sin(phse) ;
            
            rangeLabel = (rnpData[i].rdiff/sampSpacing) + (pulseLine.nx / 2) ;
            
            if (rangeLabel > NPOINTS/2 && rangeLabel < nx - NPOINTS) {
                packSinc(targ, pulseLine.data.cmpl_f, rnpData[i].rdiff, sampSpacing, pulseLine.nx, ikernel);
            }
        }
        
        // perform phase correction to account for deramped jitter in receiver timing
        //
        phasecorr = (((td->Fx0s[pulseIndex] - td->freq_centre) / td->FxSteps[pulseIndex])) * 2.0 * M_PI / pulseLine.nx;
        
        for(int x = 0; x < pulseLine.nx; x++) {
            pcorr.r = td->amp_sf0[pulseIndex] * cos(phasecorr * (x - pulseLine.nx/2)) ;
            pcorr.i = td->amp_sf0[pulseIndex] * sin(phasecorr * (x - pulseLine.nx/2)) ;
            
            CMPLX_MULT(pulseLine.data.cmpl_f[x], pcorr, tmp);
            pulseLine.data.cmpl_f[x] = tmp;
        }
        
        im_circshift(&pulseLine, -(pulseLine.nx/2), 0, &status);
        im_fftw(&pulseLine, FFT_X_ONLY+FWD+NOSCALE, &status);
        im_insert(&pulseLine, 0, pulseIndex, td->phd, &status) ;
        im_destroy(&pulseLine, &status) ;
        
#else
        
        for (int i=0; i<nrnpItems; i++){
            freqSampsToCentreStart       = (td->StartFrequency - td->Fx0s[pulseIndex]) / td->FxSteps[pulseIndex];
            ADCSampsToTargStart          = (2.0*rnpData[i].rdiff/SIPC_c)*td->ADRate;
            rnpData[i].samplingOffset    = freqSampsToCentreStart - ADCSampsToTargStart;
            rnpData[i].samplingOffsetInt = (int)(rnpData[i].samplingOffset);
            rnpData[i].indexOffset       = rnpData[i].samplingOffsetInt - rnpData[i].samplingOffset - (nx / 2) ;
        }
        
        int nbeamThreads = (int)(sysconf(_SC_NPROCESSORS_ONLN)/td->nThreads);
        int nbeamThreadsRem = nx % nbeamThreads ;
        int loops = (nbeamThreadsRem == 0) ? 1 : 2 ;
        int nSamp;
        int rc ;
        int thd ;
        
        for(int iloop = 0; iloop<loops; iloop++){
            if (iloop == 1) {
                nbeamThreads = nbeamThreadsRem ;
            }
            pthread_t *beamThreads;
            beamThreads = (pthread_t *)malloc(sizeof(pthread_t)*nbeamThreads) ;
            if (beamThreads == NULL) {
                printf("Error : Failed to malloc %d threads\n",nbeamThreads);
                exit(-1);
            }
        
            threadDataBF *threadDataArray ;
            threadDataArray = (threadDataBF *)malloc(sizeof(threadDataBF)*nbeamThreads);
            if (threadDataArray == NULL) {
                printf("Error : Failed to malloc %d threadDataArrays\n",nbeamThreads);
                exit(-1);
            }
        
            pthread_attr_t attr;
            pthread_attr_init(&attr);
            pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
        
            if(iloop==1){
                nSamp = 1 ;
            }else{
                nSamp = (int)(nx / nbeamThreads) ;
            }
           
            for (thd=0; thd<nbeamThreads; thd++) {
            
                threadDataArray[thd].A          = A ;
                threadDataArray[thd].B          = B ;
                threadDataArray[thd].nrnpItems  = nrnpItems ;
                threadDataArray[thd].nSamp      = nSamp ;
                threadDataArray[thd].nx         = (int)nx ;
                threadDataArray[thd].phd        = td->phd ;
                threadDataArray[thd].pulseIndex = pulseIndex ;
                threadDataArray[thd].rnpData    = rnpData ;
                if(iloop ==1){
                    threadDataArray[thd].startSamp  = (int)nx - nbeamThreads + thd ;
                }else{
                    threadDataArray[thd].startSamp  = thd*nSamp ;
                }
            
                // Create thread data for each device
                //
                rc = pthread_create(&beamThreads[thd], NULL, beamForm, (void *) &threadDataArray[thd]) ;
                if (rc){
                    printf("ERROR; return code from pthread_create() is %d\n", rc);
                    exit(-1);
                }
            }
        
            pthread_attr_destroy(&attr);
            void * threadStatus ;
            for (thd=0; thd<nbeamThreads; thd++) {
                rc = pthread_join(beamThreads[thd], &threadStatus);
                if (rc) {
                    printf("ERROR; return code from pthread_join() is %d\n", rc);
                    exit(-1);
                }
            }
            
            
            
            free(beamThreads);
            free(threadDataArray);
        }
        free(rnpData) ;
        
#endif
        
        free(rnp);
        
    } // end of pulse loop
    
    // Clear down OpenCL allocations
    //

    clReleaseCommandQueue(commandQ);
    clReleaseContext(context);
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseMemObject(dTriangles);
    clReleaseMemObject(dTextures);
    clReleaseMemObject(dKdTree);
    clReleaseMemObject(dtriListData);
    clReleaseMemObject(dtriListPtrs);
    clReleaseMemObject(drnp);
    
    // return to parent thread
    //
    pthread_exit(NULL) ;
}

void packSinc(SPCmplx point, SPCmplx *outData, double rdiff, double sampleSpacing, long long nxInData, double * ikernel)
{
    
    double diffFromSincCentre = rdiff / sampleSpacing + nxInData/2;
    int nearestSample = (int)(diffFromSincCentre+0.5);
    double nearestSincSample = nearestSample - diffFromSincCentre;
    int offset;
    int iidx = (int)(nearestSincSample*OVERSAMP);
    if (iidx < 0 ){
        iidx += OVERSAMP ;
        offset = 1;
    }else{
        offset = 0;
    }
    
    for(int x = 0; x < NPOINTS; x++)
    {
        outData[nearestSample+x-NPOINTS/2+offset].r += point.r * ikernel[iidx];
        outData[nearestSample+x-NPOINTS/2+offset].i += point.i * ikernel[iidx];
        iidx += OVERSAMP ;
    }

    return ;
}

double *
sinc_kernel(int oversample_factor, int num_of_points, double resolution, double sampleSpacing)
{
    int n = oversample_factor * num_of_points + 1;
    int x;
    double val;
    double * data;
    double k = M_PI / resolution ;
    double xpos ;
    
    if (oversample_factor == 0 || num_of_points == 0)
    {
        fprintf(stderr, "Either oversample_factor (= %d) or number of points (= %d) is zero in generate_sinc_kernel\n",
                oversample_factor, num_of_points);
        exit(803);
    }
    
    data = calloc(n, sizeof(double));
    if (!data)
    {
        fprintf(stderr, "Failed to calloc %s:%d\n", __FILE__, __LINE__);
        exit(802);
    }
    
    for(x = 0; x < n; x++)
    {
        xpos = (x - n/2) * (sampleSpacing/oversample_factor) ;
        val  = k * xpos ;
        data[x] = (fabs(val) < 1.0e-4 ) ? 1.0 : sin(val)/val;
    }
    
    ham1dx(data, n);
    
    return(data);
}

static void
ham1dx(double * data, int nx)
{
    int x;
    
    double ped = 0.08;
    double a = 1.0 - ped;
    double val;
    
    for(x = 0; x < nx; x++)
    {
        val = M_PI * (-0.5 + (double) x / (double) (nx - 1));
        data[x] *= ped + a * cos(val) * cos(val);
    }
}

