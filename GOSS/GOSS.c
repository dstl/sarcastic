/***************************************************************************
 *
 *       Module:    GOSS.c
 *      Program:    GOSS
 *   Created by:    Darren on 27/07/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *   <ENTER DESCRIPTION HERE>
 One set up the algorithm splits pulses across the number of available devices
 each device then runs its own pthread to process its own block of pulses.
 The pulses are interleaved across the threads (pulse n -> thread n%ndevs)
 This is so that the file can be written in parallel as each thread completes
 a loop
 Each thread loops through its pulses.
 Each pulse is then processed using OpenCL to ray trace all the rays within the 
 pulse's beam.
 Each OpenCL call returns an array of intersections (azBeam x elBeam x MAXBOUNCE
 storing (doubl)range & (double)power)
 These are then combined in the thread to form a single pulse return
 Each thread then writes its pulse to file in thread order (blocking write top the
 other threads as it  does)
 *
 *
 *   CLASSIFICATION        :  <PENDING>
 *   Date of CLASSN        :  14/03/2013
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
#include "readKdTree.h"
#include "BuildRopesAndBoxes.h"
#include "ecef2SceneCoords.h"


double TxPowerPerRay(int xRays, int yRays, double xBeamUsed, double yBeamUsed, double * raySolidAngle, double *effectiveArea);


int main (int argc, char **argv){
    
    CPHDHeader  hdr ;
    AABB        SceneBoundingBox;
    OCLPlatform platform;
    
    SPVector SRP, unitBeamAz, unitBeamEl, rVect, zHat, boxPts[8] ;
    double centreRange, maxEl, maxAz, El, Az, dAz, dEl, lambda ;
    char *KdTreeFile, *inCPHDFile, *outCPHDFile ;
    int startPulse, nPulses, bounceToShow, nTriangles, nTextures, nLeaves, nTreeNodes ;
    int useGPU, nAzBeam, nElBeam, nVec, Ropes[6], dev, rc ;
    SPImage cphd ;
    int debug, debugX=0, debugY=0;

    cl_int      err;
    cl_uint     ndevs;
    cl_ulong    memSizeTmp,devMemSize;
    
    Triangle  * Triangles     = NULL ;
    Texture   * textures      = NULL ;
    int       **triangleLists = NULL ;
    KdData    * KdTree        = NULL ;
    TriCoords * tricos        = NULL ;
    SPVector  * RxPos         = NULL ;
    SPVector  * TxPos         = NULL ;
    double    * Fx0s          = NULL ;
    double    * FxSteps       = NULL ;

    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"GOSS", argc, (char **)argv);
	CHECK_STATUS_NON_PTR(status);
    
    KdTreeFile = (char *)malloc(sizeof(char)*256);
    inCPHDFile = (char *)malloc(sizeof(char)*256);
    outCPHDFile = (char *)malloc(sizeof(char)*256);
    
    banner() ;
    
    getUserInput(&inCPHDFile, &KdTreeFile, &outCPHDFile, &startPulse, &nPulses, &bounceToShow, &nAzBeam, &nElBeam, &useGPU, &debug, &debugX, &debugY,&status) ;
    
    // Start timing after user input
    //
    Timer runTimer ;
    startTimer(&runTimer, &status) ;
    
    // Read in KdTree
    //
    readKdTree(KdTreeFile,          // Filename for KdTree
               &SceneBoundingBox,   // Bounding box for scene
               &nTriangles,         // Number of triangles in scene
               &Triangles,          // Array containing triangles
               &nTextures,          // number if textures in scene
               &textures,           // array containing textures
               &nLeaves,            // number of leaf nodes in KdTree
               &triangleLists,      // Array of int arrays containing triangles for each leaf
               &nTreeNodes,         // Number of nodes in KdTree
               &KdTree,             // KdTree returned.
               &tricos);            // Triangle coordinates
    
    
    // Build Ropes and AABBs for faster stackless traversal
    //
    for(int i=0; i<6; i++) Ropes[i] = NILROPE;
    BuildRopesAndBoxes(KdTree, Ropes, SceneBoundingBox, KdTree);
    printf("Scene Summary:\n");
    printf("  Triangles : %d\n",nTriangles);
    printf("  Textures  : %d\n",nTextures);
    printf("  Leaves    : %d\n",nLeaves);
    printf("  Nodes     : %d\n",nTreeNodes);
    printf("  Bound Box : %2.2f,%2.2f,%2.2f - %2.2f,%2.2f,%2.2f\n",
           SceneBoundingBox.AA.x,SceneBoundingBox.AA.y,SceneBoundingBox.AA.z,
           SceneBoundingBox.BB.x,SceneBoundingBox.BB.y,SceneBoundingBox.BB.z);
    
    // cannot pass a pointer to a pointer in OpenCL so we need to sort out
    // triangleLists to be two lists:
    //      triListData : array containing all indices. size=nTriIndices
    //      triPtrs     : array containing index in triIndices that matches node. - size=nLeaves
    //
    int * triPtrs = (int *)malloc(sizeof(int)*nLeaves);
    int nt=0;
    for (int n=0; n<nLeaves; n++)nt += triangleLists[n][0];
    int triListDataSize = sizeof(int)*(nt+nLeaves);
    int * triListData = (int *)malloc(triListDataSize);
    int trisInNode;
    int indCnt=0;
    for (int n=0; n<nLeaves; n++){
        trisInNode = triangleLists[n][0];
        triListData[indCnt] = trisInNode;
        triPtrs[n] = indCnt++;
        for(int m=0; m<trisInNode; m++)triListData[indCnt++] = triangleLists[n][m+1];
    }
    free(triangleLists);
    
    // Check the cphdfile
    //
    readCPHDHeader(inCPHDFile, &hdr, &status);
    nVec = hdr.num_azi ;
    
    // print out info about the collection
    //
    printCPHDCollectionInfo(&hdr, NULL, &status) ;
    
    // Rotate Rx and Tx Coords to be relative to scene centre
    //
    if((RxPos = (SPVector *)malloc(sizeof(SPVector)*nPulses))==NULL){
        printf("Error : Malloc failed allocating %ld bytes for RxPos\n",sizeof(SPVector)*nVec);
        exit(-1);
    }
    if((TxPos = (SPVector *)malloc(sizeof(SPVector)*nPulses))==NULL){
        printf("Error : Malloc failed allocating %ld bytes for TxPos\n",sizeof(SPVector)*nVec);
        exit(-1);
    }
    if((Fx0s = (double *)malloc(sizeof(double)*nPulses))==NULL){
        printf("Error : Malloc failed allocating %ld bytes for Fx0s\n",sizeof(double)*nVec);
        exit(-1);
    }
    if((FxSteps = (double *)malloc(sizeof(double)*nPulses))==NULL){
        printf("Error : Malloc failed allocating %ld bytes for FxSteps\n",sizeof(double)*nVec);
        exit(-1);
    }
    
    
    for (int p = 0; p < nPulses; p++){
        RxPos[p] = hdr.pulses[p+startPulse].sat_ps_rx ;
        TxPos[p] = hdr.pulses[p+startPulse].sat_ps_rx ;
        
        // While we are looping through pulses also load up the Fx0 and FxStepSize values
        //
        Fx0s[p]    = hdr.pulses[p+startPulse].fx0 ;
        FxSteps[p] = hdr.pulses[p+startPulse].fx_step_size ;
    }
    
    SRP = hdr.pulses[startPulse+(nPulses/2)].srp ;

    ecef2SceneCoords(nPulses, RxPos, SRP);
    ecef2SceneCoords(nPulses, TxPos, SRP);

    // Calculate azbeamwidth from scene
    // Calculate elbeamwidth from scene
    // use largest so that beam sampling is uniform
    // multiply by sceneBeamMargin
    // set to be beamwidth
    
    centreRange = VECT_MAG(TxPos[nPulses/2]);
    VECT_MINUS( TxPos[nPulses/2], rVect ) ;
    VECT_CREATE(0, 0, 1., zHat) ;
    VECT_CROSS(rVect, zHat, unitBeamAz);
    VECT_NORM(unitBeamAz, unitBeamAz) ;
    VECT_CROSS(unitBeamAz, rVect, unitBeamEl) ;
    VECT_NORM(unitBeamEl, unitBeamEl) ;
    
    maxEl = maxAz = 0.0 ;
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.AA.y, SceneBoundingBox.AA.z, boxPts[0]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.BB.y, SceneBoundingBox.AA.z, boxPts[1]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.BB.y, SceneBoundingBox.AA.z, boxPts[2]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.AA.y, SceneBoundingBox.AA.z, boxPts[3]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.AA.y, SceneBoundingBox.BB.z, boxPts[4]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.BB.y, SceneBoundingBox.BB.z, boxPts[5]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.BB.y, SceneBoundingBox.BB.z, boxPts[6]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.AA.y, SceneBoundingBox.BB.z, boxPts[7]);
    
    for( int k=0; k<8; k++){
        El = fabs(VECT_DOT(boxPts[k], unitBeamEl)) ;
        Az = fabs(VECT_DOT(boxPts[k], unitBeamAz)) ;
        maxEl = ( maxEl < El ) ? El : maxEl ;
        maxAz = ( maxAz < Az ) ? Az : maxAz ;
    }
    dAz = 2.0 * maxAz / centreRange ;
    dEl = 2.0 * maxEl / centreRange ;
    
    printf("  Beamwidth (Az,El)         : %f deg x %f deg (%f x %f metres)\n",
           GeoConsts_RADTODEG*dAz,GeoConsts_RADTODEG*dEl,centreRange*dAz,centreRange*dEl);
    dAz = dAz / nAzBeam;
    dEl = dEl / nElBeam;
    lambda = SIPC_c / hdr.freq_centre ;
    printf("  Ray density               : %f rays per wavelength cell\n",
           (lambda / (centreRange * dAz))*(lambda/ (centreRange *dEl)));
    
    double raySolidAng,TxPowPerRay, Aeff ;
    TxPowPerRay = TxPowerPerRay(nAzBeam, nElBeam, dAz, dEl,&raySolidAng, &Aeff);

    // Initialise OpenCL and load the relevent information into the
    // platform structure. OCL tasks will use this later
    //
    
    platform.clSelectedPlatformID = NULL;
    err = oclGetPlatformID (&platform.clSelectedPlatformID, &status);
    if(err != CL_SUCCESS){
        printf("Error: Failed to find a suitable OpenCL Launch platform\n");
        exit(-1);
    }
    oclPrintPlatform(platform);
    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform.clSelectedPlatformID,
        0
    };
    platform.props = props ;
        
    // Find number of devices on this platform
    //
    if (useGPU){
#ifdef NVIDIA
        err = oclGetNvidiaDevices(platform.clSelectedPlatformID, GPUCAPABILITY_MAJOR, GPUCAPABILITY_MINOR, &platform.device_ids, &ndevs, &status);
#else
        err = oclGetGPUDevices(platform.clSelectedPlatformID, 1, &platform.device_ids, &ndevs, &status) ;
#endif
        if(err == CL_SUCCESS){
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
        err = oclGetCPUDevices(platform.clSelectedPlatformID, &platform.device_ids, &ndevs, &status);
        if(err == CL_SUCCESS){
            printf("CPU DEVICES                 : %d\n",ndevs);
        }
    }

    
    // Find out the maximum amount of memory that all OpenCL devices have
    // For this task
    //
    devMemSize = 100e9;    // big number 100 GB
    for (int dev=0; dev<ndevs; dev++){
        err = clGetDeviceInfo(platform.device_ids[dev], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &memSizeTmp, NULL);
        if(err != CL_SUCCESS){
            printf("Error [%d] : couldn't obtain device info\n",err);
        }
        devMemSize = (memSizeTmp < devMemSize) ? memSizeTmp : devMemSize;
    }
    if (bounceToShow != 0) ndevs = 1;
    platform.nDevs = ndevs ;

    printf("DEVICES USED                : %d\n",ndevs);
    
    pthread_t *threads;
    threads = (pthread_t *)malloc(sizeof(pthread_t)*ndevs) ;
    if (threads == NULL) {
        printf("Error : Failed to malloc %d threads\n",ndevs);
        exit(-1);
    }
    
    threadData *threadDataArray ;
    threadDataArray = (threadData *)malloc(sizeof(threadData)*ndevs);
    if (threadDataArray == NULL) {
        printf("Error : Failed to malloc %d threadDataArrays\n",ndevs);
        exit(-1);
    }
    
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    im_init(&cphd, &status) ;
    load_cphd(&cphd, &hdr, 0, hdr.num_azi, &status);
    
    // wiping pulse data...
    //
    for (int i = 0; i<cphd.nx*cphd.ny; i++){
        cphd.data.cmpl_f[i].i = cphd.data.cmpl_f[i].r = 0.0f;
    }
    
    // Calculate how many pulses can be processed on a device at a time
    //
    int pulsesPerDevice;
   
    while (nPulses % ndevs != 0) nPulses-- ;
    pulsesPerDevice = nPulses / ndevs ;

    for (dev=0; dev<ndevs; dev++) {
        
        threadDataArray[dev].devIndex               = dev ;
        threadDataArray[dev].nThreads               = ndevs ;
        threadDataArray[dev].platform               = platform ;
        threadDataArray[dev].nTriangles             = nTriangles ;
        threadDataArray[dev].Triangles              = Triangles;
        threadDataArray[dev].nTextures              = nTextures;
        threadDataArray[dev].Textures               = textures;
        threadDataArray[dev].SceneBoundingBox       = SceneBoundingBox;
        threadDataArray[dev].nLeaves                = nLeaves ;
        threadDataArray[dev].triPtrs                = triPtrs;
        threadDataArray[dev].triListDataSize        = triListDataSize ;
        threadDataArray[dev].triListData            = triListData ;
        threadDataArray[dev].nTreeNodes             = nTreeNodes ;
        threadDataArray[dev].KdTree                 = KdTree ;
        threadDataArray[dev].startPulse             = startPulse + dev*pulsesPerDevice ; 
        threadDataArray[dev].nPulses                = pulsesPerDevice ;
        threadDataArray[dev].nAzBeam                = nAzBeam ;
        threadDataArray[dev].nElBeam                = nElBeam ;
        threadDataArray[dev].dAz                    = dAz ;
        threadDataArray[dev].dEl                    = dEl ;
        threadDataArray[dev].Aeff                   = Aeff ;            // The effective area of the Receive Antenna
        threadDataArray[dev].TxPositions            = TxPos ;           // Pointer to beginning of TxPos data
        threadDataArray[dev].RxPositions            = RxPos ;           // Pointer to beginning of RxPos data
        threadDataArray[dev].Fx0s                   = Fx0s ;            // Pointer to beginning of Fx0s data
        threadDataArray[dev].FxSteps                = FxSteps;          // Pointer to beginning of FxSteps data
        threadDataArray[dev].raySolidAng            = raySolidAng;
        threadDataArray[dev].TxPowPerRay            = TxPowPerRay ;
        threadDataArray[dev].phd                    = &cphd;            // Pointer to beginning of cphd data
        threadDataArray[dev].chirpRate              = hdr.chirp_gamma;
        threadDataArray[dev].ADRate                 = hdr.clock_speed ;
        threadDataArray[dev].oneOverLambda          = hdr.freq_centre / SIPC_c ;
        threadDataArray[dev].StartFrequency         = hdr.freq_centre - (hdr.pulse_length *hdr.chirp_gamma/2) ;
        threadDataArray[dev].bounceToShow           = bounceToShow-1;
        threadDataArray[dev].status                 = status ;
        threadDataArray[dev].debug                  = debug ;
        threadDataArray[dev].debugX                 = debugX ;
        threadDataArray[dev].debugY                 = debugY ;
        
        if (bounceToShow)printf("\n+++++++++++++++++++++++++++++++++++++++\n");
        // Create thread data for each device
        //
        rc = pthread_create(&threads[dev], NULL, devPulseBlock, (void *) &threadDataArray[dev]) ;
        if (rc){
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }
    
    pthread_attr_destroy(&attr);
    void * threadStatus ;
    for (dev=0; dev<ndevs; dev++) {
        rc = pthread_join(threads[dev], &threadStatus);
        if (rc) {
            printf("ERROR; return code from pthread_join() is %d\n", rc);
            exit(-1);
        }
    }
    if (bounceToShow)printf("\n+++++++++++++++++++++++++++++++++++++++\n");

    FILE *fp;
    fp = fopen(outCPHDFile, "w") ;
    if (fp == NULL) {
        fprintf(stderr, "Failed to open file %s\n", outCPHDFile);
        perror("Error opening file");
        exit(888);
    }
    
    endTimer(&runTimer, &status) ;
    printf("Done  is %f seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
    printf("Writing CPHD File \"%s\"....",outCPHDFile);
    hdr.data_type = ITYPE_CMPL_FLOAT ;
    writeCPHD3Header( &hdr, fp, &status ) ;
    write_cphd3_nb_vectors(&hdr, 0, fp, &status) ;
    write_cphd3_wb_vectors(&hdr, &cphd, 0, fp, &status) ;
    printf("...Done\n");
    
    im_destroy(&cphd, &status);
    im_close_lib(&status);
    return 0;

}

