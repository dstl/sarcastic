/***************************************************************************
 *
 *       Module:    SARTrace.c
 *      Program:    SARTRACE
 *   Created by:    Darren on 27/07/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Uses the engine in Sarcastic just to ray trace intersection points 
 *      in a KdTree volume
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

#include "SARTrace.h"
#include "readKdTree.h"
#include "BuildRopesAndBoxes.h"
#include "ecef2SceneCoords.h"
#include "colourCodes.h"
#include "SartraceVersion.h"

int main (int argc, char **argv){
    
    CPHDHeader  hdr ;
    AABB        SceneBoundingBox;
    OCLPlatform platform;
    
    SPVector SRP, unitBeamAz, unitBeamEl, rVect, zHat, boxPts[8], interogPt ;
    double centreRange, maxEl, maxAz, minEl, minAz, dAz, dEl, lambda, maxBeamUsedAz, maxBeamUsedEl ;
    char *KdTreeFile, *inCPHDFile ;
    int startPulse, nPulses, nTriangles, nLeaves, nTreeNodes, nRaysX, nRaysY ;
    int useGPU, nVec, Ropes[6] ;
    
    cl_int      err;
    cl_uint     ndevs;
    cl_ulong    memSizeTmp,devMemSize;
    
    ATS       * accelTriangles= NULL ;
    int       **triangleLists = NULL ;
    KdData    * KdTree        = NULL ;
    Triangle  * triangles     = NULL ;
    SPVector  * RxPos         = NULL ;
    SPVector  * TxPos         = NULL ;

    nPulses = 1;
    useGPU  = 0;
    
    SPStatus status ;
    im_init_status(status, 0) ;
    im_init_lib(&status, (char *)"sartrace", argc, (char **)argv);
	CHECK_STATUS_NON_PTR(status);
    
    SARTracebanner();
    
    char *outDir;
    getSARTraceUserInput(&inCPHDFile, &KdTreeFile, &outDir, &nRaysX, &nRaysY, &status) ;

    // Start timing after user input
    //
    Timer runTimer ;
    startTimer(&runTimer, &status) ;
    
    // Read in KdTree
    //
    readKdTree(KdTreeFile,          // Filename for KdTree
               &SceneBoundingBox,   // Bounding box for scene
               &nTriangles,         // Number of triangles in scene
               &accelTriangles,     // Array containing accelerated triangle structures
               &triangles,          // Array containing triangles
               &nLeaves,            // number of leaf nodes in KdTree
               &triangleLists,      // Array of int arrays containing triangles for each leaf
               &nTreeNodes,         // Number of nodes in KdTree
               &KdTree);            // KdTree returned.
    
    
    // Build Ropes and AABBs for faster stackless traversal
    //
    for(int i=0; i<6; i++) Ropes[i] = NILROPE;
    BuildRopesAndBoxes(KdTree, Ropes, SceneBoundingBox, KdTree);
    
    // cannot pass a pointer to a pointer in OpenCL so we need to sort out
    // triangleLists to be two lists:
    //      triListData : array containing all indices. size=nTriIndices
    //      triPtrs     : array containing index in triIndices that matches node. - size=nLeaves
    //
    int * triPtrs = (int *)malloc(sizeof(int)*nLeaves);
    int nt=0; // Number of triangles
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
    // clear down triangle lists
    //
    for (int i=0; i<nLeaves; i++) { free(triangleLists[i]); }
    free(triangleLists);
    
    KdTreeStruct KDT ;
    KDT.KdTree           = KdTree ;
    KDT.nLeaves          = nLeaves ;
    KDT.nTreeNodes       = nTreeNodes ;
    KDT.nTriangles       = nTriangles ;
    KDT.triangleListData = triListData ;
    KDT.triangleListPtrs = triPtrs ;
    KDT.accelTriangles   = accelTriangles ;
    KDT.triListDataSize  = triListDataSize ;
    
    // Check the cphdfile
    //
    readCPHDHeader(inCPHDFile, &hdr, &status);
    nVec = hdr.num_azi ;
    startPulse = nVec / 2 ;
    
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
    
    for (int p = 0; p < nPulses; p++){
        RxPos[p] = hdr.nbdata[p+startPulse].sat_ps_rx ;
        TxPos[p] = hdr.nbdata[p+startPulse].sat_ps_rx ;
    }
    
    SRP = hdr.nbdata[startPulse+(nPulses/2)].srp ;

    ecef2SceneCoords(nPulses, RxPos, SRP);
    ecef2SceneCoords(nPulses, TxPos, SRP);
    ecef2SceneCoords(1,  &interogPt, SRP);

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
    
    maxEl = maxAz = minEl = minAz = 0.0 ;
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.AA.y, SceneBoundingBox.AA.z, boxPts[0]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.BB.y, SceneBoundingBox.AA.z, boxPts[1]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.BB.y, SceneBoundingBox.AA.z, boxPts[2]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.AA.y, SceneBoundingBox.AA.z, boxPts[3]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.AA.y, SceneBoundingBox.BB.z, boxPts[4]);
    VECT_CREATE(SceneBoundingBox.AA.x, SceneBoundingBox.BB.y, SceneBoundingBox.BB.z, boxPts[5]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.BB.y, SceneBoundingBox.BB.z, boxPts[6]);
    VECT_CREATE(SceneBoundingBox.BB.x, SceneBoundingBox.AA.y, SceneBoundingBox.BB.z, boxPts[7]);
    
    for( int k=0; k<8; k++){
        double El = VECT_DOT(boxPts[k], unitBeamEl) ;
        double Az = VECT_DOT(boxPts[k], unitBeamAz) ;
        maxEl = ( maxEl < El ) ? El : maxEl ;
        maxAz = ( maxAz < Az ) ? Az : maxAz ;
        minEl = ( minEl > El ) ? El : minEl ;
        minAz = ( minAz > Az ) ? Az : minAz ;
    }
    maxBeamUsedAz = (maxAz - minAz) / centreRange ;
    maxBeamUsedEl = (maxEl - minEl) / centreRange ;
    
    collectionGeom cGeom;
    collectionGeometry(&hdr, nPulses/2, hdr.grp, &cGeom, &status);
    printf("RayTracing beam (Az,El)     : %f deg x %f deg (%3.2f x %3.2f metres (ground plane))\n",
           GeoConsts_RADTODEG*maxBeamUsedAz,GeoConsts_RADTODEG*maxBeamUsedEl,centreRange*maxBeamUsedAz,centreRange*maxBeamUsedEl/sin(cGeom.grazingRad));
    dAz = maxBeamUsedAz / nRaysX;
    dEl = maxBeamUsedEl / nRaysY;
    lambda = SIPC_c / hdr.freq_centre ;
    printf("Ray density                 : %3.1f x %3.1f [ %3.1f x %3.1f ground plane] rays per metre\n",
           1.0/(dAz * centreRange), 1.0/(dEl * centreRange),
           1.0/(dAz * centreRange), 1.0/(dEl * centreRange / sin(cGeom.grazingRad)));
    printf("Ray seperation              : %4.3f x %4.3f [ %4.3f x %4.3f ground plane] metres between adjacent rays\n",
           (dAz * centreRange), (dEl * centreRange),
           (dAz * centreRange), (dEl * centreRange / sin(cGeom.grazingRad)));

    // Initialise OpenCL and load the relevent information into the
    // platform structure. OCL tasks will use this later
    //
    
    platform.clSelectedPlatformID = NULL;
    err = oclGetPlatformID (&platform, &status);
    if(err != CL_SUCCESS){
        printf("Error: Failed to find a suitable OpenCL Launch platform\n");
        exit(-1);
    }
    cl_context_properties props[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform.clSelectedPlatformID,
        0
    };
    platform.props = props ;
        
    // Find number of devices on this platform
    //
    if (useGPU){
        err = oclGetNamedGPUDevices(&platform, "NVIDIA", "",&platform.device_ids , &ndevs, &status);
        if(err == CL_SUCCESS){
            char cbuf[1024];
            cl_uint max_compute_units;
            printf("GPU DEVICES                 : %d\n",ndevs);
            for(int d=0; d<ndevs; d++){
                printf("DEVICE                      : %d\n",d);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_VENDOR, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE VENDOR             : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_NAME, sizeof(cbuf), &cbuf, NULL);
                printf("  DEVICE NAME               : %s\n",cbuf);
                clGetDeviceInfo(platform.device_ids[d], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(cl_uint), &max_compute_units, NULL);
                printf("  DEVICE MAX COMPUTE UNITS  : %d\n",max_compute_units);
            }
        }else{
            printf("No GPU devices available with compute capability: %d.%d\n", GPUCAPABILITY_MAJOR, GPUCAPABILITY_MINOR);
            exit(1);
        }
    }else{
        err = oclGetCPUDevices(&platform, &platform.device_ids, &ndevs, &status);
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
    ndevs = 1;
    platform.nDevs = ndevs ;
    
    threadData threadDataArray ;
  

    // Calculate how many pulses can be processed on a device at a time
    //
    
    threadDataArray.devIndex               = 0 ;
    threadDataArray.nThreads               = ndevs ;
    threadDataArray.platform               = platform ;
    threadDataArray.KDT                    = KDT;
    threadDataArray.SceneBoundingBox       = SceneBoundingBox;
    threadDataArray.startPulse             = startPulse;
    threadDataArray.nPulses                = 1 ;
    threadDataArray.nAzBeam                = nRaysX ;
    threadDataArray.nElBeam                = nRaysY ;
    threadDataArray.TxPositions            = TxPos ;           // Pointer to beginning of TxPos data
    threadDataArray.RxPositions            = RxPos ;           // Pointer to beginning of RxPos data
    threadDataArray.Fx0s                   = NULL ;            // Pointer to beginning of Fx0s data
    threadDataArray.FxSteps                = NULL;             // Pointer to beginning of FxSteps data
    threadDataArray.amp_sf0                = NULL ;            // Pointer to beginning of amp_sf0 data
    threadDataArray.gainRx                 = 1.0;
    threadDataArray.PowPerRay              = 1.0 ;
    threadDataArray.phd                    = NULL ;
    threadDataArray.chirpRate              = 0.0;
    threadDataArray.ADRate                 = 0.0;
    threadDataArray.pulseDuration          = 0.0;
    threadDataArray.oneOverLambda          = 0.0;
    threadDataArray.freq_centre            = 0.0;
    threadDataArray.StartFrequency         = 0.0;
    threadDataArray.bounceToShow           = 0 ;
    threadDataArray.status                 = status ;
    threadDataArray.interrogate            = 0 ;
    VECT_CREATE(0,0,0,threadDataArray.interogPt);
    threadDataArray.interogRad             = 0 ;
    threadDataArray.interogFP              = NULL ;
    threadDataArray.triangles              = triangles ;
    
    printf("\n+++++++++++++++++++++++++++++++++++++++\n");
    
    sartraceCore(threadDataArray, outDir);
    
    printf("\n+++++++++++++++++++++++++++++++++++++++\n");

    endTimer(&runTimer, &status) ;
    printf("Raycasting Completed in %f seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
    
    im_close_lib(&status);
    free ( RxPos ) ;
    free ( TxPos ) ;
    
    // Clear down the KdTree
    //
    free(accelTriangles);
    free(KdTree);
    free(triangles);

    return 0;

}

