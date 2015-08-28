/***************************************************************************
 *
 *       Module:    threadCore.c
 *      Program:    SARTrace
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
#include "sarcastic.h"
#include <fftw3.h>

#define FASTPULSEBUILD
#define NPOINTS (32)
#define OVERSAMP (512)
#define NOINTERSECTION -1
static SPVector bouncecolours[MAXBOUNCES] = {
    {178.0, 255.0, 052.0},        // 0 = Bright Green
    {128.0, 128.0, 128.0},        // 1 = grey
    {224.0, 224.0, 224.0},        // 2 = White
    {176.0, 019.0, 035.0},        // 3 = Red
    {176.0, 094.0, 041.0},        // 4 = Orange
    {214.0, 186.0, 062.0},        // 5 = Yellow
    {166.0, 214.0, 054.0},        // 6 = Green
    {056.0, 125.0, 214.0},        // 7 = Blue
    {054.0, 039.0, 145.0},        // 8 = Indigo
    {121.0, 012.0, 165.0},        // 9 = Violet
};

void generate_gaussian_ray_distribution(double xStdev,double yStdev, SPVector txPos, SPVector GRP, int nx,int ny, Ray *rayArray);

    
void sartraceCore ( threadData td, char *outDir ) {
    
    cl_device_id devId ;
    cl_context context ;
    cl_command_queue commandQ ;

    int nAzBeam, nElBeam, bounceToShow, tid, norigrays ;
    int nbounce, nxRay, nyRay, nRays, reflectCount, iray, interrogate ;

    double gainRx, PowPerRay ;
    double bandwidth ;
    double k ;

    AABB SceneBoundingBox ;
    Ray *rayArray, *newRays, *reflectedRays;
    Hit *hitArray, *newHits;
    
    SPVector RxPos, TxPos, origin;
    
    gainRx           = td.gainRx ;
    PowPerRay        = td.PowPerRay ;
    bounceToShow     = td.bounceToShow ;
    nAzBeam          = td.nAzBeam ;
    nElBeam          = td.nElBeam ;
    SceneBoundingBox = td.SceneBoundingBox ;
    bandwidth        = td.chirpRate * td.pulseDuration ;
    tid              = td.devIndex ;
    devId            = td.platform.device_ids[tid] ;
    interrogate      = td.interrogate ;
    k                = 2 * SIPC_pi / (SIPC_c / td.freq_centre) ;
    
    VECT_CREATE(0, 0, 0, origin);
    
    // Create OpenCL context and command queue for this device
    //
    context  = CL_CHECK_ERR(clCreateContext(td.platform.props, 1, &devId, NULL, NULL, &_err));
    commandQ = CL_CHECK_ERR(clCreateCommandQueue(context, devId, 0, &_err));

    // We have the following OpenCL kernels in this thread:
    //  randomRays :  Generates a net of Gaussian distributed rays
    //  stackLessTraverse : Performs a stackless traversal of a KdTree to find the intersection points for rays
    //  reflect : Calculates the reflection ray for a net of rays hitting a surface
    //
    
    //  Build the kernels now and bail out if any fail to compile
    //
    cl_program randRaysPG,     stackTraversePG,  reflectPG ;
    cl_kernel  randRaysKL,     stackTraverseKL,  reflectKL ;
    size_t     randRaysLWS[2], stackTraverseLWS, reflectLWS ;

    static char *randRaysCode      = OCLKERNELSPATH"/randomRays.cl" ;
    static char *stackTraverseCode = OCLKERNELSPATH"/stacklessTraverse.cl" ;
    static char *reflectCode       = OCLKERNELSPATH"/reflectRays.cl" ;
    
    CL_CHECK(buildKernel(context, randRaysCode,      "randomRays",        devId, 2, &randRaysPG,      &randRaysKL,      randRaysLWS));
    CL_CHECK(buildKernel(context, stackTraverseCode, "stacklessTraverse", devId, 1, &stackTraversePG, &stackTraverseKL, &stackTraverseLWS));
    CL_CHECK(buildKernel(context, reflectCode,       "reflect",           devId, 1, &reflectPG,       &reflectKL,       &reflectLWS));

    // Allocate memory for items that are needed in all kernels
    //
    int                nTriangles       = td.KDT.nTriangles;         // number of triangles in array 'triangles'
    ATS *              accelTriangles   = td.KDT.accelTriangles;     // Array of triangles of size nTriangles
    int                nTreeNodes       = td.KDT.nTreeNodes;         // number of nodes in KdTree
    KdData *           KdTree           = td.KDT.KdTree;             // SAH - KdTree to optimise ray traversal through volume
    int                triListDataSize  = td.KDT.triListDataSize;    // size of trianglelist data
    int *              triangleListData = td.KDT.triangleListData;   // array of triangle indices into Triangles
    int                nLeaves          = td.KDT.nLeaves;            // number of leaves (nodes with triangles) in KdTree
    int *              triangleListPtrs = td.KDT.triangleListPtrs;   // array of pointers into triangleListData for each KdTree node

    cl_mem dTriangles,dKdTree, dtriListData,dtriListPtrs;
    dTriangles   = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(ATS)*nTriangles, NULL, &_err));
    dKdTree      = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(KdData)*nTreeNodes,   NULL, &_err));
    dtriListData = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  triListDataSize,             NULL, &_err));
    dtriListPtrs = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(int)*nLeaves,         NULL, &_err));
    
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dTriangles,   CL_TRUE, 0, sizeof(ATS)*nTriangles, accelTriangles,        0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dKdTree,      CL_TRUE, 0, sizeof(KdData)*nTreeNodes,   KdTree,           0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtriListData, CL_TRUE, 0, triListDataSize,             triangleListData, 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(commandQ, dtriListPtrs, CL_TRUE, 0, sizeof(int)*nLeaves,         triangleListPtrs, 0, NULL, NULL));
    
    int pulse = 0;
    int pulseIndex = (tid * td.nPulses) + pulse ;
    
    // Set correct parameters for beam to ray trace
    //
    TxPos = td.TxPositions[pulseIndex] ;
    RxPos = td.RxPositions[pulseIndex] ;
    
    // Generate a distribution of nAzbeam x nElbeam rays that originate form the TxPosition aiming at the origin. Use beamMax as the std deviation
    // for the distribution
    //
    rayArray = (Ray *)sp_malloc(nAzBeam*nElBeam * sizeof(Ray));
    
    printf("Random ray generation......");
    Timer runTimer ;
    SPStatus status;
    im_init_status(status, 0) ;

    startTimer(&runTimer, &status) ;
    oclRandomRays(context,commandQ,randRaysKL,nAzBeam,nElBeam,randRaysLWS,td.beamMaxAz,td.beamMaxEl,TxPos, origin, PowPerRay, rayArray);
    endTimer(&runTimer, &status);
    printf("Done (%8.2f ms = %8.2f rays/sec)\n", timeElapsedInMilliseconds(&runTimer, &status), nAzBeam*nElBeam/timeElapsedInMilliseconds(&runTimer, &status));
    
    nbounce = 0;
    nxRay   = nAzBeam ;
    nyRay   = nElBeam ;
    nRays   = nxRay*nyRay;
    
    while ( nbounce < MAXBOUNCES &&  nRays != 0){
        printf("Interflecting bounce %d.....",nbounce);
        startTimer(&runTimer, &status) ;
        norigrays = nRays ;
        
        // Malloc space for hits for this bounce
        //
        hitArray = (Hit *)sp_malloc(nRays * sizeof(Hit));
        
        // Cast the rays in rayArray through the KdTree using a stackless traversal technique
        // return the hit locations in hitArray
        //
        oclKdTreeHits(context, commandQ, stackTraverseKL, nRays, stackTraverseLWS, dTriangles, dKdTree, dtriListData, dtriListPtrs, SceneBoundingBox, rayArray, hitArray);
        
        // Sort out hits
        //
        reflectCount = 0 ;
        iray         = 0 ;
        
        // How many hits occured on this ray cast
        //
        for(int i=0; i<nRays; i++) {
            if ( hitArray[i].trinum != NOINTERSECTION ){
                reflectCount++ ;
            }
        }
        
        if( reflectCount == 0) {
            break ;
        }
        
        // shrink rayArray and hitArray to get rid of misses
        //
        newRays       = (Ray *)sp_malloc(reflectCount * sizeof(Ray)) ;
        newHits       = (Hit *)sp_malloc(reflectCount * sizeof(Hit)) ;
        reflectedRays = (Ray *)sp_malloc(reflectCount * sizeof(Ray)) ;
        
        for (int i=0; i<nRays; i++) {
            if ( hitArray[i].trinum != NOINTERSECTION ){
                newRays[iray] = rayArray[i] ;
                newHits[iray] = hitArray[i] ;
                iray++ ;
            }
        }
        free(rayArray) ;
        free(hitArray) ;
        rayArray = newRays ;
        hitArray = newHits ;
        nRays    = reflectCount ;
        
        // Build forward scattering rays ready for next turn round the loop
        //
        oclReflect(context, commandQ, reflectKL, dTriangles, nRays, reflectLWS, rayArray, hitArray, reflectedRays);
        endTimer(&runTimer, &status);
        
        int * trianglesHit = sp_calloc(td.KDT.nTriangles, sizeof(int));
        for (int i=0; i<nRays; i++){
            Hit h = hitArray[i] ;
            int trinum = h.trinum;
            trianglesHit[trinum] = 1;
        }

        unsigned int plyOutFlg;
        plyOutFlg = 3 ;
        if (plyOutFlg & 1 ) { // Point file out
            // Save origins of reflected rays to PLY file
            //
            char fname[256];
            sprintf(fname, "%s/bounce_%d.ply",outDir,nbounce);
            FILE *bounceFile = fopen(fname, "w");
            if (bounceFile == NULL) {
                printf("Error : could not open file %s for writing\n",fname);
                exit(1);
            }
            fprintf(bounceFile,"ply\n");
            fprintf(bounceFile,"format ascii 1.0\n");
            fprintf(bounceFile,"comment bounce %d - point intersections created by SARTrace\n",nbounce);
            fprintf(bounceFile,"element vertex %d\n",nRays);
            fprintf(bounceFile,"property float x\n");
            fprintf(bounceFile,"property float y\n");
            fprintf(bounceFile,"property float z\n");
            fprintf(bounceFile,"property uchar red\n");
            fprintf(bounceFile,"property uchar green\n");
            fprintf(bounceFile,"property uchar blue\n");
            fprintf(bounceFile,"end_header\n");
            
            for (int i=0; i<nRays; i++){
                //            printf("%f, %f, %f\n",reflectedRays[i].org.x,reflectedRays[i].org.y,reflectedRays[i].org.z);
                fprintf(bounceFile,"%f %f %f %d %d %d\n",reflectedRays[i].org.x,reflectedRays[i].org.y,reflectedRays[i].org.z,
                        (int)bouncecolours[nbounce].r,(int)bouncecolours[nbounce].g,(int)bouncecolours[nbounce].b);
            }
            fclose(bounceFile) ;

        }
        if (plyOutFlg & 2) { // Triangles that have been hit output
            char fname[256];
            sprintf(fname, "%s/trianglesHit_%d.ply",outDir,nbounce);
            FILE *fp = fopen(fname, "w");
            if (fp == NULL) {
                printf("Error : could not open file %s for writing\n",fname);
                exit(1);
            }
            // find unique vertices
            // newVerts is the new list of unique vertices and newVertCnt is the number of new vertices
            //
            int ntri=0;
            for (int i=0; i<td.KDT.nTriangles; i++) {if (trianglesHit[i] == 1) ntri++;}
            SPVector *v = (SPVector *)sp_malloc(sizeof(SPVector)*3*ntri);
            int ind = 0;
            for (int i=0; i<td.KDT.nTriangles; i++) {
                if (trianglesHit[i] == 1){
                    v[ind*3+0] = td.triangles[i].AA ;
                    v[ind*3+1] = td.triangles[i].BB ;
                    v[ind*3+2] = td.triangles[i].CC ;
                    ind++;
                }
            }
            SPVector test;
            SPVector *newVerts = (SPVector *)sp_malloc(sizeof(SPVector) * ntri * 3);
            int repeated;
            int newVertCnt = 0;
            for(int i=0; i<ntri*3; i++){
                test = v[i] ;
                repeated = 0;
                for (int j=0; j<newVertCnt; j++) {
                    if (newVerts[j].x == test.x && newVerts[j].y == test.y && newVerts[j].z == test.z ) repeated = 1;
                }
                if (!repeated) {
                    newVerts[newVertCnt++] = test ;
                }
            }
            
            // rebuild triangles using unique vertices
            //
            int t[ntri][3] ;
            ind = 0;
            for (int i=0; i<td.KDT.nTriangles; i++) {
                if(trianglesHit[i] == 1){
                    SPVector aa,bb,cc;
                    aa = td.triangles[i].AA ;
                    bb = td.triangles[i].BB ;
                    cc = td.triangles[i].CC ;
                    
                    for(int j=0; j<newVertCnt; j++){
                        if (newVerts[j].x == aa.x && newVerts[j].y == aa.y && newVerts[j].z == aa.z) t[ind][0] = j;
                        if (newVerts[j].x == bb.x && newVerts[j].y == bb.y && newVerts[j].z == bb.z) t[ind][1] = j;
                        if (newVerts[j].x == cc.x && newVerts[j].y == cc.y && newVerts[j].z == cc.z) t[ind][2] = j;
                    }
                    ind++;
                }
            }
            
            fprintf(fp,"ply\n");
            fprintf(fp,"format ascii 1.0\n");
            fprintf(fp,"comment bounce %d - Triangles hit created by SARTrace\n",nbounce);
            fprintf(fp,"element vertex %d\n",newVertCnt);
            fprintf(fp,"property float x\n");
            fprintf(fp,"property float y\n");
            fprintf(fp,"property float z\n");
            fprintf(fp,"element face %d\n",ntri);
            fprintf(fp,"property list uchar int vertex_index\n");
            fprintf(fp,"property uchar red\n");
            fprintf(fp,"property uchar green\n");
            fprintf(fp,"property uchar blue\n");
            fprintf(fp,"end_header\n");
            
            for (int i=0; i<newVertCnt; i++) {
                fprintf(fp,"%4.4f %4.4f %4.4f\n",newVerts[i].x,newVerts[i].y,newVerts[i].z);
            }
            for(int i=0; i<ntri; i++){
                fprintf(fp, "3 %d %d %d %d %d %d\n",t[i][0], t[i][1], t[i][2],
                        materialColours[accelTriangles[i].matInd][0],materialColours[accelTriangles[i].matInd][1],materialColours[accelTriangles[i].matInd][2]);
            }
            
            fclose(fp) ;
            free(v);
            free(newVerts);
            
        }
        printf("Done (%8.2f ms = %8.2f rays/sec)\n", timeElapsedInMilliseconds(&runTimer, &status), norigrays/timeElapsedInMilliseconds(&runTimer, &status));
        
        memcpy(rayArray, reflectedRays, sizeof(Ray)*reflectCount);
        free(reflectedRays);
        
        // Reset rayArrays and nRays for next time round loop
        //
        nRays = reflectCount ;
        
        nbounce++;
        
        free(hitArray);
    }
    
    free(rayArray) ;
    
    // Clear down OpenCL allocations
    //
    clReleaseKernel(randRaysKL);
    clReleaseKernel(stackTraverseKL);
    clReleaseKernel(reflectKL);
    clReleaseProgram(randRaysPG) ;
    clReleaseProgram(stackTraversePG) ;
    clReleaseProgram(reflectPG) ;
    clReleaseCommandQueue(commandQ);
    clReleaseContext(context);
    
    return ;
}
