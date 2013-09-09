/***************************************************************************
 *
 *       Module:    threadCore.c
 *      Program:    GOSS
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *   <ENTER DESCRIPTION HERE>
    This function runs on a single thread and controls a single OpenCL device
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

static char *kernelCodePath = "/Users/darren/Development/GOSS/GOSS/GOSSKernelCode.cl" ;

//typedef struct beamParams {
//    int nAzBeam ;                // Number of azimuth slices in beam
//    int nElBeam ;                // Number of elevation slices in beam
//    SPVector RxPos ;             // Receiver position for this pulse
//    SPVector TxPos ;             // Transmitter position for this pulse
//    double dAz ;                 // Azimuth slice in radians
//    double dEl ;                 // Elevation slice in radians
//    double raySolidAng ;         // Solid angle of a single ray
//    double TxPowPerRay ;         // Transmitter power per ray
//    AABB SceneBoundingBox ;      // Scene bounding box
//    double Aeff ;                // The effective area of the Receive Antenna
//    int bounceToShow;            // Which bounce to output
//
//} beamParams ;

typedef struct rndData_t {
    double range;
    double power;
    double samplingOffset;
    int samplingOffsetInt;
    int indexOffset ;
    double rdiff;
}rndData_t ;

#define CL_CHECK(_expr)                                                             \
    do {                                                                            \
        cl_int _err = _expr;                                                        \
            if (_err == CL_SUCCESS)                                                 \
            break;                                                                  \
        fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);    \
        abort();                                                                    \
} while (0)

#define CL_CHECK_ERR(_expr)                                                         \
({                                                                                  \
    cl_int _err = CL_INVALID_VALUE;                                                 \
    typeof(_expr) _ret = _expr;                                                     \
    if (_err != CL_SUCCESS) {                                                       \
        fprintf(stderr, "OpenCL Error: '%s' returned %d!\n", #_expr, (int)_err);    \
        abort();                                                                    \
    }                                                                               \
    _ret;                                                                           \
})

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
//    beamParams beam ;
    double freqSampsToCentreStart, ADCSampsToTargStart ;
    SPVector aimdir ;
    double derampRange, currentReal, currentImag, phse, power ;
    int nrnpItems;
    
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
    
    double A = 4.0*SIPC_pi*td->chirpRate / (SIPC_c * td->ADRate);
    double B = 4.0*SIPC_pi*td->oneOverLambda ;
    
    int tid ;

    tid   = td->devIndex ;
    devId = td->platform.device_ids[tid] ;
    
    context = CL_CHECK_ERR(clCreateContext(td->platform.props, 1, &devId, NULL, NULL, &_err));
    program = CL_CHECK_ERR(clCreateProgramWithSource(context, 1, (const char **) &kernelCodePath, NULL, &_err));
    
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err != CL_SUCCESS){
        size_t len;
        char buffer[32768];
        printf("[thread:%d], Error: Failed to build program executable!\n",tid);
        clGetProgramBuildInfo(program, devId, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("err: %d. Buffer:\n",err);
        printf("%s\n", buffer);
        exit(-3);
    }
    
    kernel = CL_CHECK_ERR(clCreateKernel(program, "rayTraceBeam", &_err));
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
    int reportN = 100 ;
    
    // **** loop  start here
    //
    for (int pulse=0; pulse<td->nPulses; pulse++){
        int pulseIndex = (tid * td->nPulses) + pulse ;
        if ( td->devIndex == 0 ) {
            if( pulse % reportN == 0){
                double pCentDone = 100.0*pulse/td->nPulses ;
                double sexToGo = estimatedTimeRemaining(&threadTimer, pCentDone, &status);
                printf("Processing pulses %6d - %6d out of %6d [%2d%%]",  pulse*td->nThreads, (pulse+reportN)*td->nThreads, td->nPulses*td->nThreads,(int)pCentDone);
                if(pulse != 0 ){
                    printf("  ETC in %3.0f Seconds\n",sexToGo);
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
        
        // Remove all blank data from rnp
        //
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
        rndData_t * rnpData = (rndData_t *)malloc(sizeof(rndData_t)*nrnpItems) ;
        if(rnpData == NULL){
            printf("Error : malloc fail for rnpData on device %d. request data for %d items\n",tid,nrnpItems) ;
            exit(-15);
        }
        
        VECT_MINUS(td->RxPositions[pulseIndex], aimdir) ;
        derampRange = VECT_MAG(aimdir);
        
        cnt=0;
        for (int i=0; i<nAzBeam*nElBeam*MAXBOUNCES; i++){
            if (rnp[i].power != 0 && rnp[i].range !=0){
                rnpData[cnt].power = rnp[i].power ;
                rnpData[cnt].range = rnp[i].range ;
                freqSampsToCentreStart = (td->StartFrequency - td->Fx0s[pulseIndex]) / td->FxSteps[pulseIndex];
                rnpData[cnt].rdiff = rnp[i].range - derampRange;
                ADCSampsToTargStart    = (2.0*rnpData[cnt].rdiff/SIPC_c)*td->ADRate;
                rnpData[cnt].samplingOffset = freqSampsToCentreStart - ADCSampsToTargStart;
                rnpData[cnt].samplingOffsetInt = (int)(rnpData[cnt].samplingOffset);
                rnpData[cnt].indexOffset = rnpData[cnt].samplingOffsetInt - rnpData[cnt].samplingOffset - (td->phd->nx / 2) ;
                cnt++;
            }
        }
        // Create the uncompressed (but deramped) range profile for this pulse, storing it
        // back in its original location
        //
        int samplingOffsetInt;
        for (int ksamp=0; ksamp<td->phd->nx; ksamp++) {
            for (int j=0; j<nrnpItems; j++ ){
                samplingOffsetInt = rnpData[j].samplingOffsetInt ;
                if ((ksamp+samplingOffsetInt >=0) && (ksamp+samplingOffsetInt < td->phd->nx)) {
                    currentReal = td->phd->data.cmpl_f[(pulseIndex)*td->phd->ny + (ksamp+rnpData[j].samplingOffsetInt)].r ;
                    currentImag = td->phd->data.cmpl_f[(pulseIndex)*td->phd->ny + (ksamp+rnpData[j].samplingOffsetInt)].i ;
                    
                    power = rnpData[j].power ;
                    phse  = rnpData[j].rdiff * ( A * (ksamp+rnpData[j].indexOffset) + B);
                    
                    td->phd->data.cmpl_f[(pulseIndex)*td->phd->ny + (ksamp+rnpData[j].samplingOffsetInt)].r = currentReal+power*cos(phse) ;
                    td->phd->data.cmpl_f[(pulseIndex)*td->phd->ny + (ksamp+rnpData[j].samplingOffsetInt)].i = currentImag+power*sin(phse) ;
                }
            }
        }
        
        free(rnpData) ;
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
