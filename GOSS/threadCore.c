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

typedef struct beamParams {
    int nAzBeam ;                // Number of azimuth slices in beam
    int nElBeam ;                // Number of elevation slices in beam
    SPVector RxPos ;             // Receiver position for this pulse
    SPVector TxPos ;             // Transmitter position for this pulse
    double dAz ;                 // Azimuth slice in radians
    double dEl ;                 // Elevation slice in radians
    double raySolidAng ;         // Solid angle of a single ray
    double TxPowPerRay ;         // Transmitter power per ray
    AABB SceneBoundingBox ;      // Scene bounding box
    double Aeff ;                // The effective area of the Receive Antenna
    int bounceToShow;            // Which bounce to output

} beamParams ;

typedef struct rndData_t {
    double range;
    double power;
    double samplingOffset;
    int samplingOffsetInt;
    int indexOffset ;
    double rdiff;
}rndData_t ;


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
    beamParams beam ;
    double freqSampsToCentreStart, ADCSampsToTargStart ;
    SPVector aimdir ;
    double derampRange, currentReal, currentImag, phse, power ;
    int nrnpItems;
    
    printf("ptr to threadData %p\n",threadArg);
    printf("ptr to threadData TxPositions : %p\n",td->TxPositions);
    
    beam.nAzBeam          = td->nAzBeam ;
    beam.nElBeam          = td->nElBeam ;
    beam.dAz              = td->dAz ;
    beam.dEl              = td->dEl ;
    beam.raySolidAng      = td->raySolidAng ;
    beam.TxPowPerRay      = td->TxPowPerRay ;
    beam.SceneBoundingBox = td->SceneBoundingBox ;
    beam.bounceToShow     = td->bounceToShow ;
    
    double A = 4.0*SIPC_pi*td->chirpRate / (SIPC_c * td->ADRate);
    double B = 4.0*SIPC_pi*td->oneOverLambda ;
    
    int tid ;

    tid   = td->devIndex ;
    devId = td->platform.device_ids[tid] ;
    
    context = clCreateContext(td->platform.props, 1, &devId, NULL, NULL, &err);
    if (!context){
        printf("[thread:%d], Error [%d] : Failed to create a compute context! \n",tid,err);
        exit(-1);
    }
    
    program = clCreateProgramWithSource(context, 1, (const char **) &kernelCodePath, NULL, &err);
    if (!program){
        printf("[thread:%d], Error [%d]: Failed to create compute program for device!\n",tid,err);
        exit(-2);
    }
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
    
    kernel = clCreateKernel(program, "rayTraceBeam", &err);
    if (!kernel || err != CL_SUCCESS){
        printf("[thread:%d], Error: Failed to build kernel \"%s\" ! : %d\n",tid,"rayTraceBeam",err);
        exit(-4);
    }
    commandQ = clCreateCommandQueue(context, devId, 0, &err);
    if (!commandQ){
        printf("[thread:%d], Error: Failed to create a command commands!\n",tid);
        exit(-5);
    }
    
    // Calculate global and local work sizes
    //
    int tx,ty ;
    size_t localWorkSize [2] ;
    size_t globalWorkSize[2] ;
    globalWorkSize[0] = beam.nAzBeam ;
    globalWorkSize[1] = beam.nElBeam ;
    best2DWorkSize(kernel, devId, globalWorkSize[0], globalWorkSize[1], &tx, &ty, &td->status) ;
    localWorkSize[0]  = tx ;
    localWorkSize[1]  = ty ;
    
    // Create device buffers
    //
    dTriangles = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Triangle)*td->nTriangles, NULL, &err);
    if (err != CL_SUCCESS) { printf("Error : failed to create buffer for device %d : [%d]\n",tid,err); exit(-6); }
    dTextures = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Texture)*td->nTextures, NULL, &err);
    if (err != CL_SUCCESS) { printf("Error : failed to create buffer for device %d : [%d]\n",tid,err); exit(-7); }
    dKdTree = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(KdData)*td->nTreeNodes, NULL, &err);
    if (err != CL_SUCCESS) { printf("Error : failed to create buffer for device %d : [%d]\n",tid,err); exit(-8); }
    dtriListData = clCreateBuffer(context, CL_MEM_READ_ONLY, td->triListDataSize, NULL, &err);
    if (err != CL_SUCCESS) { printf("Error : failed to create buffer for device %d : [%d]\n",tid,err); exit(-9); }
    dtriListPtrs = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int)*td->nLeaves, NULL, &err);
    if (err != CL_SUCCESS) { printf("Error : failed to create buffer for device %d : [%d]\n",tid,err); exit(-10); }
    drnp = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(rangeAndPower)*beam.nAzBeam*beam.nElBeam*MAXBOUNCES, NULL, &err);
    if (err != CL_SUCCESS) { printf("Error : failed to create buffer for device %d : [%d]\n",tid,err); exit(-11); }

    // Load the device buffers
    //
    err  = clEnqueueWriteBuffer(commandQ, dTriangles,   CL_TRUE, 0, sizeof(Triangle)*td->nTriangles, td->Triangles,    0, NULL, NULL);
    err |= clEnqueueWriteBuffer(commandQ, dTextures,    CL_TRUE, 0, sizeof(Texture)*td->nTextures,td->Textures,        0, NULL, NULL);
    err |= clEnqueueWriteBuffer(commandQ, dKdTree,      CL_TRUE, 0, sizeof(KdData)*td->nTreeNodes, td->KdTree,         0, NULL, NULL);
    err |= clEnqueueWriteBuffer(commandQ, dtriListData, CL_TRUE, 0, td->triListDataSize, td->triListData,              0, NULL, NULL);
    err |= clEnqueueWriteBuffer(commandQ, dtriListPtrs, CL_TRUE, 0, sizeof(int)*td->nLeaves, td->triPtrs,              0, NULL, NULL);
    if (err != CL_SUCCESS){
        printf("Error: Failed to set write buffer for device %d : [%d]\n",tid, err);
        exit(-12);
    }
    
    // Set up the kernel arguments
    //
    err  = clSetKernelArg(kernel, 0,  sizeof(beamParams), &beam);
    err |= clSetKernelArg(kernel, 1,  sizeof(int), &td->bounceToShow);
    err |= clSetKernelArg(kernel, 2,  sizeof(cl_mem), &dTriangles);
    err |= clSetKernelArg(kernel, 3,  sizeof(cl_mem), &dTextures);
    err |= clSetKernelArg(kernel, 4,  sizeof(cl_mem), &dKdTree);
    err |= clSetKernelArg(kernel, 5,  sizeof(cl_mem), &dtriListData);
    err |= clSetKernelArg(kernel, 6,  sizeof(cl_mem), &dtriListPtrs);
    err |= clSetKernelArg(kernel, 7,  sizeof(cl_mem), &drnp);
    if (err != CL_SUCCESS){
        printf("Error: Failed to set kernel arguments! %d\n", err);
        exit(1);
    }
    
    // **** loop  start here
    //
    for (int pulse=0; pulse<td->nPulses; pulse++){
        int pulseIndex = td->startPulse+pulse ;
        
        VECT_MINUS(td->RxPositions[pulseIndex], aimdir) ;
        derampRange = VECT_MAG(aimdir);

        // allocate buffer for answers and make sure its zeroed (using calloc)
        //
        rnp = (rangeAndPower *)calloc(beam.nAzBeam*beam.nElBeam*MAXBOUNCES, sizeof(rangeAndPower));
        if (rnp == NULL){
            printf("Error : failed to malloc rangeAndPower array in thread %d (size %ld)\n"
                   ,tid, sizeof(rangeAndPower)*beam.nAzBeam*beam.nElBeam*MAXBOUNCES);
            exit(-6);
        }
        
        err = clEnqueueWriteBuffer(commandQ,drnp,CL_TRUE,0,sizeof(rangeAndPower)*beam.nAzBeam*beam.nElBeam*MAXBOUNCES,rnp,0,NULL,NULL);
        if (err != CL_SUCCESS){
            printf("Error: Failed to set write buffer for device %d : [%d]\n",tid, err);
            exit(-13);
        }
        
        // Set correct parameters for beam to ray trace
        //
        beam.TxPos.x = td->TxPositions[pulseIndex].x ;
        beam.TxPos.y = td->TxPositions[pulseIndex].y ;
        beam.TxPos.z = td->TxPositions[pulseIndex].z ;

        printf("Tx position is %f,%f,%f\n",td->TxPositions[pulseIndex].x,td->TxPositions[pulseIndex].y,td->TxPositions[pulseIndex].z );
        printf("beam.TxPos is %f,%f,%f\n", beam.TxPos.x, beam.TxPos.y, beam.TxPos.z) ;
        beam.RxPos = td->RxPositions[pulseIndex] ;
                
        err = clEnqueueNDRangeKernel(commandQ, kernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
        if (err){
            printf("Error: Failed to execute kernel! %d \n", err);
            exit(EXIT_FAILURE);
        }
        
        // Read results from device
        //
        err = clEnqueueReadBuffer(commandQ,drnp, CL_TRUE,0,sizeof(rangeAndPower)*beam.nAzBeam*beam.nElBeam*MAXBOUNCES,rnp,0,NULL,NULL);
        if (err != CL_SUCCESS){
            printf("Error: Failed to read rnp array from device %d! %d\n",tid , err);
            exit(-14);
        }
        
        // Remove all blank data from rnp
        //
        int cnt=0;
        for (int i=0; i<beam.nAzBeam*beam.nElBeam*MAXBOUNCES; i++){
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
        cnt=0;
        for (int i=0; i<beam.nAzBeam*beam.nElBeam*MAXBOUNCES; i++){
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
