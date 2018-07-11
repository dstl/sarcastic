/***************************************************************************
 *
 *       Module:    oclSRPRanges.c
 *      Program:    tdpocl
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Function to use opencl to precalculate platform to SRP ranges
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  05/05/2014
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

#include <sarclib/sarclib.h>
#include "OpenCLUtils.h"

char * SRPRangesKernelCode = "                                                                                                  \n \
#pragma OPENCL EXTENSION cl_khr_fp64: enable                                                                                    \n \
                                                                                                                                \n \
    typedef struct { double x, y, z; } VectorH ;                                                                                \n \
                                                                                                                                \n \
__kernel void SRPRanges(int n, __global VectorH *TXs, __global VectorH *RXs, __global VectorH *SRPs, __global double * ranges){ \n \
                                                                                                                                \n \
    int pnt = get_global_id(0);                                                                                                 \n \
    double rng_to_srp,r1,r2 ;                                                                                                   \n \
    double3 txp,rxp,srp;                                                                                                        \n \
                                                                                                                                \n \
    if (pnt >= 0 && pnt < n ){                                                                                                  \n \
                                                                                                                                \n \
        txp.x = TXs[pnt].x  ;                                                                                                   \n \
        txp.y = TXs[pnt].y  ;                                                                                                   \n \
        txp.z = TXs[pnt].z  ;                                                                                                   \n \
        rxp.x = RXs[pnt].x  ;                                                                                                   \n \
        rxp.y = RXs[pnt].y  ;                                                                                                   \n \
        rxp.z = RXs[pnt].z  ;                                                                                                   \n \
        srp.x = SRPs[pnt].x ;                                                                                                   \n \
        srp.y = SRPs[pnt].y ;                                                                                                   \n \
        srp.z = SRPs[pnt].z ;                                                                                                   \n \
                                                                                                                                \n \
        r1 = distance( rxp, srp) ;                                                                                              \n \
        r2 = distance( txp, srp) ;                                                                                              \n \
        rng_to_srp = r1+r2;                                                                                                     \n \
                                                                                                                                \n \
        ranges[pnt] = rng_to_srp ;                                                                                              \n \
    }                                                                                                                           \n \
                                                                                                                                \n \
    return ;                                                                                                                    \n \
}                                                                                                                               \n \
" ;

//void oclDestroyTask(OCLTask *task){
//    
//    free(task->device_ids);
//    free(task->contexts);
//    free(task->programs);
//    free(task->kernels);
//    free(task->commands);
//    free((char *)task->kernelName);
//    
//    
//    return ;
//}

void oclSRPRanges(OCLPlatform platform,     // OpenCl Platform definition structure
                  int nPulses,              // Number of pulses to convert
                  SPVector *TXs,            // Array of transmit locations. One for each pulse
                  SPVector *RXs,            // Array of receive locations. One for each pulse
                  SPVector *SRP,            // Array of SRP locations. One for each pulse
                  double **ranges,          // Pointer to array of ranges to be filled. One for each pulse
                  SPStatus *status          // Status array
)
{
    
    CHECK_STATUS(status);

    OCLTask task ;
    oclCreateTask(&platform, platform.nDevs, "SRPRanges", SRPRangesKernelCode, &task, status);
    
    size_t localWorkSize ;
    size_t globalWorkSize ;

    // As this is a pretty quick routine and only gets called once then we will only use the first
    // OpenCL device. This means we dont have to worry about splitting the pulses across several devices
    // and ensuring that nPulses is an integer multiple of teh number of devices.
    //
    const int nDevs     = 1 ;
    const int devPulses = nPulses / nDevs ;
    
    
    cl_mem *dTXs, *dRXs, *dSRP, *dRanges ;
    
    dTXs    = (cl_mem *)sp_malloc(sizeof(cl_mem)*nDevs);
    dRXs    = (cl_mem *)sp_malloc(sizeof(cl_mem)*nDevs);
    dSRP    = (cl_mem *)sp_malloc(sizeof(cl_mem)*nDevs);
    dRanges = (cl_mem *)sp_malloc(sizeof(cl_mem)*nDevs);
    
    Timer timer ;
    startTimer(&timer, status);
    
    size_t lws_t ;
    size_t warpSize;
    
    // Calculate maximum local work size that this kernel can support on this device
    //
    CL_CHECK(clGetKernelWorkGroupInfo(task.kernels[0], task.device_ids[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(lws_t), &lws_t, NULL));
    
    // What do we think is the best performance size though for this kernel, based upon this device's properties
    //
    CL_CHECK(clGetKernelWorkGroupInfo(task.kernels[0], task.device_ids[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(warpSize), &warpSize, NULL));
    
    // For efficiency reduce the local workgroup size until it is a multiple of warpSize
    // As RandomRays requires 2D then do this by shrinking y dimension
    //
    int x ;
    x = (int)lws_t ;
    while ( x % warpSize ) x--;
    localWorkSize = x ;
    
    // Make sure globalWorkSize is a multiple of localWorkSize
    //
    globalWorkSize = devPulses ;
    while (globalWorkSize % localWorkSize)globalWorkSize++;
    
    
    dTXs[0]    = CL_CHECK_ERR(clCreateBuffer(task.contexts[0], CL_MEM_READ_ONLY,  sizeof(SPVector)*devPulses, NULL, &_err));
    dRXs[0]    = CL_CHECK_ERR(clCreateBuffer(task.contexts[0], CL_MEM_READ_ONLY,  sizeof(SPVector)*devPulses, NULL, &_err));
    dSRP[0]    = CL_CHECK_ERR(clCreateBuffer(task.contexts[0], CL_MEM_READ_ONLY,  sizeof(SPVector)*devPulses, NULL, &_err));
    dRanges[0] = CL_CHECK_ERR(clCreateBuffer(task.contexts[0], CL_MEM_READ_WRITE, sizeof(double)*devPulses,   NULL, &_err));
    
    CL_CHECK(clEnqueueWriteBuffer(task.commands[0], dTXs[0],   CL_TRUE, 0, sizeof(SPVector)*devPulses, &(TXs[0*devPulses]), 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(task.commands[0], dRXs[0],   CL_TRUE, 0, sizeof(SPVector)*devPulses, &(RXs[0*devPulses]), 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(task.commands[0], dSRP[0],   CL_TRUE, 0, sizeof(SPVector)*devPulses, &(SRP[0*devPulses]), 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(task.commands[0], dRanges[0],CL_TRUE, 0, sizeof(double)*devPulses, &(*ranges[0*devPulses]), 0, NULL, NULL));
    
    
    CL_CHECK(clSetKernelArg(task.kernels[0], 0,   sizeof(int),    &devPulses));
    CL_CHECK(clSetKernelArg(task.kernels[0], 1,   sizeof(cl_mem), &(dTXs[0])));
    CL_CHECK(clSetKernelArg(task.kernels[0], 2,   sizeof(cl_mem), &(dRXs[0])));
    CL_CHECK(clSetKernelArg(task.kernels[0], 3,   sizeof(cl_mem), &(dSRP[0])));
    CL_CHECK(clSetKernelArg(task.kernels[0], 4,   sizeof(cl_mem), &(dRanges[0])));
    
    CL_CHECK(clEnqueueNDRangeKernel(task.commands[0], task.kernels[0], 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL));
    
    
    CL_CHECK(clEnqueueReadBuffer(task.commands[0],dRanges[0], CL_TRUE,0,sizeof(double)*devPulses,&(*ranges[0*devPulses]),0,NULL,NULL));
    
    endTimer(&timer, status);
    if(status->debug >= 10)
    printf("TX/RX to SRP ranges calculated in %f seconds\n",timeElapsedInSeconds(&timer, status));
    
    clReleaseMemObject(dTXs[0]);
    clReleaseMemObject(dRXs[0]);
    clReleaseMemObject(dSRP[0]);
    clReleaseMemObject(dRanges[0]);
    free(dTXs);
    free(dRXs);
    free(dSRP);
    free(dRanges);
    
    oclDestroyTask(&task, status);
    
    return ;
}

