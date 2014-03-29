/***************************************************************************
 *
 *       Module:    oclReflect.c
 *      Program:    GOSS
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Function to use opencl to reflect an array of rays
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  09/03/2014
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

void oclReflect(cl_context          context,            // OpenCL context - alrready built
                cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                cl_kernel           kernel,             // OpenCl kernel for this routine to call
                cl_mem              dTriangles,
                cl_mem              dTextures,
                int                 nTextures,
                int                 nRays,              // Number of rays to reflect
                size_t              localWorkSize,      // Local workgroupsize to use for OpenCL Kernel
                Ray                 *rays,              // Array of rays to consider
                Hit                 *hits,              // Array of hit points for each ray
                Ray                 *reflectedRays      // output array of reflected rays
                )

{
    cl_mem dRays ;
    cl_mem dHits ;
    cl_mem dReflected ;
    
    size_t globalWorkSize ;
    
    // Make sure globalWorkSize is a multiple of localWorkSize
    //
    globalWorkSize = nRays ;
    while (globalWorkSize % localWorkSize)globalWorkSize++;
    
    // Create OpenCl device buffers to hold data for this ray cast
    //
    dRays        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Ray)*nRays,                      NULL, &_err));
    dHits        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Hit)*nRays,                      NULL, &_err));
    dReflected   = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Ray)*nRays,                      NULL, &_err));
    
    // Add write buffer commands to the command queue ready for execution
    //
    CL_CHECK(clEnqueueWriteBuffer(Q, dRays,        CL_TRUE, 0, sizeof(Ray)*nRays,           rays,             0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dHits,        CL_TRUE, 0, sizeof(Hit)*nRays,           hits,             0, NULL, NULL));
    
    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(kernel, 0,   sizeof(cl_mem), &dTriangles));
    CL_CHECK(clSetKernelArg(kernel, 1,   sizeof(cl_mem), &dTextures));
    CL_CHECK(clSetKernelArg(kernel, 2,   sizeof(int),    &nTextures));
    CL_CHECK(clSetKernelArg(kernel, 3,   sizeof(cl_mem), &dRays));
    CL_CHECK(clSetKernelArg(kernel, 4,   sizeof(cl_mem), &dHits));
    CL_CHECK(clSetKernelArg(kernel, 5,   sizeof(cl_mem), &dReflected));
    CL_CHECK(clSetKernelArg(kernel, 6,   sizeof(int),    &nRays));

    
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(Q, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL));
    
    // Read results from device
    //
    CL_CHECK(clEnqueueReadBuffer(Q,dReflected, CL_TRUE,0,sizeof(Ray)*nRays,reflectedRays,0,NULL,NULL));
    
    // Clear down OpenCL allocations
    //
    clReleaseMemObject(dRays);
    clReleaseMemObject(dHits);
    clReleaseMemObject(dReflected);
    return ;
}