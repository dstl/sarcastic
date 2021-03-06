/***************************************************************************
 * 
 *           Module :  oclReflect.cpp
 *          Program :  bircs
 *       Created by :  Darren Muff on 12/06/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *     Function to build an array of reflection rays for bircs using OpenCL
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
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
 ***************************************************************************/
#include "oclReflect.hpp"
#include "threadCore.hpp"

void oclReflect(cl_context          context,            // OpenCL context - alrready built
                cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                cl_kernel           kernel,             // OpenCl kernel for this routine to call
                cl_mem              dTriangles,         // Device array of triangles
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
    CL_CHECK(clSetKernelArg(kernel, 1,   sizeof(cl_mem), &dRays));
    CL_CHECK(clSetKernelArg(kernel, 2,   sizeof(cl_mem), &dHits));
    CL_CHECK(clSetKernelArg(kernel, 3,   sizeof(cl_mem), &dReflected));
    CL_CHECK(clSetKernelArg(kernel, 4,   sizeof(int),    &nRays));
    
    
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
