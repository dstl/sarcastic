/***************************************************************************
 * 
 *           Module :  oclBuildShadowRays.cpp
 *          Program :  bircs
 *       Created by :  Darren Muff on  12/06/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *     Function to build shadow rays for bircs using OpenCL
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

#include "threadCore.hpp"
#include "oclBuildShadowRays.hpp"


void oclBuildShadowRays(cl_context          context,            // OpenCL context - alrready built
                        cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                        cl_kernel           kernel,             // OpenCl kernel for this routine to call
                        size_t              localWorkSize,      // Local workgroupsize to use for OpenCL Kernel
                        int                 nRays,              // The number of reflected rays being considered
                        SPVector            RxPos,              // The Receiver location in x,y,z
                        Ray                 *reflectedRays,     // Array of reflected rays - used for their origin as its the reflection point to Rx
                        Ray                 *shadowRays,        // Output - shadow rays to be tested for occlusion later
                        double              *ranges             // Output - the range to the receiver for each shadowray. Calculated here to save calcs later
)

{
    cl_mem dReflectedRays ;
    cl_mem dShadowRays ;
    cl_mem dRanges ;
    
    size_t globalWorkSize ;
    
    // Make sure globalWorkSize is a multiple of localWorkSize
    //
    globalWorkSize = nRays ;
    while (globalWorkSize % localWorkSize)globalWorkSize++;
    
    // Create OpenCl device buffers to hold data for this ray cast
    //
    dReflectedRays = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Ray)*nRays,    NULL, &_err));
    dShadowRays    = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Ray)*nRays,    NULL, &_err)); // Should this be CL_MEM_READ_WRITE
    dRanges        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(double)*nRays, NULL, &_err)); // Should this be CL_MEM_READ_WRITE
    
    
    // Add write buffer commands to the command queue ready for execution
    //
    CL_CHECK(clEnqueueWriteBuffer(Q, dReflectedRays, CL_TRUE, 0, sizeof(Ray)*nRays, reflectedRays, 0, NULL, NULL));
    
    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(kernel, 0,   sizeof(int),      &nRays));
    CL_CHECK(clSetKernelArg(kernel, 1,   sizeof(SPVector), &RxPos));
    CL_CHECK(clSetKernelArg(kernel, 2,   sizeof(cl_mem),   &dReflectedRays));
    CL_CHECK(clSetKernelArg(kernel, 3,   sizeof(cl_mem),   &dShadowRays));
    CL_CHECK(clSetKernelArg(kernel, 4,   sizeof(cl_mem),   &dRanges));
    
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(Q, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL));
    
    // Read results from device
    //
    CL_CHECK(clEnqueueReadBuffer(Q,dShadowRays, CL_TRUE,0,sizeof(Ray)*nRays,   shadowRays, 0,NULL,NULL));
    CL_CHECK(clEnqueueReadBuffer(Q,dRanges,     CL_TRUE,0,sizeof(double)*nRays,ranges,     0,NULL,NULL));
    
    // Clear down OpenCL allocations
    //
    clReleaseMemObject(dReflectedRays);
    clReleaseMemObject(dShadowRays);
    clReleaseMemObject(dRanges);
    return ;
}
