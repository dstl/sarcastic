/***************************************************************************
 *
 *       Module:    oclKdTreeHits.cl
 *      Program:    SARCASTIC
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Function to use opencl to traverse kdTree
 *      Inputs are the KdTree and rays
 *      Outputs are hit points
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

void oclKdTreeHits(cl_context         context,            // OpenCL context - already instantiated
                   cl_command_queue   Q,                  // OpenCl command Q - already instatiated
                   cl_kernel          STkernel,           // stacklessTraverse kernel, already compiled and created from a program
                   int                nRays,              // Total number of rays to cast
                   size_t             localWorkSize,      // Work dimensions for this device
                   cl_mem             dTriangles,
                   cl_mem             dTextures,
                   cl_mem             dKdTree,
                   cl_mem             dtriListData,
                   cl_mem             dtriListPtrs,
                   AABB               SceneBoundingBox,   // Bounding box of scene - required for ray traversal optimisation
                   Ray *              rays,               // Array of rays (size nRays). Each ray will be cast through KdTree
                   Hit *              hits                // output array of hit locations
)
{
    cl_mem dHits ;
    cl_mem dRays ;
    size_t globalWorkSize ;
    
    // Make sure globalWorkSize is a multiple of localWorkSize
    //
    globalWorkSize = nRays ;
    while (globalWorkSize % localWorkSize)globalWorkSize++;
    
    // zero hits array so that we can see which rays hit something
    //
    memset(hits, 0, sizeof(Hit)*nRays);
    
    // Create OpenCl device buffers to hold data for this ray cast
    //
    dHits        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Hit)*nRays,           NULL, &_err));
    dRays        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(Ray)*nRays,           NULL, &_err));
    
    // Add write buffer commands to the command queue ready for execution
    //
    CL_CHECK(clEnqueueWriteBuffer(Q, dRays,        CL_TRUE, 0, sizeof(Ray)*nRays,           rays,             0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dHits,        CL_TRUE, 0, sizeof(Hit)*nRays,           hits,             0, NULL,NULL));
    
    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(STkernel, 0,  sizeof(cl_mem), &dKdTree));
    CL_CHECK(clSetKernelArg(STkernel, 1,  sizeof(cl_mem), &dtriListData));
    CL_CHECK(clSetKernelArg(STkernel, 2,  sizeof(cl_mem), &dtriListPtrs));
    CL_CHECK(clSetKernelArg(STkernel, 3,  sizeof(cl_mem), &dTriangles));
    CL_CHECK(clSetKernelArg(STkernel, 4,  sizeof(AABB),   &SceneBoundingBox));
    CL_CHECK(clSetKernelArg(STkernel, 5,  sizeof(int),    &nRays));
    CL_CHECK(clSetKernelArg(STkernel, 6,  sizeof(cl_mem), &dRays));
    CL_CHECK(clSetKernelArg(STkernel, 7,  sizeof(cl_mem), &dHits));
    
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(Q, STkernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL));
    
    // Read results from device
    //
    CL_CHECK(clEnqueueReadBuffer(Q,dRays, CL_TRUE,0,sizeof(Ray)*nRays,rays,0,NULL,NULL));
    CL_CHECK(clEnqueueReadBuffer(Q,dHits, CL_TRUE,0,sizeof(Hit)*nRays,hits,0,NULL,NULL));
    
    // Clear down OpenCL allocations
    //
    clReleaseMemObject(dHits);
    clReleaseMemObject(dRays);
    
    return ;
}