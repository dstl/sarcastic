/***************************************************************************
 *
 *       Module:    oclReflectPower.c
 *      Program:    GOSS
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Function to use opencl to calculate the power reflect to the reciever
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



void oclReflectPower(cl_context          context,            // OpenCL context - alrready built
                     cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                     cl_kernel           kernel,             // OpenCl kernel for this routine to call
                     size_t              localWorkSize,      // Local workgroupsize to use for OpenCL Kernel
                     KdTreeStruct        KDT,                // Structure containing all KDTree info
                     Hit                 *hits,              // Array of hit locations to x-ref with triangles (and then Textures) for material props
                     SPVector            RxPos,              // Location of Receiver in x,y,z
                     double              GrxOverFourPi,      // Receiver antenna gain / 4Pi.
                     int                 nRays,              // The number of reflected rays being considered
                     Ray                 *rays,              // unit vector rays arriving at hitpoint
                     Ray                 *reflectedRays,     // Array of reflected rays - used for their origin as its the reflection point to Rx
                     Ray                 *shadowRays,        // Viewpoint unit vector. The vector direction to calculate power for (usually to receiver)
                     double              *ranges,            // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                     rangeAndPower       *rnp                // Output array of ray power at, and range to reciever
)

{
    cl_mem dTriangles ;
    cl_mem dTextures ;
    cl_mem dHits ;
    cl_mem dRays ;
    cl_mem dReflectedRays ;
    cl_mem dShadowRays ;
    cl_mem dRanges ;
    cl_mem dRnp ;
    
    size_t globalWorkSize ;
    
    // Make sure globalWorkSize is a multiple of localWorkSize
    //
    globalWorkSize = nRays ;
    while (globalWorkSize % localWorkSize)globalWorkSize++;
    
    int                nTriangles       = KDT.nTriangles;         // number of triangles in array 'triangles'
    Triangle *         triangles        = KDT.triangles;          // Array of triangles of size nTriangles
    int                nTextures        = KDT.nTextures;          // number of textures in array 'textures'
    Texture *          textures         = KDT.textures;           // Array of textures of size nTextures
   
    
    // Create OpenCl device buffers to hold data for this ray cast
    //
    
    dTriangles     = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Triangle)*nTriangles, NULL, &_err));
    dTextures      = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Texture)*nTextures,   NULL, &_err));
    dHits          = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Hit)*nRays,           NULL, &_err));
    dRays          = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Ray)*nRays,           NULL, &_err));
    dReflectedRays = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Ray)*nRays,           NULL, &_err));
    dShadowRays    = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Ray)*nRays,           NULL, &_err));
    dRanges        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(double)*nRays,        NULL, &_err));
    dRnp           = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(rangeAndPower)*nRays, NULL, &_err));
    
    
    // Add write buffer commands to the command queue ready for execution
    //
    CL_CHECK(clEnqueueWriteBuffer(Q, dTriangles,     CL_TRUE, 0, sizeof(Triangle)*nTriangles, triangles,     0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dTextures,      CL_TRUE, 0, sizeof(Texture)*nTextures,   textures,      0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dHits,          CL_TRUE, 0, sizeof(Hit)*nRays,           hits,          0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dRays,          CL_TRUE, 0, sizeof(Ray)*nRays,           rays,          0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dReflectedRays, CL_TRUE, 0, sizeof(Ray)*nRays,           reflectedRays, 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dShadowRays,    CL_TRUE, 0, sizeof(Ray)*nRays,           shadowRays,    0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dRanges,        CL_TRUE, 0, sizeof(double)*nRays,        ranges,        0, NULL, NULL));
    
    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(kernel, 0,   sizeof(cl_mem),      &dTriangles));
    CL_CHECK(clSetKernelArg(kernel, 1,   sizeof(cl_mem), &dTextures));
    CL_CHECK(clSetKernelArg(kernel, 2,   sizeof(cl_mem),   &dHits));
    CL_CHECK(clSetKernelArg(kernel, 3,   sizeof(SPVector),   &RxPos));
    CL_CHECK(clSetKernelArg(kernel, 4,   sizeof(double),   &GrxOverFourPi));
    CL_CHECK(clSetKernelArg(kernel, 5,   sizeof(int),   &nRays));
    CL_CHECK(clSetKernelArg(kernel, 6,   sizeof(cl_mem),   &dRays));
    CL_CHECK(clSetKernelArg(kernel, 7,   sizeof(cl_mem),   &dReflectedRays));
    CL_CHECK(clSetKernelArg(kernel, 8,   sizeof(cl_mem),   &dShadowRays));
    CL_CHECK(clSetKernelArg(kernel, 9,   sizeof(cl_mem),   &dRanges));
    CL_CHECK(clSetKernelArg(kernel, 10,  sizeof(cl_mem),   &dRnp));
                 
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(Q, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL));
    
    // Read results from device
    //
    CL_CHECK(clEnqueueReadBuffer(Q,dRnp, CL_TRUE,0,sizeof(rangeAndPower)*nRays,   rnp, 0,NULL,NULL));
    
    // Clear down OpenCL allocations
    //
    clReleaseMemObject(dTriangles);
    clReleaseMemObject(dTextures);
    clReleaseMemObject(dHits);
    clReleaseMemObject(dRays);
    clReleaseMemObject(dReflectedRays);
    clReleaseMemObject(dShadowRays);
    clReleaseMemObject(dRanges);
    clReleaseMemObject(dRnp);

    return ;
}