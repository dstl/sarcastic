/***************************************************************************
 *
 *       Module:    oclRayTrace.c
 *      Program:    GOSS
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Function to use opencl to ray trace a single pulse
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

#include "GOSS.h"

void oclRayTrace(cl_context         context,            // OpenCL context - already instantiated
                 cl_command_queue   Q,                  // OpenCl command Q - already instatiated
                 cl_kernel          RTkernel,           // Ray Tracing kernel, already compiled and created from a program
                 size_t             globalWorkSize[2],  // Total number of rays in each dimension to cast
                 size_t             localWorkSize[2],   // Work dimensions for this device
                 
                 int                nTriangles,         // number of triangles in array 'triangles'
                 Triangle *         triangles,          // Array of triangles of size nTriangles
                 int                nTextures,          // number of textures in array 'textures'
                 Texture *          textures,           // Array of textures of size nTextures
                 int                nTreeNodes,         // number of nodes in KdTree
                 KdData *           KdTree,             // SAH - KdTree to optimise ray traversal through volume
                 int                triListDataSize,    // size of trianglelist data
                 int *              triangleListData,   // array of triangle indices into Triangles
                 int                nLeaves,            // number of leaves (nodes with triangles) in KdTree
                 int *              triangleListPtrs,   // array of pointers into triangleListData for each KdTree node
                 
                 SPVector           RxPos,              // Receive position of radar (used to calculate range)
                 double             raySolidAng,        // ray solid angle
                 double             TxPowPerRay,        // TxPowerPerRay
                 double             Aeff,               // Antenna efficiancy required 
                 AABB               SceneBoundingBox,   // Bounding box of scene - required for ray traversal optimisation
                 
                 int                bounceToShow,       // Useful for debugging
                 int                pulseIndex,          // Useful for debugging
                 
                 Ray *              rays,               // Array of rays (size nAzBeam*nElBeam). Each ray will be cast through KdTree
                 rangeAndPower      *rnp                // output array of ranges and powers for each ray intersection
                )
{
    cl_mem dTriangles ;
    cl_mem dTextures ;
    cl_mem dKdTree ;
    cl_mem dtriListData ;
    cl_mem dtriListPtrs ;
    cl_mem drnp ;
    cl_mem dRays ;
    
    int nRays;
    int nAzBeam;            // Number of rays in 'rays' in azimuth direction
    int nElBeam;            // Number of rays in 'rays' in elevation direction

    nAzBeam = (int)globalWorkSize[0];
    nElBeam = (int)globalWorkSize[1];
    nRays   = nAzBeam * nElBeam ;
    
    // zero rnp array so that we can see which rays hit something
    //
    memset(rnp, 0, sizeof(rangeAndPower)*nRays*MAXBOUNCES);
    
    // Create OpenCl device buffers to hold data for this ray cast
    //
  
    dTriangles   = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Triangle)*nTriangles,            NULL, &_err));
    dTextures    = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Texture)*nTextures,              NULL, &_err));
    dKdTree      = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(KdData)*nTreeNodes,              NULL, &_err));
    dtriListData = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  triListDataSize,                        NULL, &_err));
    dtriListPtrs = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(int)*nLeaves,                    NULL, &_err));
    drnp         = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(rangeAndPower)*nRays*MAXBOUNCES, NULL, &_err));
    dRays        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Ray)*nRays,                      NULL, &_err));

    // Add write buffer commands to the command queue ready for execution
    //
    CL_CHECK(clEnqueueWriteBuffer(Q, dTriangles,   CL_TRUE, 0, sizeof(Triangle)*nTriangles, triangles,        0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dTextures,    CL_TRUE, 0, sizeof(Texture)*nTextures,   textures,         0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dKdTree,      CL_TRUE, 0, sizeof(KdData)*nTreeNodes,   KdTree,           0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dtriListData, CL_TRUE, 0, triListDataSize,             triangleListData, 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dtriListPtrs, CL_TRUE, 0, sizeof(int)*nLeaves,         triangleListPtrs, 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dRays,        CL_TRUE, 0, sizeof(Ray)*nRays,           rays,             0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, drnp,         CL_TRUE, 0, sizeof(rangeAndPower)*nRays*MAXBOUNCES,rnp,    0,NULL,NULL));

    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(RTkernel, 0,   sizeof(int), &nAzBeam));
    CL_CHECK(clSetKernelArg(RTkernel, 1,   sizeof(int), &nElBeam));
    CL_CHECK(clSetKernelArg(RTkernel, 2,   sizeof(SPVector), &RxPos));
    CL_CHECK(clSetKernelArg(RTkernel, 3,   sizeof(double), &raySolidAng));
    CL_CHECK(clSetKernelArg(RTkernel, 4,   sizeof(double), &TxPowPerRay));
    CL_CHECK(clSetKernelArg(RTkernel, 5,   sizeof(AABB), &SceneBoundingBox));
    CL_CHECK(clSetKernelArg(RTkernel, 6,   sizeof(double), &Aeff));
    CL_CHECK(clSetKernelArg(RTkernel, 7,   sizeof(int), &bounceToShow));
    CL_CHECK(clSetKernelArg(RTkernel, 8,   sizeof(cl_mem), &dRays));
    CL_CHECK(clSetKernelArg(RTkernel, 9,   sizeof(cl_mem), &dTriangles));
    CL_CHECK(clSetKernelArg(RTkernel, 10,  sizeof(cl_mem), &dTextures));
    CL_CHECK(clSetKernelArg(RTkernel, 11,  sizeof(cl_mem), &dKdTree));
    CL_CHECK(clSetKernelArg(RTkernel, 12,  sizeof(cl_mem), &dtriListData));
    CL_CHECK(clSetKernelArg(RTkernel, 13,  sizeof(cl_mem), &dtriListPtrs));
    CL_CHECK(clSetKernelArg(RTkernel, 14,  sizeof(cl_mem), &drnp));
    CL_CHECK(clSetKernelArg(RTkernel, 15,  sizeof(int), &pulseIndex));
    
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(Q, RTkernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL));
    
    // Read results from device
    //
    CL_CHECK(clEnqueueReadBuffer(Q,drnp, CL_TRUE,0,sizeof(rangeAndPower)*nRays*MAXBOUNCES,rnp,0,NULL,NULL));
    
    // Clear down OpenCL allocations
    //
    clReleaseMemObject(dTriangles);
    clReleaseMemObject(dTextures);
    clReleaseMemObject(dKdTree);
    clReleaseMemObject(dtriListData);
    clReleaseMemObject(dtriListPtrs);
    clReleaseMemObject(dRays);
    clReleaseMemObject(drnp);
    return ;
}