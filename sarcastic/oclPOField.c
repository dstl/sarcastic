/***************************************************************************
 *
 *       Module:    oclPOField.c
 *      Program:    SARCASTIC
 *   Created by:    Darren on 25/01/2015.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Function to use opencl to calculate the E field at the reciever
 *      using Physical Optics
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  25/01/2015
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

typedef struct HitPoint {
    SPVector hit;       // Location of hitpoint in x,y,z
    int tri;            // index of triangle that this hit is on
} HitPoint ;

void oclPOField(cl_context          context,            // OpenCL context - already built
                cl_command_queue    Q,                  // OpenCl command Q - already instatiated
                cl_kernel           kernel,             // OpenCl kernel for this routine to call
                size_t              localWorkSize,      // Local workgroupsize to use for OpenCL Kernel
                Triangle            *triangles,         // Array of triangles
                int                 ntris,              // Number of triangles
                Hit                 *hits,              // Array of hit locations to x-ref with triangles for material props
                int                 nRays,              // The number of reflected rays being considered
                Ray                 *rays,              // unit vector rays arriving at hitpoint
                Ray                 *shadowRays,        // Array of reflected rays - used for their origin as its the reflection point to Rx
                SPVector            RxPos,              // Location of Receiver in x,y,z
                double              k,                  // Wavenumber constant k = 2 * PI / Lambda
                double              *ranges,            // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                double              gainRx,             // Receiver gain used for power calculations
                rangeAndPower       *rnp                // Output array of ray power at, and range to reciever
)

{
    cl_mem dHitpoints ;
    cl_mem dRays ;
    cl_mem dReflectedRays ;
    cl_mem dRanges ;
    cl_mem dEsVs;
    cl_mem dEsHs;
    cl_mem dTriangles ;
    
    size_t globalWorkSize ;
    
    // Make sure globalWorkSize is a multiple of localWorkSize
    //
    globalWorkSize = nRays ;
    while (globalWorkSize % localWorkSize)globalWorkSize++;
    
    // Create OpenCl device buffers to hold data for this ray cast
    //
    dTriangles     = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Triangle)*ntris,      NULL, &_err));
    dHitpoints     = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(HitPoint)*nRays,      NULL, &_err));
    dRays          = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Ray)*nRays,           NULL, &_err));
    dReflectedRays = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(Ray)*nRays,           NULL, &_err));
    dRanges        = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_ONLY,  sizeof(double)*nRays,        NULL, &_err));
    dEsVs          = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(SPCmplx)*nRays,       NULL, &_err));
    dEsHs          = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(SPCmplx)*nRays,       NULL, &_err));
    
    // Prepare the hitpoints from the input 'Hits' (which contain the index of the triangle assiciated with each hit) and
    // the x,y,z location of the hit (taken from the input shadow ray origin).
    //
    HitPoint *hitpoints = (HitPoint *)sp_malloc(sizeof(HitPoint)*nRays);
    int *hitsOnEachTri  = (int *)sp_calloc(ntris, sizeof(int));
    
    for(int i=0; i<nRays; i++){
        hitpoints[i].hit = shadowRays[i].org ;
        hitpoints[i].tri = hits[i].trinum ;
        hitsOnEachTri[hitpoints[i].tri]++;
        rays[i].pow = rays[i].pow; // (4 * SIPC_pi * hits[i].dist * hits[i].dist) ;
    }
    
    // Add write buffer commands to the command queue ready for execution
    //
    CL_CHECK(clEnqueueWriteBuffer(Q, dTriangles,     CL_TRUE, 0, sizeof(Triangle)*ntris,      triangles,     0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dHitpoints,     CL_TRUE, 0, sizeof(HitPoint)*nRays,      hitpoints,     0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dRays,          CL_TRUE, 0, sizeof(Ray)*nRays,           rays,          0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dReflectedRays, CL_TRUE, 0, sizeof(Ray)*nRays,           shadowRays, 0, NULL, NULL));
    CL_CHECK(clEnqueueWriteBuffer(Q, dRanges,        CL_TRUE, 0, sizeof(double)*nRays,        ranges,        0, NULL, NULL));
    
    // Define unit vectors for V & H fields
    // We do this as the definition of H and V from the sensor may not be truly horizontal
    // or vertical from the sensor. By allowing us to define the V & H directions we can
    // accurately model the V & H at the sensor.
    // We also do it here as we dont want the OCL kernel calculating this each time it is called
    //
    SPVector zhat, RXVdir, RXHdir, obsDir ;
    VECT_CREATE(0, 0, 1, zhat) ;
    VECT_NORM(RxPos, obsDir);
    VECT_CROSS(zhat, obsDir, RXHdir) ;
    VECT_CROSS(obsDir, RXHdir, RXVdir) ;
    VECT_NORM(RXVdir, RXVdir) ;
    VECT_NORM(RXHdir, RXHdir) ;
    
    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(kernel, 0,   sizeof(cl_mem),   &dTriangles));   // input array of triangles
    CL_CHECK(clSetKernelArg(kernel, 1,   sizeof(cl_mem),   &dRays));        // input array of rays
    CL_CHECK(clSetKernelArg(kernel, 2,   sizeof(int),      &nRays));        // input number of incident rays to process
    CL_CHECK(clSetKernelArg(kernel, 3,   sizeof(cl_mem),   &dHitpoints));   // input array of hit points for each ray
    CL_CHECK(clSetKernelArg(kernel, 4,   sizeof(SPVector), &RxPos));        // input receiver location to calc field at
    CL_CHECK(clSetKernelArg(kernel, 5,   sizeof(double),   &k));            // input k=2*PI/lambda
    CL_CHECK(clSetKernelArg(kernel, 6,   sizeof(SPVector), &RXVdir));       // input unit vector defining V pol at receiver
    CL_CHECK(clSetKernelArg(kernel, 7,   sizeof(SPVector), &RXHdir));       // input unit vector defining H pol at receiver
    CL_CHECK(clSetKernelArg(kernel, 8,   sizeof(cl_mem),   &dEsVs));        // output array of scattered field strengths for V pol
    CL_CHECK(clSetKernelArg(kernel, 9,   sizeof(cl_mem),   &dEsHs));        // output array of scattered field strengths for H pol
    
    
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(Q, kernel, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL));
    
    // Read results from device
    //
    SPCmplx *VsTmp, *HsTmp;
    VsTmp = sp_malloc(sizeof(SPCmplx) * nRays);
    HsTmp = sp_malloc(sizeof(SPCmplx) * nRays);
    
    CL_CHECK(clEnqueueReadBuffer(Q,dEsVs, CL_TRUE,0,sizeof(SPCmplx)*nRays, VsTmp, 0,NULL,NULL));
    CL_CHECK(clEnqueueReadBuffer(Q,dEsHs, CL_TRUE,0,sizeof(SPCmplx)*nRays, HsTmp, 0,NULL,NULL));
    
    // The PO Field equation calculations take into account the inverse square law path loss
    // back to the receiver. We do need to multiply by the receiver gain though.
    // We also need to take into account how many times each triangle has been hit as the PO
    // calculations assume RCS from one point. 
    //
    for (int r=0; r<nRays; r++ ){
        rnp[r].range = (ranges[r] + rays[r].len + hits[r].dist) / 2;
        CMPLX_SCMULT(1.0 / hitsOnEachTri[hitpoints[r].tri], VsTmp[r], rnp[r].Es) ;
//        printf("Range to Rx: %f, ray.len: %f, hit.dist: %f (Range : %f) E: %e (normalised E: %e %d tris)\n",
//               ranges[r],rays[r].len,hits[r].dist,rnp[r].range,CMPLX_MAG(VsTmp[r]),CMPLX_MAG(rnp[r].Es),hitsOnEachTri[hitpoints[r].tri]);

    }
    
    // Clear down OpenCL allocations
    //
    clReleaseMemObject(dTriangles);
    clReleaseMemObject(dHitpoints);
    clReleaseMemObject(dRays);
    clReleaseMemObject(dReflectedRays);
    clReleaseMemObject(dRanges);
    clReleaseMemObject(dEsVs);
    clReleaseMemObject(dEsHs);
    
    free(hitpoints);
    free(VsTmp);
    free(HsTmp);
    free(hitsOnEachTri);
    
    return ;
}