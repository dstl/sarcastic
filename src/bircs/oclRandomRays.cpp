/***************************************************************************
 * 
 *           Module :  oclRandomRays.cpp
 *          Program :  bircs
 *       Created by :  Darren Muff on 12/06/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *     Function to build an array of rays for bircs using OpenCL
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

#include "oclRandomRays.hpp"
#include "threadCore.hpp"

void oclRandomRays(cl_context context,          // OpenCL context - alrready built
                   cl_command_queue Q,          // OpenCl command Q - already instatiated
                   cl_kernel RRkernel,          // OpenCl kernel for this routine to call
                   int nAzBeam,                 // Number of azimuth rays in beam
                   int nElBeam,                 // Number of elevation rays in beam
                   size_t localWorkSize[2],     // local work size for this device - prealculated
                   double azStdDev,             // Az standard deviation of rays
                   double elStdDev,             // El standard deviation of rays
                   SPVector origin,             // Location of origin or rays
                   SPVector aimpoint,           // Mean aimpoint for rays
                   double raypow,               // Ray power (Pp = (Pt * Gtx / (4*PI)) * dAz * dEl)
                   Ray *rayArray                // Array that will be returned
){
    cl_mem dRays;
    size_t globalWorkSize[2];
    int x,y,nRays ;
    unsigned long seed ;
    
    x = nAzBeam ;
    y = nElBeam ;
    nRays   = nAzBeam * nElBeam ;
    
    // Seed a random number with the coordinates of the Transmitter
    //
    srandom(origin.x+origin.y+origin.z);
    seed = random() ;
    
    // grow globalworksize until its a multiple of the local workgroup size
    //
    while ( x % localWorkSize[0] )x++;
    while ( y % localWorkSize[1] )y++;
    globalWorkSize[0] = x ;
    globalWorkSize[1] = y ;
    
    dRays = CL_CHECK_ERR(clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(Ray)*nRays, NULL, &_err));
    
    // Set up the kernel arguments
    //
    CL_CHECK(clSetKernelArg(RRkernel, 0,   sizeof(unsigned long), &seed));
    CL_CHECK(clSetKernelArg(RRkernel, 1,   sizeof(double), &azStdDev));
    CL_CHECK(clSetKernelArg(RRkernel, 2,   sizeof(double), &elStdDev));
    CL_CHECK(clSetKernelArg(RRkernel, 3,   sizeof(SPVector), &origin));
    CL_CHECK(clSetKernelArg(RRkernel, 4,   sizeof(SPVector), &aimpoint));
    CL_CHECK(clSetKernelArg(RRkernel, 5,   sizeof(int), &nAzBeam));
    CL_CHECK(clSetKernelArg(RRkernel, 6,   sizeof(int), &nElBeam));
    CL_CHECK(clSetKernelArg(RRkernel, 7,   sizeof(double), &raypow));
    CL_CHECK(clSetKernelArg(RRkernel, 8,   sizeof(cl_mem), &dRays));
    
    // Queue up the commands on the device to be executed as soon as possible
    //
    CL_CHECK(clEnqueueNDRangeKernel(Q, RRkernel, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL));
    
    // Read results from device
    //
    CL_CHECK(clEnqueueReadBuffer(Q,dRays, CL_TRUE,0,sizeof(Ray)*nRays,rayArray,0,NULL,NULL));
    
    clReleaseMemObject(dRays) ;
    
    return ;
    
}


