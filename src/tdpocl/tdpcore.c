/***************************************************************************
 * 
 *           Module :  tdpcore.c
 *          Program :  tdpocl
 *       Created by :  Darren Muff on 24/05/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *         Core thread to process TDP pulses
 *
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
#include <stdio.h>
#include "tdpcore.h"

void tdpcore(TDPCoreStruct *tdp_s) {
    
    SPStatus status;
    status.status = 0;
    status.debug  = 0;
    OCLTask task ;
    CHECK_STATUS_NON_PTR(status);
    int x,y, imNx,imNy ;

    task = tdp_s->task_s ;
    
    
    size_t lws_t ;
    size_t warpSize;
    
    imNx = (int)tdp_s->images[0].nx ;
    imNy = (int)tdp_s->images[0].ny ;

    // Calculate maximum local work size that this kernel can support on this device
    //
    CL_CHECK(clGetKernelWorkGroupInfo(task.kernels[0], task.device_ids[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(lws_t), &lws_t, NULL));
    
    // What do we think is the best performance size though for this kernel, based upon this device's properties
    //
    CL_CHECK(clGetKernelWorkGroupInfo(task.kernels[0], task.device_ids[0], CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(warpSize), &warpSize, NULL));

    // For experimentation to find a good Local worksize then uncomment these lines and
    // insert your own values
    //
    //    tdp_s->localWorkSize[0] = PUT_X_LWS_HERE ;
    //    tdp_s->localWorkSize[1] = PUT_Y_LWS_HERE ;
    
    tdp_s->globalWorkSize[0] = imNx;
    tdp_s->globalWorkSize[1] = imNy;
    int calcLWS = 1;
    if (tdp_s->localWorkSize[0] > 0 && tdp_s->localWorkSize[1] > 0){
        calcLWS = 0;
    }

    if (calcLWS){
        // Build work sizes on devices and check that devices are identical.
        //
        int tx,ty;
        for(int dev=0; dev<task.devsToUse; dev++){
            best2DWorkSize(task.kernels[dev], task.device_ids[dev], imNx, imNy,&tx,&ty, &status);
            if( dev !=0){
                if(tx != tdp_s->localWorkSize[0] || ty != tdp_s->localWorkSize[1]){
                    printf("ERROR : device[%d] is different to device[0] - d[%d]localWorkSize: %d, %d; d[0]localWorkSize: %zd %zd\n",dev,dev,tx,ty,tdp_s->localWorkSize[0],tdp_s->localWorkSize[1]);
                    exit(-1);
                }
            }
            tdp_s->localWorkSize[0] = tx;
            tdp_s->localWorkSize[1] = ty;
        }
    }
    
    x = imNx;
    y = imNy;
    
    // grow globalworksize until its a multiple of the local workgroup size
    // This is for efficiency on the device. The device kernel doesnt process any extra bits.
    //
    while ( x % tdp_s->localWorkSize[0] )x++;
    while ( y % tdp_s->localWorkSize[1] )y++;
    tdp_s->globalWorkSize[0] = x ;
    tdp_s->globalWorkSize[1] = y ;
    
    if(tdp_s->run & TDPINIT){  // Initialise those buffer items that are based upon size of tdp_s->im
        
        // Uncomment the following for working out warp sizes
        //
        printf("  Device max local work size         : %zu\n",lws_t);
        printf("  Device prefered warp size          : %zu\n", warpSize);
        printf("  Calculated device work size        : %zu,%zu\n",tdp_s->localWorkSize[0], tdp_s->localWorkSize[1]);
        printf("  Calculated device global data size : %zu,%zu\n",tdp_s->globalWorkSize[0], tdp_s->globalWorkSize[1]);

        // Each Chunk goes through the following OCL steps, each one for each device:
        //   Create kernels, context, program and command queues
        //   Allocate memory for device buffer
        //   Enqueue data transmission to each device
        //   Execute kernel on each device
        //   Read results from each device
        //   Release memory objects for each device
        
        /*
         // Handy for debug so keep in for the moment
         //
         printf("Memory reqs: %lldMB\n",(sizeof(SPCmplx)*tdp_s->phd.nx*tdp_s->phd.ny+
         sizeof(SPVector)*tdp_s->phd.ny+
         sizeof(SPVector)*tdp_s->phd.ny+
         sizeof(SPVector)*tdp_s->phd.ny+
         sizeof(float)*tdp_s->phd.ny+
         sizeof(SPVector)*tdp_s->surface.nx*tdp_s->surface.ny+
         sizeof(SPCmplx)*imNx*imNy+
         sizeof(float)*(OVERSAMP*NPOINTS+1))/1024/1024);*/
        
        // Allocate buffers on devices
        //
        for (int dev=0; dev<task.devsToUse; dev++){
            
            tdp_s->dPHD[dev]       = CL_CHECK_ERR(clCreateBuffer(task.contexts[dev], CL_MEM_READ_ONLY, sizeof(SPCmplx)*tdp_s->phd.nx*tdp_s->phd.ny, NULL, &_err));
            tdp_s->dRXs[dev]       = CL_CHECK_ERR(clCreateBuffer(task.contexts[dev], CL_MEM_READ_ONLY, sizeof(SPVector)*tdp_s->phd.ny, NULL, &_err));
            tdp_s->dTXs[dev]       = CL_CHECK_ERR(clCreateBuffer(task.contexts[dev], CL_MEM_READ_ONLY, sizeof(SPVector)*tdp_s->phd.ny, NULL, &_err));
            tdp_s->dSRPRanges[dev] = CL_CHECK_ERR(clCreateBuffer(task.contexts[dev], CL_MEM_READ_ONLY, sizeof(double)*tdp_s->phd.ny,   NULL, &_err));
            tdp_s->dFXStep[dev]    = CL_CHECK_ERR(clCreateBuffer(task.contexts[dev], CL_MEM_READ_ONLY, sizeof(float)*tdp_s->phd.ny,    NULL, &_err));
            tdp_s->dSurf[dev]      = CL_CHECK_ERR(clCreateBuffer(task.contexts[dev], CL_MEM_READ_ONLY, sizeof(SPVector)*imNx*imNy,     NULL, &_err));
            tdp_s->dIm[dev]        = CL_CHECK_ERR(clCreateBuffer(task.contexts[dev], CL_MEM_READ_WRITE,sizeof(SPCmplx)*imNx*imNy,      NULL, &_err));
            tdp_s->diKernel[dev]   = CL_CHECK_ERR(clCreateBuffer(task.contexts[dev], CL_MEM_READ_ONLY, sizeof(float)*(OVERSAMP*NPOINTS+1),NULL, &_err));
        }
        
        tdp_s->run -= TDPINIT ;
        
    }
    
    // Enqueue data transmission for each device
    // perform this each time to load up new data
    //
    for (int dev=0; dev<task.devsToUse; dev++){
        
        CL_CHECK(clEnqueueWriteBuffer(task.commands[dev],  tdp_s->dPHD[dev],    CL_TRUE, 0, sizeof(SPCmplx)*tdp_s->phd.nx*tdp_s->phd.ny, tdp_s->phd.data.cmpl_f, 0, NULL, NULL));
        CL_CHECK(clEnqueueWriteBuffer(task.commands[dev],  tdp_s->dRXs[dev],    CL_TRUE, 0, sizeof(SPVector)*tdp_s->phd.ny, tdp_s->RXs,                0, NULL, NULL));
        CL_CHECK(clEnqueueWriteBuffer(task.commands[dev],  tdp_s->dTXs[dev],    CL_TRUE, 0, sizeof(SPVector)*tdp_s->phd.ny, tdp_s->TXs,                0, NULL, NULL));
        CL_CHECK(clEnqueueWriteBuffer(task.commands[dev],  tdp_s->dSRPRanges[dev], CL_TRUE, 0, sizeof(double)*tdp_s->phd.ny, tdp_s->SRPRanges,         0, NULL, NULL));
        CL_CHECK(clEnqueueWriteBuffer(task.commands[dev],  tdp_s->dFXStep[dev], CL_TRUE, 0, sizeof(float)*tdp_s->phd.ny,    tdp_s->FXStep,             0, NULL, NULL));
        CL_CHECK(clEnqueueWriteBuffer(task.commands[dev],  tdp_s->dSurf[dev],   CL_TRUE, 0, sizeof(SPVector)*imNx*imNy,tdp_s->surfaces[dev].data.vect, 0, NULL, NULL));
        CL_CHECK(clEnqueueWriteBuffer(task.commands[dev],  tdp_s->diKernel[dev],CL_TRUE, 0, sizeof(float)*(OVERSAMP*NPOINTS+1),tdp_s->ikernel,         0, NULL, NULL));
        CL_CHECK(clEnqueueWriteBuffer(task.commands[dev],  tdp_s->dIm[dev],     CL_TRUE, 0, sizeof(SPCmplx)*imNx*imNy,tdp_s->images[dev].data.cmpl_f,  0, NULL, NULL));
    }
    
    // Now queue up the kernel
    // Perform on each run
    //
    for (int dev=0; dev<task.devsToUse; dev++){
        
        CL_CHECK(clSetKernelArg(task.kernels[dev], 0,  sizeof(cl_mem), &(tdp_s->dPHD[dev])));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 1,  sizeof(cl_mem), &(tdp_s->dIm[dev])));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 2,  sizeof(int), &(imNx)));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 3,  sizeof(int), &(imNy)));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 4,  sizeof(cl_mem), &(tdp_s->dRXs[dev])));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 5,  sizeof(cl_mem), &(tdp_s->dTXs[dev])));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 6,  sizeof(cl_mem), &(tdp_s->dSRPRanges[dev])));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 7,  sizeof(cl_mem), &(tdp_s->dSurf[dev])));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 8,  sizeof(int), &(tdp_s->phseSgn)));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 9,  sizeof(cl_mem), &(tdp_s->dFXStep[dev])));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 10,  sizeof(int), &(tdp_s->phd.nx)));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 11,  sizeof(double), &(tdp_s->fcent)));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 12, sizeof(cl_mem), &(tdp_s->diKernel[dev])));
        CL_CHECK(clSetKernelArg(task.kernels[dev], 13, sizeof(int), &(tdp_s->phd.ny)));
        
        CL_CHECK(clEnqueueNDRangeKernel(task.commands[dev], task.kernels[dev], 2, NULL, tdp_s->globalWorkSize, tdp_s->localWorkSize, 0, NULL, NULL));
    }
    // Read results form image buffers
    //
    for (int dev=0; dev<task.devsToUse; dev++){
        
        CL_CHECK(clEnqueueReadBuffer(task.commands[dev], tdp_s->dIm[dev], CL_TRUE, 0, sizeof(SPCmplx)*imNx*imNy, &(tdp_s->images[dev].data.cmpl_f[0]), 0, NULL, NULL));
    }
    
    //    for(int y=0; y<imNy; y++){
    //        for(int x=0; x<imNx; x++){
    //            printf("results[%d,%d] : %f,%f\n",x,y,tdp_s->im.data.cmpl_f[y*imNx+x].r,tdp_s->im.data.cmpl_f[y*imNx+x].i);
    //        }
    //    }
    
    if(tdp_s->run & TDPCLOSE ){
        
        for (int dev=0; dev<task.devsToUse; dev++){
            clReleaseMemObject(tdp_s->dPHD[dev]);
            clReleaseMemObject(tdp_s->dRXs[dev]);
            clReleaseMemObject(tdp_s->dTXs[dev]);
            clReleaseMemObject(tdp_s->dSRPRanges[dev]);
            clReleaseMemObject(tdp_s->dFXStep[dev]);
            clReleaseMemObject(tdp_s->dSurf[dev]);
            clReleaseMemObject(tdp_s->dIm[dev]);
            clReleaseMemObject(tdp_s->diKernel[dev]);
        }
        
        
    }
    return ;
}
