/***************************************************************************
 *
 *       Module:    buildKernel.hpp
 *      Program:    bircs
 *   Created by:    Darren Muff on 12/06/2017.
 *                  Copyright (c) 2017 [Dstl]. All rights reserved.
 *
 *   Description:
 *     Function to build an openCL kernel for BIRCS
 *
 *
 *   CLASSIFICATION        :  OFFICIAL
 *   Date of CLASSN        :  12/06/2017
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
#include "buildKernel.hpp"

int buildKernel(cl_context context,             // OpenCL Context. Already created
                const char *kernelCodePath,     // String defining FQPN for kernel code
                const char *kernelCodeName,     // OpenCL kernel name for kernel code
                cl_device_id devID,             // Device Id on this platform to compile for
                int        workDims,            // WorkDimensions to optimise localWorkSize for
                cl_program *program,            // Output - Compiled OCL programme
                cl_kernel  *kernel,             // Output - OpenCL Kernel
                size_t     *localWorkSize       // Output - Prefered local worksize for kernel
)
{
    char compilerOptions[255];
    cl_int err ;
    size_t warpSize;
    
    //    sprintf(compilerOptions, "-Werror -I/Users/darren/Development/sarcastic/sarcastic/kernels/") ;
    sprintf(compilerOptions, "-I%s",KERNELDIR) ;
    
    *program = CL_CHECK_ERR(clCreateProgramWithSource(context, 1, (const char **) &kernelCodePath, NULL, &_err));
    
    err = clBuildProgram(*program, 0, NULL, compilerOptions, NULL, NULL);
    if (err != CL_SUCCESS){
        size_t len;
        int buffsize = 32768;
        char buffer[buffsize] ;
        printf("*** Error: Failed to build program executable for kernel \"%s\" \n",kernelCodeName);
        clGetProgramBuildInfo(*program, devID, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("err: %d. Buffer length is %zu bytes:\n",err,len);
        if (len > buffsize) {
            char bigbuffer[len] ;
            clGetProgramBuildInfo(*program, devID, CL_PROGRAM_BUILD_LOG, sizeof(bigbuffer), bigbuffer, &len);
            printf("%s\n", bigbuffer);
        }else{
            printf("%s\n", buffer);
        }
        return (err) ;
    }
    
    *kernel   = CL_CHECK_ERR(clCreateKernel(*program, kernelCodeName, &_err));
    
    // Calculate maximum local work size that this kernel can support on this device
    //
    size_t lws_t ;
    CL_CHECK(clGetKernelWorkGroupInfo(*kernel, devID, CL_KERNEL_WORK_GROUP_SIZE, sizeof(lws_t), &lws_t, NULL));
    
    // What do we think is the best performance size though for this kernel, based upon this device's properties
    //
    CL_CHECK(clGetKernelWorkGroupInfo(*kernel, devID, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE, sizeof(warpSize), &warpSize, NULL));
    
    // For efficiency reduce the local workgroup size until it is a multiple of warpSize
    // As RandomRays requires 2D then do this by shrinking y dimension
    //
    int x,y ;
    switch (workDims){
        case 1:
            x = (int)lws_t ;
            while ( x % warpSize ) x--;
            *localWorkSize = x ;
            break ;
            
        case 2:
            x = (int)sqrt(lws_t);
            y = x;
            while ( (y > 1) && ((x * y) % warpSize) ) {       // Shrink y until it fits optimised warp size
                y--;
            }
            if ( y == 1 && ((x * y) % warpSize) ) {           // If we shrank y to 1 and still didnt optimise shrink x
                while ( (x > 1) && ((x * y) % warpSize) ) {
                    x--;
                }
            }
            localWorkSize[0] = x ;
            localWorkSize[1] = y ;
            break;
            
        default:
            printf("Unsupported work dimensions trying to build kernel \"%s\" \n",kernelCodeName);
            return CL_INVALID_WORK_GROUP_SIZE ;
    }
    
    return CL_SUCCESS ;
}



