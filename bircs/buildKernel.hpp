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

#ifndef buildKernel_hpp
#define buildKernel_hpp

#include <stdio.h>
#include "OpenCLUtils.h"

int buildKernel(cl_context context,             // OpenCL Context. Already created
                const char *kernelCodePath,     // String defining FQPN for kernel code
                const char *kernelCodeName,     // OpenCL kernel name for kernel code
                cl_device_id devID,             // Device Id on this platform to compile for
                int        workDims,            // WorkDimensions to optimise localWorkSize for
                cl_program *program,            // Output - Compiled OCL programme
                cl_kernel  *kernel,             // Output - OpenCL Kernel
                size_t     *localWorkSize       // Output - Prefered local worksize for kernel
) ;


#endif /* buildKernel_hpp */
