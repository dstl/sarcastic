/***************************************************************************
 *
 *       Module:    tdpcore.h
 *      Program:    tdpocl
 *   Created by:    Darren on 24/05/2013.
 *                  Copyright (c) 2013 [Dstl]. All rights reserved.
 *
 *   Description:
 *   <ENTER DESCRIPTION HERE>
 *
 *
 *   CLASSIFICATION        :  <PENDING>
 *   Date of CLASSN        :  14/03/2013
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

#ifndef tdpocl_tdpcore_h
#define tdpocl_tdpcore_h
#include <sarclib/sarclib.h>
#if defined (__APPLE__) || defined(MACOSX)
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif
#define NPOINTS (8)
#define OVERSAMP (512)
#include "OpenCLUtils.h"

typedef int Runtype;
#define RUNONLY  0  // Call kernel without building device buffers
#define TDPINIT  1  // Initialise device buffers
#define TDPCLOSE 2  // Close device buffers and recover memory

typedef struct {
    Runtype run ;                   // 0 to process image; 1 to initialise (or reininitalise) OCL device buffers; 2 to shut down
    size_t globalWorkSize [2] ;     // Device global work size for this run
    size_t localWorkSize  [2] ;     // Device local work size for this run
    OCLTask task_s;                 // Structure containing ocl contexts, programs, kernels, command queues etc
    SPImage *surfaces;              // array of arrays (one per device) of ECEF positions defining the output image
    SPImage phd;                    // Range compressed phase history data used to build the image ( x is samples, y is pulses)
    SPImage *images;                // array of output images (one per device) - same size as surface
    SPVector *TXs;                  // array of transmit positions for this run (phd->ny elements)
    SPVector *RXs;                  // array of receieve positions for this run (phd->ny elements)
    SPVector *SRPs;                 // array of sensor aim points for each line in PHD (phd->ny elements)
    double *SRPRanges ;             // array of pre-calculated ranges from pulse Tx to SRO to pulseRX (phd->ny elements)
    float *FXStep;                  // array of ADC sample sizes in frequency (phd->ny elements)
    int phseSgn;                    // value defining direction of increasing range
    double fcent;                   // value defining the SAR centre frequency
    float * ikernel;                // An interpolation kernel of length (NPOINTS*OVERSAMP)+1
    cl_mem *dSurf;                  // Pointer to device mem buffer array for surface
    cl_mem *dPHD;                   // Pointer to device mem buffer array for phd
    cl_mem *dIm;                    // Pointer to device mem buffer array for output image
    cl_mem *dTXs;                   // Pointer to device mem buffer array for TXs
    cl_mem *dRXs;                   // Pointer to device mem buffer array for RXs
    cl_mem *dSRPs;                  // Pointer to device mem buffer array for SRPs
    cl_mem *dSRPRanges;             // Pointer to device mem buffer array for pulse to SRP ranges
    cl_mem *dFXStep;                // Pointer to device mem buffer array for FXstep
    cl_mem *diKernel;               // Pointer to device mem buffer array for interpolation kernel
} TDPCoreStruct;

void tdpcore(TDPCoreStruct *tdp_s); // Structure for configuring the run
            

#endif
