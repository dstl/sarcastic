/***************************************************************************
 *
 *       Module:    fft.h
 *      Program:    SIlib2
 *   Created by:    Matt Nottingham on 04/01/2005.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains information about FFT constants
 *      The FFT is controlled by bitwise addition ('&') of mode parameters
 *      Available modes for FFT's are:
 *         
 *      FWD, REV                       - forward or reverse.
 *      NOSCALE, SCALE_N, SCALE_R      - No scaling, Scale by N points or SQRT(N points)
 *      FFT_X_ONLY, FFT_Y_ONLY, FFT_2D - Perform FFT in X direction, Y direction or both
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  19/07/2013
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


#ifndef SILib2_FFT_H__
#define SILib2_FFT_H__

#include "SIlib2.h"

typedef enum
{
	FWD		    = 0, 	REV 		= 4,
   	NOSCALE		= 0, 	SCALE_N 	= 16,  SCALE_R	= 32,
	FFT_X_ONLY 	= 64,	FFT_Y_ONLY 	= 128, FFT_2D	= 192,
	FFT_FFTW_ESTIMATE    = 0,    FFT_FFTW_MEASURE = 256, FFT_FFTW_PATIENT = 512
} FFTMODE; 

int fft_c8(SPCmplx *data, int ntot, int n, int nspan, int isn);

#endif

