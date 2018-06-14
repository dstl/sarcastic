/** @file*************************************************************************
 *
 *       Module:    SIlib2.h
 *      Program:    SIlib2
 *   Created by:    Emma Griffiths 01/02/2005.
 *                  Darren Muff 21/07/2013
 *                  Darren Muff 29 Aug 2016
 *                  Copyright  (c) 2013 [dstl]. All rights reserved. 
 *
 *   Description:
 *      Originally mtlib.h. This file is the main header file for SIlib2
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

#ifndef SILib2_SILIB2_H__
#define SILib2_SILIB2_H__

#include <pthread.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/types.h>
#include <string.h>
#include <inttypes.h>
#include <unistd.h>

#if defined __CUFILE__ ||  defined __cplusplus
extern "C"
{
#endif

#include "Version.h"                    // Version info for SIlib2
#include "error_defn.h"                 // all the error definitions
#include "mtmaths.h"                    // macros that do some basic maths functions, square, max, min, sign
#include "cmplx.h"                      // complex functions and definition
#include "vector.h"                     // vector functions and definition
#include "matrix.h"                     // defines a matrix structure and some macros for working on them
#include "fft.h"                        // contains constants for the fft mode
#include "image_stats.h"                // the structure & defines used in statics & histograms
#include "image.h"                      // image structure and functions
#include "enum2string.h"                // turns enums into useful text strings
#include "dataio.h"                     // saving and loading and inputting of data from command line
#include "latlon.h"                     // Converts between ECEF & Lat/Longs
#include "image_interp.h"               // include file for Mark's code
#include "window.h"                     // kaiser-bessel window function
#include "math_macro.h"                 // wonderful math macros
#include "interpolate.h"                // Sinc interpolation routines
#include "physical.h"                   // Physical constants and useful units
#include "GeoConst.h"                   // Constants associated with the Earth
#include "GeoVector.h"                  // Routines for converting between Earth coordinate systems
#include "collectionGeometry.h"         // Header file for calculating collection geometries
#include "Timer.h"                      // Header file for timer routines
#include "alloc.h"                      // Header file for safe malloc / calloc macros
    
#define BETWEEN_PIS 1000
    
#define PATH "/local_storage/"
    
#ifndef TRUE
#define TRUE (1)
#endif
    
#ifndef FALSE
#define FALSE (0)
#endif
    
#if defined __CUFILE__ || defined __cplusplus
}
#endif

#endif

