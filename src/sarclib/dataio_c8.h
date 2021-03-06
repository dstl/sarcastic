/***************************************************************************
 * 
 *           Module :  dataio_c8.h
 *          Program :  sarclib
 *       Created by :  Matt Nottingham on 05/10/2006
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      This file contains functions for reading c8 files made by IFP4
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

#ifndef sarclib_DATAIO_C8_H__
#define sarclib_DATAIO_C8_H__

#include "sarclib.h"

//// Parse the c8 image filename
///
SPStatus * im_c8_parse_fname(const char * fname, int * nx, int * ny, SPStatus * status);

/// Read in a complete c8 file from the given file name
///
SPStatus * im_load_c8(SPImage * im, const char * filename, SPStatus * status);

/// reads in a im->nx, im->ny subset of image from offset ox, oy
///
SPStatus * im_load_c8_subset(SPImage *im, const char * fname, long ox, long oy, SPStatus * status);

#endif

