/** @file********************************************************************
 *
 *       Module:    dataio_envi.h
 *      Program:    SIlib2
 *   Created by:    Matt Nottingham on 21/09/2007.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains the prototypes for the functions (def'ed in dataio_envi.c)
 *      for reading and writing ENVI files
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

#ifndef SILib2_DATAIO_ENVI_H__
#define SILib2_DATAIO_ENVI_H__

#include "SIlib2.h"

typedef enum {ENVI_ITYPE_UINT8 = 1,        ENVI_ITYPE_INT16 = 2,   ENVI_ITYPE_INT32 = 3,
              ENVI_ITYPE_FLOAT = 4,        ENVI_ITYPE_DOUBLE = 5,  ENVI_ITYPE_CMPL_FLOAT = 6,
              ENVI_ITYPE_CMPL_DOUBLE =  9, ENVI_ITYPE_UINT16 = 12, ENVI_ITYPE_UINT32 = 13,
              ENVI_ITYPE_INT64 = 14,       ENVI_ITYPE_UINT64 = 15, ENVI_ITYPE_UNKNOWN = 0} EnviImageType;

/// Loads up a complete image from file
///
SPStatus * im_load_envi (SPImage *a, const char *fname, SPStatus *status);

/// Loads up a subset of an image from a file
///
SPStatus * im_load_envi_subset (SPImage *a, const char *fname, int64_t ox, int64_t oy, SPStatus *status);

/// Loads up just the metadata from file
///
SPStatus* im_load_envi_metadata (SPImage *a, const char *fname, SPStatus *status); 

SPImageType im_conv_from_envi_type(EnviImageType itype);
EnviImageType im_conv_to_envi_type(SPImageType itype);

#endif

