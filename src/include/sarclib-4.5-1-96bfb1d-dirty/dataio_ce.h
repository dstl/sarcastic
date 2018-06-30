/** @file********************************************************************
 *
 *       Module:    dataio_ce.h
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 06/01/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions for reading from CaseExec files
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

#ifndef sarclib_DATAIO_CE_H__
#define sarclib_DATAIO_CE_H__

#include "sarclib.h"

/// Structure that contains infomation about a CaseExec dataset
///
typedef struct
{
    int do_swap;
    char * sd_name;
    char * ci_name;
    long nx;
    long ny;
    double sx;
    double sy;
} SPCEinfo;

/// Reads in the metadata from a case exec dataset
///
SPStatus* im_info_CEdataset(SPCEinfo * ce_info, const char * dir, SPStatus * status); 

/// Read in a complete case exec dataset from the given directory
///
SPStatus* im_load_CEdataset(SPImage * im, const char * dir, SPStatus * status);

/// reads in a im->nx, im->ny subset of image from offset ox, oy
/// using the already filled in ce_info structure
///
SPStatus* im_load_CEdataset_subset(SPImage *im, SPCEinfo * ce_info, long ox, long oy, SPStatus * status);

/// Clean up the SPCEinfo structure
///
SPStatus* im_destroy_SPCEinfo(SPCEinfo * info, SPStatus * status);

#endif


