/** @file********************************************************************
 *
 *       Module:    dataio_sicd.h
 *      Program:    sarclib
 *   Created by:    Matt Nottingham in 2011.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions needed for reading and writing the
 *      SICD (Sensor Independent Complex Data) format.
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

#ifndef sarclib_DATAIO_SICD_H__
#define sarclib_DATAIO_SICD_H__

#include "sarclib.h"

typedef struct {
  int numi;

} SICDmetadata;


/// Save a complete image to a file
///
SPStatus * im_save_sicd (SPImage *a, const char *fname, SICDmetadata * mdata, SPStatus *status);


/// Initialises the line save info structure so that a file can be
/// written a few lines at a time
///
SPStatus * im_save_sicd_line_init(SPImage * im, const char * fname, SPImageLineSaveInfo * info, SICDmetadata * mdata, SPStatus * status);


/// Saves an image to a file that is already partially written. This is useful for writing
/// an image to a file a piece at a time. The output file must be already initialised with
/// im_save_sicd_line_init().
///
SPStatus * im_save_sicd_line(SPImage * im, SPImageLineSaveInfo * info, SPStatus * status);


/// Closes down the line save info structure and file
///
SPStatus * im_save_sicd_line_close(SPImageLineSaveInfo * info, SPStatus * status);


/// Save just the header information - currently just a non-working stub
///
SPStatus * im_save_sicd_header(SPImage *a, SICDmetadata * mdata, FILE * fp, SPStatus * status);


/// Loads up just the metadata from file
///
SPStatus * im_load_sicd_metadata (SPImage *a, const char *fname, SICDmetadata * mdata, SPStatus *status); 


/// Loads up a complete image from file
///
SPStatus * im_load_sicd (SPImage *a, const char *fname, SICDmetadata * mdata, SPStatus *status); 


/// Loads up a subset of an image from a file
///
SPStatus * im_load_sicd_subset (SPImage *a, const char *fname, int64_t ox, int64_t oy, SPStatus *status);  /* Loads up a subset of an image from a file */

#endif
