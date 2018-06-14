/** @file********************************************************************
 *
 *       Module:    dataio_im.h
 *      Program:    SIlib2
 *   Created by:    Emma Griffiths on 19/07/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains the prottypes for the functions that read and write
 *      our file format.
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

#ifndef SILib2_DATAIO_IM_H__
#define SILib2_DATAIO_IM_H__

#include "SIlib2.h"

/// Structure that contains info about the image being saved a line at a time
///
typedef struct
{
    int64_t nx;         ///< nx dimension of image being saved
    char * fname;       ///< filename that is being written to
    FILE * fp;          ///< file pointer for file being written to
} SPImageLineSaveInfo;  

/// Saves the header for an image
///
/// The header conists of:
///
///  - a uint32  magic number 0x87654321
///  - a uint64 pair of numbers describing the size of the image (x then y axis)
///  - a pair of doubles describing the pixels spacing of each axis (x then y axis)
///  - an int32 describing the image type (eg. complex float = 2, double = 4,
///
///  look in image.h for the declaration of SPImageType for the values of all the
///     other types
///
SPStatus * im_save_header(SPImage * a, FILE * fp, SPStatus *status);


/// Save a complete image to a file
///
SPStatus * im_save (SPImage *a, const char *fname, SPStatus *status);


/// Save a complete image to a file but also specifying the
/// corner coordinates of the image. Corner coordinates are specified
/// as ECEF vectors of teh corners of the image starting with the top left
/// and going clockwise around the image (tl,tr,br,bl)
///
SPStatus * im_save_with_geo_ecf(SPImage *a, SPVector tl, SPVector tr, SPVector br, SPVector bl, const char *fname, SPStatus *status);


/// Initialises the line save info structure so that a file can be
/// written a few lines at a time
///
SPStatus * im_save_line_init (SPImage *a, const char *fname, SPImageLineSaveInfo * info, SPStatus *status);


/// Saves an image to a file that is already partially written. This is useful for writing
/// an image to a file a piece at a time. The output file must be already initialised with
/// im_save_line_init().
///
SPStatus * im_save_line (SPImage *a, SPImageLineSaveInfo * info, SPStatus *status);


/// Closes down the line save info structure and file
///
SPStatus * im_save_line_close (SPImageLineSaveInfo * info, SPStatus *status);


/// Closes down the line save info structure and file
///
SPStatus * im_save_line_close_with_geo_ecf (SPImageLineSaveInfo * info, SPVector p1, SPVector p2, SPVector p3, SPVector p4, SPStatus *status);


/// Loads up a complete image from file
///
SPStatus * im_load (SPImage *a, const char *fname, SPStatus *status);    /* Loads up a complete image from file */


/// Loads up just the metadata from file
///
SPStatus * im_load_metadata (SPImage *a, const char *fname, SPStatus *status);  /* Loads up just the metadata from file */


/// Loads up a subset of an image from a file
///
SPStatus * im_load_subset (SPImage *a, const char *fname, int64_t ox, int64_t oy, SPStatus *status);   /* Loads up a subset of an image from a file */


/// Perform byte swapping on the data within an image
///
void swap_data(SPImage * a, int do_swap);

#endif
