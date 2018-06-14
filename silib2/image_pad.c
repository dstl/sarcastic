/***************************************************************************
 *
 *       Module:    image_pad.c
 *      Program:    SIlib2
 *   Created by:    Emma Griffiths on 22/06/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains a function to pad an image.
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

#include "SIlib2.h"

// zero pad an image structure to a new size
//
SPStatus* im_zpad (SPImage *a, int64_t new_size_x, int64_t new_size_y, SPImage *out, SPStatus *status)
{
    CHECK_STATUS(status);
    CHECK_IMINTEG(a, status);
    
    if ( (new_size_x < a->nx) || (new_size_y < a->ny) )
    {
        status->status = ARRAY_SIZES_WRONG;
        fprintf(stderr, "%s:%d new x %ld should be > nx %ld or new y %ld should be > ny %ld\n", __FILE__, __LINE__, (long) new_size_x, (long) a->nx, (long) new_size_y, (long) a->ny);
        CHECK_STATUS(status);
    }
    
    // If out is a empty image, then create it to be the correct size otherwise make sure the image is zeroed
    //
    if (out->nx == 0 && out->ny == 0 && out->data.v == NULL)
    {
        im_create(out, a->image_type, new_size_x, new_size_y, a->xspc, a->yspc, status);
    }
    else
    {
        im_zero_data(out, status);
    }
    
    CHECK_IMINTEG(out, status);
    
    im_insert(a, 0, 0, out, status);  // inserts the old image from 0,0 into the new padded image
    
    CHECK_STATUS(status);
    
    return (status);
}
