/***************************************************************************
 *
 *       Module:    image_unwrap.c
 *      Program:    SIlib2
 *   Created by:    Matt Nottingham on 23/03/2007.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This function takes a complex or floating point line of phase, then unwraps
 *      it and returns a floating point line.
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

#define sign(a)  (((a) >= 0) ? 1 : 0)

// Unwrap phase from the supplied 1D array
//
SPStatus * im_unwrap (SPImage *phas, SPImage * unwrapped, SPStatus *status)
{
    int64_t i;
    SPImage actual_phase;
    
    CHECK_STATUS(status);
    CHECK_PTR(phas,status);
    CHECK_IMCONT(phas,status);
    CHECK_IMINTEG(phas,status);
    
    if (phas->nx != 1 && phas->ny != 1) {
        status->status = ARRAY_SIZES_WRONG;
        fprintf(stderr, "im_unwrap: Expected a 1D array, but got %ld x %ld\n", (long) phas->nx, (long)phas->ny);
        return status;
    }
    
    switch (phas->image_type)
    {
        case ITYPE_CMPL_FLOAT:
            im_create(&actual_phase, ITYPE_FLOAT, phas->nx, phas->ny, phas->xspc, phas->yspc, status);
            im_cast(phas, CAST_PHASE, &actual_phase, status);
            CHECK_STATUS(status);
            break;
            
        case ITYPE_FLOAT:
            im_clone(phas, &actual_phase, status);
            break;
            
        default:
            fprintf(stderr, "Image type (%d) not supported in im_unwrap\n", phas->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    if (unwrapped->nx == 0 && unwrapped->ny == 0 && unwrapped->data.v == NULL) {
        im_create(unwrapped, ITYPE_FLOAT, phas->nx, phas->ny, phas->xspc, phas->yspc, status);
    }
    
    unwrapped->data.f[0] = actual_phase.data.f[0];
    
    for(i = 1; i < (phas->nx * phas->ny); i++) {
        if (sign(actual_phase.data.f[i-1]) != sign(actual_phase.data.f[i]) &&  (fabs(actual_phase.data.f[i] - actual_phase.data.f[i-1]) > M_PI)) {
            if (actual_phase.data.f[i] < 0 && actual_phase.data.f[i-1] > 0) {
                unwrapped->data.f[i] = unwrapped->data.f[i-1] + (actual_phase.data.f[i] - actual_phase.data.f[i-1] + 2.0 * M_PI);
            } else {
                unwrapped->data.f[i] = unwrapped->data.f[i-1] + (actual_phase.data.f[i] - actual_phase.data.f[i-1] - 2.0 * M_PI);;
            }
        } else {
            unwrapped->data.f[i] = unwrapped->data.f[i-1] + (actual_phase.data.f[i] - actual_phase.data.f[i-1]);
        }
    }
    
    if (phas->image_type == ITYPE_CMPL_FLOAT) {
        im_destroy(&actual_phase, status);
    }
    
    return (status);
}

