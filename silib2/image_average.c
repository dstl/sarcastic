/***************************************************************************
 *
 *       Module:    image_average.c
 *      Program:    SIlib2
 *   Created by:    Matt Nottingham on 30/06/2005.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Averages down an image into a smaller image (coherently).
 *
 *      (Lots of fun with macros.....if only you could for a foreach with the
 *      C pre-processor!)
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

#define AVE(zero, accum, divide) {                          \
    for(y = 0; y < src->ny; y += step_y) {                  \
        for(x = 0; x < src->nx; x += step_x) {              \
            zero ;                                          \
            for(suby = 0; suby < step_y; suby++) {          \
                for(subx = 0; subx < step_x; subx++) {      \
                    accum ;                                 \
                }                                           \
            }                                               \
            divide ;                                        \
        }                                                   \
    }                                                       \
}

#define S_AVE(type)                                                                        \
    AVE(total = 0, total += src->data.type [(x + subx) + (y + suby) * src->nx],            \
    dest->data.type [(x / step_x) + dest->nx * y / step_y] = total / (step_x * step_y) )


// Averages down the src image to be the dest image size.
//
SPStatus* im_average(SPImage * src, SPImage * dest, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_PTR(dest,status);
    CHECK_IMINTEG(src,status);
    int64_t x;
    int64_t y;
    int64_t step_x;
    int64_t step_y;
    int64_t suby;
    int64_t subx;
    double total;
    
    if (dest->nx <= 0 || dest->ny <= 0)
    {
        fprintf(stderr, "im_average: Output image not allowable size (nx = %ld, ny = %ld)\n", (long)dest->nx, (long)dest->ny);
        status->status = OUTPUT_SIZES_MISMATCHED;
        return (status);
    }
    
    if ((src->nx % dest->nx) != 0 || (src->ny % dest->ny) != 0)
    {
        fprintf(stderr, "im_average: Output image not allowable size (src(nx = %ld, ny = %ld), dest(nx = %ld, ny = %ld))\n", (long)src->nx, (long)src->ny, (long)dest->nx, (long)dest->ny);
        status->status = OUTPUT_SIZES_MISMATCHED;
        return (status);
    }
    
    step_x = src->nx / dest->nx;
    step_y = src->ny / dest->ny;
    
    if (dest->data.v == NULL)
    {
        im_create(dest, src->image_type, dest->nx, dest->ny, src->xspc * step_x, src->yspc * step_y, status);
    }
    else
    {
        dest->xspc = src->xspc * step_x;
        dest->yspc = src->yspc * step_y;
    }
    
    CHECK_STATUS(status);
    
    switch(src->image_type)
    {
        case ITYPE_POLAR:
            AVE({dest->data.pol[(x / step_x) + dest->nx * y / step_y].pabs = 0.0;dest->data.pol[(x / step_x) + dest->nx * y / step_y].parg = 0.0;},
                CMPLX_PADD(src->data.pol[(x + subx) + (y + suby) * src->nx], dest->data.pol[(x / step_x) + dest->nx * y / step_y],dest->data.pol[(x / step_x) + dest->nx * y / step_y] ),
                dest->data.pol[(x / step_x) + dest->nx * y / step_y].pabs /= (step_x * step_y); );
            break;
            
        case ITYPE_CMPL_FLOAT:
            AVE({dest->data.cmpl_f[(x / step_x) + dest->nx * y / step_y].r = 0.0;dest->data.cmpl_f[(x / step_x) + dest->nx * y / step_y].i = 0.0;},
                CMPLX_ADD(src->data.cmpl_f[(x + subx) + (y + suby) * src->nx], dest->data.cmpl_f[(x / step_x) + dest->nx * y / step_y],dest->data.cmpl_f[(x / step_x) + dest->nx * y / step_y] ),
                {dest->data.cmpl_f[(x / step_x) + dest->nx * y / step_y].r /= (step_x * step_y); dest->data.cmpl_f[(x / step_x) + dest->nx * y / step_y].i /= (step_x * step_y);} );
            break;
            
        case ITYPE_VECTOR:
            AVE({dest->data.vect[(x / step_x) + dest->nx * y / step_y].x = 0.0;dest->data.vect[(x / step_x) + dest->nx * y / step_y].y = 0.0;dest->data.vect[(x / step_x) + dest->nx * y / step_y].z = 0.0;},
                VECT_ADD(src->data.vect[(x + subx) + (y + suby) * src->nx], dest->data.vect[(x / step_x) + dest->nx * y / step_y],dest->data.vect[(x / step_x) + dest->nx * y / step_y] ),
                {dest->data.vect[(x / step_x) + dest->nx * y / step_y].x /= (step_x * step_y); dest->data.vect[(x / step_x) + dest->nx * y / step_y].y /= (step_x * step_y);dest->data.vect[(x / step_x) + dest->nx * y / step_y].z /= (step_x * step_y);} );
            
            break;
            
        case ITYPE_DOUBLE:
            S_AVE(d);
            break;
            
        case ITYPE_FLOAT:
            S_AVE(f);
            break;
            
        case ITYPE_INT64:
            S_AVE(i64);
            break;
            
        case ITYPE_INT32:
            S_AVE(i32);
            break;
            
        case ITYPE_INT16:
            S_AVE(i16);
            break;
            
        case ITYPE_INT8:
            S_AVE(i8);
            break;
            
        case ITYPE_UINT64:
            S_AVE(ui64);
            break;
            
        case ITYPE_UINT32:
            S_AVE(ui32);
            break;
            
        case ITYPE_UINT16:
            S_AVE(ui16);
            break;
            
        case ITYPE_UINT8:
            S_AVE(ui8);
            break;
            
        default:
            fprintf(stderr, "im_average: Src image type (%d) not supported\n", src->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    
    return (status);
}
