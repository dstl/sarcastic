/***************************************************************************
 *
 *       Module:    image_extract.c
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 30/06/2005.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This function extracts a smaller part of an image from a larger original
 *      image.
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

#include "sarclib.h"

// This is an "insert loop" copies data from one image to another in the function im_insert
//
#define INS_LOOP(type) {                                                                                \
    for(curr_y = y; curr_y < (y + src->ny); curr_y++) {                                                 \
        memcpy(&dest->data. type [x + curr_y * dest->nx], &src->data. type [src->nx * (curr_y - y)],    \
        src->nx * im_getsizeoftype(src->image_type));                                                   \
    }                                                                                                   \
}

// This is an "extract loop" copies data from one image to another in the function im_extract
//
#define FLOOP(type) {                                                                                   \
    for(y = starty; y < starty+dest->ny; y++) {                                                         \
        memcpy(&dest->data. type [dest->nx * (y - starty)], &src->data. type [startx + src->nx * y],    \
        dest->nx * im_getsizeoftype(src->image_type));                                                  \
    }                                                                                                   \
}

// This shifts data in the x-direction
//
#define XSHIFT(type) {                                                                                  \
    for (y = 0; y<a->ny; y++) {                                                                         \
        memcpy(&(temp.data. type [y*a->nx+shift_x]), &(a->data. type [y*a->nx]),                        \
        im_getsizeoftype(a->image_type)*(a->nx-shift_x));                                               \
        memcpy(&(temp.data. type [y*a->nx]), &(a->data. type [y*a->nx+(a->nx-shift_x)]),                \
        im_getsizeoftype(a->image_type)*shift_x);                                                       \
    }                                                                                                   \
    im_destroy(a, status);                                                                              \
    im_clone(&temp, a, status);                                                                         \
}

// This shifts data in the y-direction
//
#define YSHIFT(type) {                                                                                  \
    if (shift_x != 0) {                                                                                 \
        im_init(&temp, status);                                                                         \
        im_create (&temp, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);                       \
    }                                                                                                   \
    memcpy(&(temp.data. type [shift_y*a->nx]), &(a->data. type [0]),                                    \
    im_getsizeoftype(a->image_type)*(a->ny-shift_y)*a->nx);                                             \
    memcpy(&(temp.data. type[0]), &(a->data. type[(a->ny-shift_y)*a->nx]),                              \
    im_getsizeoftype(a->image_type)*a->nx*shift_y);                                                     \
    im_destroy(a, status);                                                                              \
    im_clone(&temp, a, status);                                                                         \
}

// This shifts data
//
#define SHIFT(type) {                                                                                   \
    if (shift_x != 0)    /* for a shift in the x-direction */                                           \
    {                                                                                                   \
        XSHIFT(type);                                                                                   \
    }                                                                                                   \
    if (shift_y != 0)   /* for a shift in the y-direction */                                            \
    {                                                                                                   \
        YSHIFT(type);                                                                                   \
    }                                                                                                   \
    break;                                                                                              \
}

// Extracts a sub-image from larger original image. The sub-image will be dest.(nx, ny)
// in size and its (0,0) pixel will be at startx,starty.
//
SPStatus * im_extract (SPImage *src, int64_t startx, int64_t starty, SPImage * dest, SPStatus *status)
{
    int64_t y;
    CHECK_STATUS(status);
    CHECK_IMINTEG(src, status);
    
    if (startx < 0 || starty < 0)
    {
        fprintf(stderr, "im_extract: start position negative! (startx = %ld, starty = %ld)\n", (long)startx, (long)starty);
        status->status = INPUT_SIZES_MISMATCHED;
        return (status);
    }
    
    if (dest->nx + startx > src->nx)
    {
        fprintf(stderr, "im_extract: Requested image bigger than original image! (src->nx = %ld, x = %ld, dest->nx = %ld)\n",
                (long)src->nx, (long)startx, (long)dest->nx);
        status->status = INPUT_SIZES_MISMATCHED;
        return (status);
    }
    
    if (dest->ny + starty > src->ny)
    {
        fprintf(stderr, "im_extract: Requested image bigger than original image! (src->ny = %ld, y = %ld, dest->ny = %ld)\n",
                (long)src->ny, (long)starty, (long)dest->ny);
        status->status = INPUT_SIZES_MISMATCHED;
        return (status);
    }
    
    if (dest->data.v == NULL)
    {
        im_create(dest, src->image_type, dest->nx, dest->ny, src->xspc, src->yspc, status);
    }
    
    CHECK_TYPES_EQUAL(src, dest, status);
    
    switch (dest->image_type)
    {
        case ITYPE_POLAR:
            FLOOP(pol);
            break;
            
        case ITYPE_CMPL_FLOAT:
            FLOOP(cmpl_f);
            break;
            
        case ITYPE_VECTOR:
            FLOOP(vect);
            break;
            
        case ITYPE_DOUBLE:
            FLOOP(d);
            break;
            
        case ITYPE_FLOAT:
            FLOOP(f);
            break;
            
        case ITYPE_INT64:
            FLOOP(i64);
            break;
            
        case ITYPE_INT32:
            FLOOP(i32);
            break;
            
        case ITYPE_INT16:
            FLOOP(i16);
            break;
            
        case ITYPE_INT8:
            FLOOP(i8);
            break;
            
        case ITYPE_UINT64:
            FLOOP(ui64);
            break;
            
        case ITYPE_UINT32:
            FLOOP(ui32);
            break;
            
        case ITYPE_UINT16:
            FLOOP(ui16);
            break;
            
        case ITYPE_UINT8:
            FLOOP(ui8);
            break;
            
        case ITYPE_CMPL_INT64:
            FLOOP(cmpl_i64);
            break;
            
        case ITYPE_CMPL_INT32:
            FLOOP(cmpl_i32);
            break;
            
        case ITYPE_CMPL_INT16:
            FLOOP(cmpl_i16);
            break;
            
        case ITYPE_CMPL_INT8:
            FLOOP(cmpl_i8);
            break;
            
        default:
            fprintf(stderr, "Image type (%d) not supported in im_extract\n", dest->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    return (status);
}

// This shifts the data in the image by the given number of shift_x and shift_y
// pixels.
//
SPStatus* im_circshift (SPImage *a, int64_t shift_x, int64_t shift_y, SPStatus *status)
{
    int64_t y;
    SPImage temp;
    
    CHECK_STATUS(status);
    CHECK_PTR(a,status);
    CHECK_IMCONT(a,status);
    CHECK_IMINTEG(a,status);
    
    shift_x = (shift_x % a->nx);                        // for shifts larger than the array work out the resulting shift
    shift_y = (shift_y % a->ny);
    
    if ((shift_x != 0) || (shift_y != 0))
    {
        im_create (&temp, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
        
        CHECK_IMSIZEOUT(a, &temp, status);
        
        if (shift_x<0)                                  // calculating equivalent positive shift for negative shifts 
        {
            shift_x = a->nx+shift_x;
        }
        if (shift_y<0)
        {
            shift_y = a->ny+shift_y;
        }
        
        switch (a->image_type)
        {
            case ITYPE_CMPL_FLOAT:
                SHIFT(cmpl_f);
            case ITYPE_CMPL_INT64:
                SHIFT(cmpl_i64);
            case ITYPE_CMPL_INT32:
                SHIFT(cmpl_i32);
            case ITYPE_CMPL_INT16:
                SHIFT(cmpl_i16);
            case ITYPE_CMPL_INT8:
                SHIFT(cmpl_i8);
            case ITYPE_POLAR:
                SHIFT(pol);
            case ITYPE_DOUBLE:
                SHIFT(d);
            case ITYPE_FLOAT:
                SHIFT(f);
            case ITYPE_INT64:
                SHIFT(i64);
            case ITYPE_INT32:
                SHIFT(i32);
            case ITYPE_INT16:
                SHIFT(i16);
            case ITYPE_INT8:
                SHIFT(i8);
            case ITYPE_UINT64:
                SHIFT(ui64);
            case ITYPE_UINT32:
                SHIFT(ui32);
            case ITYPE_UINT16:
                SHIFT(ui16);
            case ITYPE_UINT8:
                SHIFT(ui8);
            case ITYPE_VECTOR:
                SHIFT(vect);
            default:
                fprintf(stderr, "Src image type (%d) not supported in im_circshift\n", a->image_type);
                status->status = INVALID_TYPE;
                break;
        }
    }
    
    return(status);
}

// This inserts the given src image into the dest image at position x,y.
// The data is overwritten in the dest image.
//
SPStatus* im_insert(SPImage * src, int64_t x, int64_t y, SPImage * dest, SPStatus * status)
{
    int64_t curr_y;
    
    CHECK_STATUS(status);
    CHECK_PTR(src,status);
    CHECK_IMCONT(src,status);
    CHECK_IMINTEG(src,status);
    
    CHECK_PTR(dest,status);
    CHECK_IMCONT(dest,status);
    CHECK_IMINTEG(dest,status);
    
    if (src->nx + x > dest->nx)
    {
        fprintf(stderr, "Image won't fit in destination image! (im_insert) (src->nx = %ld, x = %ld, dest->nx = %ld)\n",
                (long)src->nx, (long)x, (long)dest->nx);
        status->status = INPUT_SIZES_MISMATCHED;
        return (status);
    }
    
    if (src->ny + y > dest->ny)
    {
        fprintf(stderr, "Image won't fit in destination image! (im_insert) (src->ny = %ld, y = %ld, dest->ny = %ld)\n",
                (long)src->ny, (long)y, (long)dest->ny);
        
        status->status = INPUT_SIZES_MISMATCHED;
        return (status);
    }
    
    CHECK_TYPES_EQUAL(src, dest, status);
    
    switch (src->image_type)
    {
        case ITYPE_POLAR:
            INS_LOOP(pol);
            break;
            
        case ITYPE_CMPL_FLOAT:
            INS_LOOP(cmpl_f);
            break;
            
        case ITYPE_VECTOR:
            INS_LOOP(vect);
            break;
            
        case ITYPE_DOUBLE:
            INS_LOOP(d);
            break;
            
        case ITYPE_FLOAT:
            INS_LOOP(f);
            break;
            
        case ITYPE_INT64:
            INS_LOOP(i64);
            break;
            
        case ITYPE_INT32:
            INS_LOOP(i32);
            break;
            
        case ITYPE_INT16:
            INS_LOOP(i16);
            break;
            
        case ITYPE_INT8:
            INS_LOOP(i8);
            break;
            
        case ITYPE_UINT64:
            INS_LOOP(ui64);
            break;
            
        case ITYPE_UINT32:
            INS_LOOP(ui32);
            break;
            
        case ITYPE_UINT16:
            INS_LOOP(ui16);
            break;
            
        case ITYPE_UINT8:
            INS_LOOP(ui8);
            break;
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_insert\n", src->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    return (status);
}
