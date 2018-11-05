/***************************************************************************
 * 
 *           Module :  image_transpose.c
 *          Program :  sarclib
 *       Created by :  Emma Griffiths on 30/06/2005
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Function to transpose an image
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

#include "sarclib.h"

// This macro simply transposes data
//
#define FLOOP(type) {                                                       \
    for (y=0; y<a->ny; y++) {                                               \
        for (x=0; x<a->nx; x++) {                                           \
            tmp.data. type [x * a->ny + y] = a->data. type [y * a->nx + x]; \
        }                                                                   \
    }                                                                       \
}


// Transposes the supplied image.
//
SPStatus* im_transp (SPImage *a, SPStatus *status)
{
    int64_t x;
    int64_t y;
    SPImage tmp;
    CHECK_STATUS(status);
    CHECK_IMINTEG(a, status);
    
    im_create(&tmp, a->image_type, a->ny, a->nx, a->yspc, a->xspc, status);
    
    switch (a->image_type)
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
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_transp\n", a->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    im_destroy(a, status);
    im_clone(&tmp, a, status);
    
    return (status);
}

