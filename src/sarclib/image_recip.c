/***************************************************************************
 *
 *       Module:    image_recip.c
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 20/09/2006.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      The function below takes the reciprocal of an image.
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

#define RECIP(type){                                                    \
    for (i = 0; i <in->nx*in->ny; i++) {                                \
        out->data. type [i] = 1/in->data. type [i];                     \
    }                                                                   \
}

#define RECIP_CMPL_INTXX(no){                                           \
    case ITYPE_CMPL_INT##no :                                           \
        for (i = 0; i<in->nx*in->ny; i++) {                             \
            out->data.cmpl_i##no [i].r = 1.0/in->data.cmpl_i##no [i].r; \
            out->data.cmpl_i##no [i].i = 1.0/in->data.cmpl_i##no [i].i; \
        }                                                               \
        break;                                                          \
}

#define ALL_CASES(macro){                                               \
    case ITYPE_FLOAT:                                                   \
        macro(f);                                                       \
        break;                                                          \
    case ITYPE_DOUBLE:                                                  \
        macro(d);                                                       \
        break;                                                          \
    case ITYPE_INT64:                                                   \
        macro(i64);                                                     \
        break;                                                          \
    case ITYPE_INT32:                                                   \
        macro(i32);                                                     \
        break;                                                          \
    case ITYPE_INT16:                                                   \
        macro(i16);                                                     \
        break;                                                          \
    case ITYPE_INT8:                                                    \
        macro(i8);                                                      \
        break;                                                          \
    case ITYPE_UINT64:                                                  \
        macro(ui64);                                                    \
        break;                                                          \
    case ITYPE_UINT32:                                                  \
        macro(ui32);                                                    \
        break;                                                          \
    case ITYPE_UINT16:                                                  \
        macro(ui16);                                                    \
        break;                                                          \
    case ITYPE_UINT8:                                                   \
        macro(ui8);                                                     \
        break;                                                          \
}

// This is a function to take the reciprocal of an image, in and out can be the same
// image if necessary
//
SPStatus *im_recip(SPImage *in, SPImage *out, SPStatus *status)
{
    int64_t i;
    CHECK_STATUS(status);
    CHECK_IMINTEG(in, status);
    CHECK_PTR(out, status);
    CHECK_SPC(in, out, status);
    CHECK_IMSIZEIN(in, out, status);
    CHECK_TYPES_EQUAL(in, out, status);
    
    switch(in->image_type)
    {
        case ITYPE_POLAR:
            for (i = 0; i<in->nx*in->ny; i++)
            {
                out->data.pol[i].pabs = 1.0/in->data.pol[i].pabs;
                out->data.pol[i].parg = 1.0/in->data.pol[i].parg;
            }
            break;
        case ITYPE_CMPL_FLOAT:
            for (i = 0; i<in->nx*in->ny; i++)
            {
                out->data.cmpl_f[i].r = 1.0/in->data.cmpl_f[i].r;
                out->data.cmpl_f[i].i = 1.0/in->data.cmpl_f[i].i;
            }
            break;
            RECIP_CMPL_INTXX(64);
            RECIP_CMPL_INTXX(32);
            RECIP_CMPL_INTXX(16);
            RECIP_CMPL_INTXX(8);
        case ITYPE_VECTOR:
            for (i = 0; i<in->nx*in->ny; i++)
            {
                out->data.vect[i].x = 1.0/in->data.vect[i].x;
                out->data.vect[i].y = 1.0/in->data.vect[i].y;
                out->data.vect[i].z = 1.0/in->data.vect[i].z;
            }
            break;
            ALL_CASES(RECIP);
        case ITYPE_CMPL_NIBBLE:
            fprintf(stderr, "im_recip does not support complex nibbles\n");
            status->status = INVALID_TYPE;
            break;
        case ITYPE_UNKNOWN:
            fprintf(stderr, "im_recip in %s:%d expects a recognised ITYPE\n", __FILE__, __LINE__);
            status->status = INVALID_TYPE;
            break;
    }
    return(status);
}





