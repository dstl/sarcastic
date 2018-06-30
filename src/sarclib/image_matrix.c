/***************************************************************************
 *
 *       Module:    image_matrix.c
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 22/03/2007.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains functions needed to operate on complex polar arrays.
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

#define MMULT(type) {                                                                                                           \
    for(y = 0; y < out->ny; y++) {                                                                                              \
        for(x = 0; x < out->nx; x++) {                                                                                          \
            out->data. type [x + y * out->nx] = 0;                                                                              \
            for (n = 0; n < in1->nx; n++)  {                                                                                    \
                out->data. type [x + y * out->nx] +=  in1->data. type [n + y * in1->nx] * in2->data. type [x + n * in2->nx];    \
            }                                                                                                                   \
        }                                                                                                                       \
    }                                                                                                                           \
}

#define ALLTYPES()      \
    case ITYPE_DOUBLE:  \
        MMULT(d);       \
        break;          \
                        \
    case ITYPE_FLOAT:   \
        MMULT(f);       \
        break;          \
                        \
    case ITYPE_INT64:   \
        MMULT(i64);     \
        break;          \
                        \
    case ITYPE_INT32:   \
        MMULT(i32);     \
        break;          \
                        \
    case ITYPE_INT16:   \
        MMULT(i16);     \
        break;          \
                        \
    case ITYPE_INT8:    \
        MMULT(i8);      \
        break;          \
                        \
    case ITYPE_UINT64:  \
        MMULT(ui64);    \
        break;          \
                        \
    case ITYPE_UINT32:  \
        MMULT(ui32);    \
        break;          \
                        \
    case ITYPE_UINT16:  \
        MMULT(ui16);    \
        break;          \
                        \
    case ITYPE_UINT8:   \
        MMULT(ui8);     \
        break;          \


// Multiplies two images and outputs to a third - using matrix style multiplication.
//
SPStatus* im_matrix_mult (SPImage *in1, SPImage *in2, SPImage *out, SPStatus *status)
{
    int64_t n;
    int64_t x;
    int64_t y;
    SPCmplx tmp;
    
    CHECK_STATUS(status);
    CHECK_PTR(out,status);
    CHECK_IMINTEG(in1, status);
    CHECK_IMINTEG(in2, status);
    CHECK_TYPES_EQUAL(in1, in2, status);
    if (in1->nx  != in2->ny) {
        printf("Matrices not of allowable size (%"PRId64", %"PRId64") x (%"PRId64", %"PRId64") (row, col)\n",
               in1->ny, in1->nx, in2->ny, in2->nx);
        status->status = INPUT_SIZES_MISMATCHED;
        return status;
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0)) {
        im_create (out, in1->image_type, in2->nx, in1->ny, in2->xspc, in1->yspc, status);
        CHECK_STATUS(status);
    }
    
    CHECK_TYPES_EQUAL(in1, out, status);
    
    if (out->nx != in2->nx || out->ny != in1->ny) {
        printf("Output matrix not of allowable size (%"PRId64", %"PRId64") x (%"PRId64", %"PRId64") = (%"PRId64", %"PRId64")(row, col)\n   However output matrix is (%"PRId64", %"PRId64")\n\n",
               in1->ny, in1->nx, in2->ny, in2->nx, in1->ny, in2->nx, out->ny, out->nx);
        status->status = OUTPUT_SIZES_MISMATCHED;
        return status;
    }
    
    switch (in1->image_type)
    {
            
        case ITYPE_CMPL_FLOAT:
            for(y = 0; y < out->ny; y++) {
                for(x = 0; x < out->nx; x++) {
                    out->data.cmpl_f[x + y * out->nx].r = 0.0;
                    out->data.cmpl_f[x + y * out->nx].i = 0.0;
                    for (n=0; n < in1->nx; n++)
                    {
                        CMPLX_MULT(in1->data.cmpl_f[n + y * in1->nx], in2->data.cmpl_f[x + n * in2->nx], tmp);
                        CMPLX_ADD(tmp, out->data.cmpl_f[x + y * out->nx], out->data.cmpl_f[x + y * out->nx])
                    }
                }
            }
            break;
            
            ALLTYPES();
            
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_mult\n", in1->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    return(status);
}

