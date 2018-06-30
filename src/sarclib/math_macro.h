/** @file********************************************************************
 *
 *       Module:    math_macro.h
 *      Program:    sarclib
 *   Created by:    Matt Nottingham on 30/06/2005.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *   This file contains some wonderful macro definitions! Use with caution!
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


#ifndef sarclib_MATH_MACRO__
#define sarclib_MATH_MACRO__

/// the macro below performs a given operation on every element of an
/// array of a given type where the input and output location are the
/// same
///
#define FOREACH(sub, oper) {                                            \
    for (n=0; n<((a->nx)*(a->ny)); n++) {                               \
        out->data. sub [n] = a->data. sub [n] oper b->data. sub [n];  	\
    }                                                               	\
}

/// This macro sets the specified element of the specified data type
/// to the specified value
///
#define SETELEMENT(sub, pos, val) {     \
    out->data. sub [pos] = val;         \
}

/// his macro performs a specified operation on two specified input elements of a
/// given type and outputs it to a partiular element
///
#define FORELEMENT(sub, oper, in1, posin1, in2, posin2, out, posout) {              \
    out->data. sub [posout] = in1->data. sub [posin1] oper in2->data. sub [posin2]; \
}

/// This macro performs the "FOREACH" macro to perform an operation whichever
/// of all the currently used image types the input and output images are
///
#define FORALLTYPES(oper)       \
    case ITYPE_DOUBLE:          \
        FOREACH(d, oper);       \
        break;                  \
                                \
    case ITYPE_FLOAT:           \
        FOREACH(f, oper);       \
        break;                  \
                                \
    case ITYPE_INT64:           \
        FOREACH(i64, oper);     \
        break;                  \
                                \
    case ITYPE_INT32:           \
        FOREACH(i32, oper);     \
        break;                  \
                                \
    case ITYPE_INT16:           \
        FOREACH(i16, oper);     \
        break;                  \
                                \
    case ITYPE_INT8:            \
        FOREACH(i8, oper);      \
        break;                  \
                                \
    case ITYPE_UINT64:          \
        FOREACH(ui64, oper);    \
        break;                  \
                                \
    case ITYPE_UINT32:          \
        FOREACH(ui32, oper);    \
        break;                  \
                                \
    case ITYPE_UINT16:          \
        FOREACH(ui16, oper);    \
        break;                  \
                                \
    case ITYPE_UINT8:           \
        FOREACH(ui8, oper);     \
        break;                  \

/// This macro uses the "SETELEMENT" macro to set each element of whichever
/// of each possible image type the input and output images are
///
#define SETALLTYPES(pos, val)       \
    case ITYPE_DOUBLE:              \
        SETELEMENT(d, pos, val);    \
        break;                      \
                                    \
    case ITYPE_FLOAT:               \
        SETELEMENT(f, pos, val);    \
        break;                      \
                                    \
    case ITYPE_INT64:               \
        SETELEMENT(i64, pos, val);  \
        break;                      \
                                    \
    case ITYPE_INT32:               \
        SETELEMENT(i32, pos, val);  \
        break;                      \
                                    \
    case ITYPE_INT16:               \
        SETELEMENT(i16, pos, val);  \
        break;                      \
                                    \
    case ITYPE_INT8:                \
        SETELEMENT(i8, pos, val);   \
        break;                      \
                                    \
    case ITYPE_UINT64:              \
        SETELEMENT(ui64, pos, val); \
        break;                      \
                                    \
    case ITYPE_UINT32:              \
        SETELEMENT(ui32, pos, val); \
        break;                      \
                                    \
    case ITYPE_UINT16:              \
        SETELEMENT(ui16, pos, val); \
        break;                      \
                                    \
    case ITYPE_UINT8:               \
        SETELEMENT(ui8, pos, val);  \
    break;                          \

/// this macro performs an operation on any image type where the input
/// and output indeces may not all be the same */
///
#define OPERALLTYPESX(oper, in1, in2, out)                              \
    case ITYPE_DOUBLE:                                                  \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(d, oper, in1, x+y*in1->nx, in2, y, out, y);      \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_FLOAT:                                                   \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(f, oper, in1, x+y*in1->nx, in2, y, out, y);      \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_INT64:                                                   \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(i64, oper, in1, x+y*in1->nx, in2, y, out, y);    \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_INT32:                                                   \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(i32, oper, in1, x+y*in1->nx, in2, y, out, y);    \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_INT16:                                                   \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(i16, oper, in1, x+y*in1->nx, in2, y, out, y);    \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_INT8:                                                    \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(i8, oper, in1, x+y*in1->nx, in2, y, out, y);     \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_UINT64:                                                  \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(ui64, oper, in1, x+y*in1->nx, in2, y, out, y);   \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_UINT32:                                                  \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(ui32, oper, in1, x+y*in1->nx, in2, y, out, y);   \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_UINT16:                                                  \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(ui16, oper, in1, x+y*in1->nx, in2, y, out, y);   \
        }                                                               \
        break;                                                          \
                                                                        \
    case ITYPE_UINT8:                                                   \
        for (x = 0; x < a->nx; x++)                                     \
        {                                                               \
            FORELEMENT(ui16, oper, in1, x+y*in1->nx, in2, y, out, y);   \
        }                                                               \
        break;                                                          \


/// This macro performs a specified operation on the rows of an image
///
#define FORROWS(a, out, oper_type) {                                                                            \
    int x = 0, y;                                                                                               \
    for (y = 0; y < a->ny; y++)                                                                                 \
    {                                                                                                           \
        switch(a->image_type) {                                                                                 \
            case ITYPE_POLAR:                                                                                   \
                out->data.pol[y].pabs = 0;                                                                      \
                out->data.pol[y].parg = 0;                                                                      \
                break;                                                                                          \
            case ITYPE_CMPL_FLOAT:                                                                              \
                out->data.cmpl_f[y].r = 0;                                                                      \
                out->data.cmpl_f[y].i = 0;                                                                      \
                break;                                                                                          \
            case ITYPE_VECTOR:                                                                                  \
                out->data.vect[y].x = 0;                                                                        \
                out->data.vect[y].y = 0;                                                                        \
                out->data.vect[y].z = 0;                                                                        \
                break;                                                                                          \
            SETALLTYPES(y, 0);                                                                                  \
            default:                                                                                            \
                fprintf(stderr, "Src image type (%d) not supported in FORROWS\n", a->image_type);               \
                status->status = INVALID_TYPE;                                                                  \
                break;                                                                                          \
            }                                                                                                   \
        switch(oper_type)                                                                                       \
        {                                                                                                       \
            case OPER_ADD:                                                                                      \
                switch(a->image_type)                                                                           \
                {                                                                                               \
                    case ITYPE_POLAR:                                                                           \
                        for (x = 0; x < a->nx; x++)                                                             \
                        {                                                                                       \
                            CMPLX_PADD(a->data.pol[x+y*a->nx], out->data.pol[y], out->data.pol[y]);             \
                        }                                                                                       \
                        break;                                                                                  \
                                                                                                                \
                    case ITYPE_CMPL_FLOAT:                                                                      \
                        for (x = 0; x < a->nx; x++)                                                             \
                        {                                                                                       \
                            CMPLX_ADD(a->data.cmpl_f[x+y*a->nx], out->data.cmpl_f[y], out->data.cmpl_f[y]);     \
                        }                                                                                       \
                        break;                                                                                  \
                                                                                                                \
                    case ITYPE_VECTOR:                                                                          \
                        for (x = 0; x < a->nx; x++)                                                             \
                        {                                                                                       \
                            VECT_ADD(a->data.vect[x+y*a->nx], out->data.vect[y], out->data.vect[y]);  		    \
                        }                                                                                       \
                        break;                                                                                  \
                                                                                                                \
                        OPERALLTYPESX(+, a, out, out);                                                          \
                                                                                                                \
                        break;                                                                                  \
                                                                                                                \
                    default:                                                                                    \
                        fprintf(stderr, "Src image type (%d) not supported in FORROWS\n", a->image_type);       \
                        status->status = INVALID_TYPE;                                                          \
                        break;                                                                                  \
                }                                                                                               \
                                                                                                                \
                break;                                                                                          \
                                                                                                                \
            default:                                                                                            \
                fprintf(stderr, "Operator type %d not supperted in FORROWS macro\n", oper_type);                \
                status->status = UNSUPPORTED_OPERATOR;                                                          \
                break;                                                                                          \
        }                                                                                                       \
    }                                                                                                           \
}

typedef enum {OPER_ADD, OPER_SUB, OPER_MULT, OPER_DIV, OPER_SET} SPOperator;


/// Some simple macros to convert between degrees and radians
///
#define DEG2RAD(a)   ((a) * M_PI / 180.0)
#define RAD2DEG(a)   ((a) * 180.0 / M_PI)

#endif
