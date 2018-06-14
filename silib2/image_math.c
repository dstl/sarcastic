/***************************************************************************
 *
 *       Module:    image_math.c
 *      Program:    SIlib2
 *   Created by:    Emma Griffiths on 22/03/2007.
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

#include "SIlib2.h"

// this defines a number of checks commonly used
//
#define CHECKS() {                                                                  \
    CHECK_STATUS(status);                                                           \
    CHECK_PTR(out,status);                                                          \
    CHECK_SPC(a,b,status);                                                          \
    CHECK_IMINTEG(a,status);                                                        \
    CHECK_IMINTEG(b,status);                                                        \
    CHECK_IMSIZEIN(a,b,status);                                                     \
    CHECK_TYPES_EQUAL(a, b, status);                                                \
                                                                                    \
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0))                   \
    {                                                                               \
        im_create (out, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);     \
    }                                                                               \
    CHECK_TYPES_EQUAL(a,out,status);                                                \
    CHECK_IMSIZEOUT(a,out,status);                                                  \
}

// this macro multiplies a polar array by a constant 
//
#define MULT_POL_CONST(type, con){                                                  \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        im->data.pol[i].pabs *= c;                                                  \
    }                                                                               \
}

// this macro multiplies a caratestian array by a constant 
//
#define MULT_CMPL_CONST(type, con){                                                 \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx * im->ny; i++) {                                         \
        im->data.cmpl_f[i].r *= c;                                                  \
        im->data.cmpl_f[i].i *= c;                                                  \
    }                                                                               \
}

// this macro multiplies a float array by a constant 
//
#define MULT_FLOAT(type, con){                                                      \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        im->data.f[i] *= c;                                                         \
    }                                                                               \
}

// this macro multiplies a double array by a constant 
//
#define MULT_DOUBLE(type, con){                                                     \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        im->data.d[i] *= c;                                                         \
    }                                                                               \
}

// this macro multiplies a int32 array by a constant 
//
#define MULT_INT32(type, con){                                                      \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        im->data.i32[i] *= c;                                                       \
    }                                                                               \
}

// this macro multiplies a vector array by a constant 
//
#define MULT_VECT(type, con){                                                       \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        im->data.vect[i].x *= c;                                                    \
        im->data.vect[i].y *= c;                                                    \
        im->data.vect[i].z *= c;                                                    \
    }                                                                               \
}

// this macro multiplies a polar array by a constant 
//
#define MULT_POL_CONST_O(type, con){                                                \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        out->data.pol[i].pabs = im->data.pol[i].pabs * c;                           \
    }                                                                               \
}

// this macro multiplies a caratestian array by a constant 
//
#define MULT_CMPL_CONST_O(type, con){                                               \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx * im->ny; i++) {                                         \
        out->data.cmpl_f[i].r = im->data.cmpl_f[i].r * c;                           \
        out->data.cmpl_f[i].r = im->data.cmpl_f[i].i * c;                           \
    }                                                                               \
}

// this macro multiplies a float array by a constant 
//
#define MULT_FLOAT_O(type, con){                                                    \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        out->data.f[i] = im->data.f[i] * c;                                         \
    }                                                                               \
}

// this macro multiplies a double array by a constant 
//
#define MULT_DOUBLE_O(type, con){                                                   \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        out->data.d[i] = im->data.d[i] * c;                                         \
    }                                                                               \
}

// this macro multiplies a int32 array by a constant 
//
#define MULT_INT32_O(type, con){                                                    \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        out->data.i32[i] = im->data.i32[i] * c;                                     \
    }                                                                               \
}

// this macro multiplies a vector array by a constant 
//
#define MULT_VECT_O(type, con){                                                     \
    type c = *(type *)con;                                                          \
    for (i = 0; i < im->nx*im->ny; i++) {                                           \
        out->data.vect[i].x =  im->data.vect[i].x * c;                              \
        out->data.vect[i].y =  im->data.vect[i].y * c;                              \
        out->data.vect[i].z =  im->data.vect[i].z * c;                              \
    }                                                                               \
}

// this function performs a macro on a number of cases 
//
#define TYPE_CASES(macro){                                                          \
    switch(c_type)                                                                  \
    {                                                                               \
        case ITYPE_FLOAT:                                                           \
            macro(float, con);                                                      \
            break;                                                                  \
        case ITYPE_DOUBLE:                                                          \
            macro(double, con);                                                     \
            break;                                                                  \
        case ITYPE_INT32:                                                           \
            macro(int32_t, con);                                                    \
            break;                                                                  \
        default:                                                                    \
            fprintf(stderr, "im_mult_scalar: The constant must be a scalar!\n");    \
            status->status = INVALID_TYPE;                                          \
            break;                                                                  \
    }                                                                               \
}

#define TRIG(fn){                                                                   \
    for (i = 0; i <in->nx*in->ny; i++)                                              \
    {                                                                               \
        out->data.d[i] = fn (in->data.d[i]);                                        \
    }                                                                               \
}


// Adds two images and outputs to a third.
//
SPStatus* im_add (SPImage *a, SPImage *b, SPImage *out, SPStatus *status)               // adds the data part of two image structures 
{
    int64_t n;
    
    CHECKS();
    
    switch (a->image_type)
    {
        case ITYPE_POLAR:
            for (n=0; n<((a->nx)*(a->ny)); n++)
            {
                CMPLX_PADD(a->data.pol[n], b->data.pol[n], out->data.pol[n]);
            }
            break;
            
        case ITYPE_CMPL_FLOAT:
            for (n=0; n<((a->nx)*(a->ny)); n++)
            {
                CMPLX_ADD(a->data.cmpl_f[n], b->data.cmpl_f[n], out->data.cmpl_f[n]);
            }
            break;
            
        case ITYPE_VECTOR:
            for(n = 0; n<((a->nx)*(a->ny)); n++)
            {
                VECT_ADD(a->data.vect[n], b->data.vect[n], out->data.vect[n]);
            }
            break;
            
            FORALLTYPES(+);
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_add\n", a->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    return(status);
}

// Subtracts two images and outputs to a third.
//
SPStatus* im_sub (SPImage *a, SPImage *b, SPImage *out, SPStatus *status) 
{
    int64_t n;
    CHECKS();
    
    switch (a->image_type)
    {
        case ITYPE_POLAR:
            for (n=0; n<(a->nx*a->ny); n++)
            {
                CMPLX_PSUB(a->data.pol[n], b->data.pol[n], out->data.pol[n]);
            }
            break;
            
        case ITYPE_CMPL_FLOAT:
            for (n=0; n<(a->nx*a->ny); n++)
            {
                CMPLX_SUB(a->data.cmpl_f[n], b->data.cmpl_f[n], out->data.cmpl_f[n]);
            }
            break;
            
        case ITYPE_VECTOR:
            for(n = 0; n<((a->nx)*(a->ny)); n++)
            {
                VECT_SUB(a->data.vect[n], b->data.vect[n], out->data.vect[n]);
            }
            break;
            
            FORALLTYPES(-);
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_sub\n", a->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    return(status);
}

// Multiplies (pointwise) two images and outputs to a third.
//
SPStatus* im_mult (SPImage *a, SPImage *b, SPImage *out, SPStatus *status)
{
    int64_t n;
    CHECKS();
    
    switch (a->image_type)
    {
        case ITYPE_POLAR:
            for (n=0; n<(a->nx*a->ny); n++)
            {
                CMPLX_PMULT(a->data.pol[n], b->data.pol[n], out->data.pol[n]);
            }
            break;
            
        case ITYPE_CMPL_FLOAT:
        {
            SPCmplx tmp;
            for (n=0; n<(a->nx*a->ny); n++) {
                CMPLX_MULT(a->data.cmpl_f[n], b->data.cmpl_f[n], tmp);
                out->data.cmpl_f[n] = tmp;
            }
        }
            break;
            
            FORALLTYPES(*);
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_mult\n", a->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    return(status);
}

// Divides two images and outputs to a third.
//
SPStatus* im_div (SPImage *a, SPImage *b, SPImage *out, SPStatus *status)
{
    int64_t n;
    CHECKS();
    
    switch(a->image_type)
    {
        case ITYPE_POLAR:
            for (n=0; n<(a->nx*a->ny); n++)
            {
                CMPLX_PDIV(a->data.pol[n], b->data.pol[n], out->data.pol[n]);
            }
            break;
            
        case ITYPE_CMPL_FLOAT:
            for (n=0; n<(a->nx*a->ny); n++)
            {
                CMPLX_DIV(a->data.cmpl_f[n], b->data.cmpl_f[n], out->data.cmpl_f[n]);
            }
            break;
            
            FORALLTYPES(/);
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_div\n", a->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    
    return(status);
}

// Conjgates the input image and outputs to another image.
//
SPStatus* im_conjg (SPImage *a, SPImage *out, SPStatus *status)
{
    int64_t n;
    CHECK_STATUS(status);
    CHECK_PTR(out,status);
    CHECK_IMINTEG(a,status);
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0))
    {
        im_create (out, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_IMSIZEOUT(a, out, status);
    CHECK_TYPES_EQUAL(a, out, status);
    
    switch(a->image_type)
    {
        case ITYPE_POLAR:
            for (n=0; n<(a->nx*a->ny); n++)
            {
                CMPLX_PCONJG(a->data.pol[n], out->data.pol[n]);
            }
            break;
            
        case ITYPE_CMPL_FLOAT:
            for (n=0; n<(a->nx*a->ny); n++)
            {
                CMPLX_CONJG(a->data.cmpl_f[n], out->data.cmpl_f[n]);
            }
            break;
            
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_conjg\n", a->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    return(status);
}

// Multplies the input image by a supplied scalar, output goes to supplied image.
//
SPStatus* im_mult_scalar(SPImage *im, SPImageType c_type, void *con, SPStatus *status)
{
    int64_t i;
    
    CHECK_STATUS(status);
    CHECK_IMINTEG(im,status);
    
    switch (im->image_type)
    {
        case ITYPE_POLAR:
            TYPE_CASES(MULT_POL_CONST);
            break;
        case ITYPE_CMPL_FLOAT:
            TYPE_CASES(MULT_CMPL_CONST);
            break;
        case ITYPE_FLOAT:
            TYPE_CASES(MULT_FLOAT);
            break;
        case ITYPE_DOUBLE:
            TYPE_CASES(MULT_DOUBLE);
            break;
        case ITYPE_INT32:
            TYPE_CASES(MULT_INT32);
            break;
        case ITYPE_VECTOR:
            TYPE_CASES(MULT_VECT);
            break;
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_mult_scalar\n", im->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    return(status);
}

// Multplies the input image by a supplied scalar, output goes to output image.
//
SPStatus* im_mult_scalar_o(SPImage *im, SPImageType c_type, void *con, SPImage *out, SPStatus *status)
{
    int64_t i;
    
    CHECK_STATUS(status);
    CHECK_IMINTEG(im,status);
    
    switch (im->image_type)
    {
        case ITYPE_POLAR:
            TYPE_CASES(MULT_POL_CONST_O);
            break;
        case ITYPE_CMPL_FLOAT:
            TYPE_CASES(MULT_CMPL_CONST_O);
            break;
        case ITYPE_FLOAT:
            TYPE_CASES(MULT_FLOAT_O);
            break;
        case ITYPE_DOUBLE:
            TYPE_CASES(MULT_DOUBLE_O);
            break;
        case ITYPE_INT32:
            TYPE_CASES(MULT_INT32_O);
            break;
        case ITYPE_VECTOR:
            TYPE_CASES(MULT_VECT_O);
            break;
        default:
            fprintf(stderr, "Src image type (%d) not supported in im_mult_scalar\n", im->image_type);
            status->status = INVALID_TYPE;
            break;
    }
    return(status);
}

// subtracts each row of a col of doubles to each row of the data part of an image structure
//
SPStatus* im_colsub (SPImage *a, double *col, SPImage *out, SPStatus *status)
{
    int64_t x;
    int64_t y;
    SPCmplx tmp;
    CHECK_STATUS(status);
    CHECK_PTR(out, status);
    CHECK_IMINTEG(a, status);
    
    if (a->image_type != ITYPE_POLAR) {
        status->status = INVALID_TYPE;
        fprintf(stderr, "Only polar type images ca be used in im_colsub\n");
        return (status);
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny == 0)) {
        im_create (out, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_IMSIZEOUT(a, out, status);
    
    for (y=0; y<a->ny; y++) {
        for (x=0; x<a->nx; x++) {
            CMPLX_POL2CART (a->data.pol[x+y*a->nx], tmp);
            tmp.r = tmp.r - col[y];
            CMPLX_CART2POL (tmp, out->data.pol[x+y*a->nx]);
        }
    }
    return(status);
}

// adds each row of a column of doubles to each row of the data part of an image structure
//
SPStatus* im_coladd (SPImage *a, double *col, SPImage *out, SPStatus *status)
{
    int64_t x;
    int64_t y;
    SPCmplx tmp;
    CHECK_STATUS(status);
    CHECK_PTR(out, status);
    CHECK_IMINTEG(a, status);
    
    if (a->image_type != ITYPE_POLAR) {
        status->status = INVALID_TYPE;
        fprintf(stderr, "Only polar type images ca be used in im_coladd\n");
        return (status);
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0))  {
        im_create (out, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_IMSIZEOUT(a, out, status);
    
    for (y=0; y<a->ny; y++) {
        for (x=0; x<a->nx; x++) {
            CMPLX_POL2CART (a->data.pol[x+y*a->nx], tmp);
            tmp.r = tmp.r - col[y];
            CMPLX_CART2POL (tmp, out->data.pol[x+y*a->nx]);
        }
    }
    
    return(status);
}

// correlate the data part of two image structures
//
SPStatus* im_corr (SPImage *a, SPImage *b, SPImage *out, SPStatus *status)
{
    int64_t x;
    int64_t y;
    CHECK_STATUS(status);
    CHECK_SPC(a, b, status);
    CHECK_IMINTEG(a, status);
    CHECK_IMINTEG(b, status);
    CHECK_ROWS(a, b, status);
    
    if (a->image_type != ITYPE_POLAR) {
        status->status = INVALID_TYPE;
        fprintf(stderr, "Only polar type images ca be used in im_coladd\n");
        return (status);
    }
    
    im_conjg (a, a, status);
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny == 0)) {
        im_create (out, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_IMSIZEOUT(b, out, status);
    
    for (y=0; y<b->ny; y++) {
        for (x=0; x<b->nx; x++) {
            CMPLX_PMULT(a->data.pol[x], b->data.pol[x+b->nx*y], out->data.pol[x+b->nx*y]);
        }
    }
    
    return (status);
}

// add up each row in the image, output a column vector
//
SPStatus* im_sumrows(SPImage *a, SPImage *out, SPStatus *status)                              
{
    CHECK_STATUS(status);
    CHECK_IMINTEG(a, status);
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0)) {
        im_create (out, a->image_type, a->ny, 1, 1, 1, status);
    }
    
    CHECK_ROWS(a, out, status);
    CHECK_TYPES_EQUAL(a, out, status);
    
    FORROWS(a, out, OPER_ADD);
    return(status);
}

// Adds the constant onto the supplied image of some types.  The "scalar" is generic and may be a single vector
// or complex value or doubles
//
SPStatus* im_add_scalar(SPImage * im, SPImageType c_type, void * con, SPStatus * status)
{
    CHECK_STATUS(status);
    CHECK_IMINTEG(im, status);
    int64_t i;
    
    if ((im->image_type == ITYPE_POLAR || im->image_type == ITYPE_CMPL_FLOAT) ^
        (c_type == ITYPE_POLAR || c_type == ITYPE_CMPL_FLOAT)) {
        fprintf(stderr, "im_add_scalar: Either the image type is complex or the constant is complex but not both!\n");
        status->status = INVALID_TYPE;
        return (status);
    }
    
    switch (im->image_type)
    {
        case ITYPE_CMPL_FLOAT:
            switch (c_type)
        {
            case ITYPE_CMPL_FLOAT:
            {
                SPCmplx v = *(SPCmplx *) con;
                for(i = 0; i < im->nx * im->ny; i++)  {
                    im->data.cmpl_f[i].r += v.r;
                    im->data.cmpl_f[i].i += v.i;
                }
            }
                break;
            default:
                fprintf(stderr, "im_add_scalar: Constant type (%s) not supported!\n", itype2string(c_type));
                status->status = INVALID_TYPE;
                break;
        }
            break;
        case ITYPE_VECTOR:
            switch (c_type)
        {
            case ITYPE_VECTOR:
            {
                SPVector v = *(SPVector *) con;
                for(i = 0; i < im->nx * im->ny; i++)  {
                    im->data.vect[i].x += v.x;
                    im->data.vect[i].y += v.y;
                    im->data.vect[i].z += v.z;
                }
            }
                break;
            default:
                fprintf(stderr, "im_add_scalar: Constant type (%s) not supported!\n", itype2string(c_type));
                status->status = INVALID_TYPE;
                break;
        }
            break;
        case ITYPE_DOUBLE:
            switch (c_type)
        {
            case ITYPE_DOUBLE:
            {
                double v = *(double *) con;
                for (i = 0; i < im->nx * im->ny; i++) {
                    im->data.d[i] += v;
                }
            }
                break;
            default:
                fprintf(stderr, "im_add_scalar: Constant type (%s) not supported!\n", itype2string(c_type));
                status->status = INVALID_TYPE;
                break;
        }
            break;
        case ITYPE_FLOAT:
            switch (c_type)
        {
            case ITYPE_FLOAT:
            {
                float v = *(float *) con;
                for (i = 0; i < im->nx * im->ny; i++) {
                    im->data.f[i] += v;
                }
            }
                break;
            default:
                fprintf(stderr, "im_add_scalar: Constant type (%s) not supported!\n", itype2string(c_type));
                status->status = INVALID_TYPE;
                break;
        }
            break;
        default:
            fprintf(stderr, "im_add_scalar: Image type (%s) not supported!\n", itype2string(im->image_type));
            status->status = INVALID_TYPE;
            break;
    }
    return(status);
}

// Adds the constant onto the supplied image of some types.  The "scalar" is generic and may be a single vector
// or complex value or doubles.  This function outputs an image.
//
SPStatus* im_add_scalar_o(SPImage * im, SPImageType c_type, void * con, SPImage *out, SPStatus * status)
{
    int64_t i;
    
    CHECK_STATUS(status);
    CHECK_IMINTEG(im, status);
    
    if ((im->image_type == ITYPE_POLAR || im->image_type == ITYPE_CMPL_FLOAT) ^
        (c_type == ITYPE_POLAR || c_type == ITYPE_CMPL_FLOAT)) {
        fprintf(stderr, "im_add_scalar: Either the image type is complex or the constant is complex but not both!\n");
        status->status = INVALID_TYPE;
        return (status);
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny == 0)) {
        im_create (out, im->image_type, im->ny, 1, 1, 1, status);
    }
    
    CHECK_TYPES_EQUAL(im, out, status);
    
    switch (im->image_type)
    {
        case ITYPE_CMPL_FLOAT:
            switch (c_type)
        {
            case ITYPE_CMPL_FLOAT:
            {
                SPCmplx v = *(SPCmplx *) con;
                for(i = 0; i < im->nx * im->ny; i++)
                {
                    out->data.cmpl_f[i].r = im->data.cmpl_f[i].r + v.r;
                    out->data.cmpl_f[i].i = im->data.cmpl_f[i].r + v.i;
                }
            }
                break;
            default:
                fprintf(stderr, "im_add_scalar: Constant type (%s) not supported!\n", itype2string(c_type));
                status->status = INVALID_TYPE;
                break;
        }
            break;
        case ITYPE_VECTOR:
            switch (c_type)
        {
            case ITYPE_VECTOR:
            {
                SPVector v = *(SPVector *) con;
                for(i = 0; i < im->nx * im->ny; i++)
                {
                    out->data.vect[i].x = im->data.vect[i].x + v.x;
                    out->data.vect[i].y = im->data.vect[i].y + v.y;
                    out->data.vect[i].z = im->data.vect[i].z + v.z;
                }
            }
                break;
            default:
                fprintf(stderr, "im_add_scalar: Constant type (%s) not supported!\n", itype2string(c_type));
                status->status = INVALID_TYPE;
                break;
        }
            break;
        case ITYPE_DOUBLE:
            switch (c_type)
        {
            case ITYPE_DOUBLE:
            {
                double v = *(double *) con;
                for (i = 0; i < im->nx; i++)
                {
                    out->data.d[i] = im->data.d[i] + v;
                }
            }
                break;
            default:
                fprintf(stderr, "im_add_scalar: Constant type (%s) not supported!\n", itype2string(c_type));
                status->status = INVALID_TYPE;
                break;
        }
            break;
        default:
            fprintf(stderr, "im_add_scalar: Image type (%s) not supported!\n", itype2string(im->image_type));
            status->status = INVALID_TYPE;
            break;
    }
    return(status);
}

int64_t
next_power_2(int64_t num)
{
    int64_t val = 1;
    
    if (num >= 1)
    {
        while ( val < num)
        {
            val *= 2;
        }
    }
    
    return val;
    
}


// Multiplies (pointwise) a 2D image by a 1D image.
//
SPStatus* im_mult_line (SPImage *a, SPImage *b, SPDirection dir, SPImage *out, SPStatus *status)
{
    int64_t x, y;
    
    CHECK_STATUS(status);
    CHECK_PTR(out,status);
    CHECK_IMINTEG(a,status);
    CHECK_IMINTEG(b,status);
    CHECK_TYPES_EQUAL(a, b, status);
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0))
    {
        im_create (out, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    CHECK_TYPES_EQUAL(a,out,status);
    CHECK_STATUS(status);
    
    if(status->debug >=10)
    {
        printf("%s:%d\n",__FILE__,__LINE__);
    }
    
    switch (a->image_type)
    {
        case ITYPE_POLAR:
        {
            switch (dir)
            {
                case X_DIR:
                {
                    if (a->nx != b->nx)
                    {
                        status->status = INPUT_SIZES_MISMATCHED;
                        fprintf(stderr, "%s:%d Multiplying 2D array %"PRId64" long by 1D array %"PRId64" long\n", __FILE__, __LINE__, a->nx, b->nx);
                        CHECK_STATUS(status);
                    }
                    for (y = 0; y < a->ny; y++)
                    {
                        for (x = 0; x < a->nx; x++)
                        {
                            CMPLX_PMULT(a->data.pol[x + a->nx*y], b->data.pol[x], out->data.pol[x + a->nx*y]);
                        }
                    }
                }
                    break;
                    
                case Y_DIR:
                {
                    if (a->ny != b->nx)
                    {
                        status->status = INPUT_SIZES_MISMATCHED;
                        fprintf(stderr, "%s:%d Multiplying 2D array %"PRId64" long by 1D array %"PRId64" long\n", __FILE__, __LINE__, a->ny, b->nx);
                        CHECK_STATUS(status);
                    }
                    for (x = 0; x < a->nx; x++)
                    {
                        for (y = 0; y < a->ny; y++)
                        {
                            CMPLX_PMULT(a->data.pol[y + a->ny*x], b->data.pol[x], out->data.pol[y + a->ny*x]);
                        }
                    }
                }
                    break;
                    
                default:
                    fprintf(stderr, "Direction (%s) not supported in im_mult_line\n", dir2string(dir));
                    status->status = INVALID_TYPE;
                    break;
            }
        }
            break;
            
        case ITYPE_CMPL_FLOAT:
        {
            if(status->debug >=10)
            {
                printf("%s:%d\n",__FILE__,__LINE__);
            }
            switch (dir)
            {
                case X_DIR:
                {
                    if(status->debug >=10)
                    {
                        printf("%s:%d a->nx %"PRId64" a->ny %"PRId64" b->nx %"PRId64" b->ny %"PRId64"\n",__FILE__,__LINE__, a->nx, a->ny, b->nx, b->ny);
                    }
                    if (a->nx != b->nx)
                    {
                        status->status = INPUT_SIZES_MISMATCHED;
                        fprintf(stderr, "%s:%d Multiplying 2D array %"PRId64" long by 1D array %"PRId64" long\n", __FILE__, __LINE__, a->nx, b->nx);
                        CHECK_STATUS(status);
                    }
                    for (y = 0; y < a->ny; y++)
                    {
                        for (x = 0; x < a->nx; x++)
                        {
                            CMPLX_MULT(a->data.cmpl_f[x + a->nx*y], b->data.cmpl_f[x], out->data.cmpl_f[x + a->nx*y]);
                            // 		      printf("x %ld y %ld\n",x,y); 
                        }
                    }
                    if(status->debug >=10)
                    {
                        printf("%s:%d\n",__FILE__,__LINE__);
                    }
                }
                    break;
                    
                case Y_DIR:
                {
                    if (a->ny != b->nx)
                    {
                        status->status = INPUT_SIZES_MISMATCHED;
                        fprintf(stderr, "%s:%d Multiplying 2D array %"PRId64" long by 1D array %"PRId64" long\n", __FILE__, __LINE__, a->ny, b->nx);
                        CHECK_STATUS(status);
                    }
                    for (x = 0; x < a->nx; x++)
                    {
                        for (y = 0; y < a->ny; y++)
                        {
                            CMPLX_MULT(a->data.cmpl_f[y + a->ny*x], b->data.cmpl_f[x], out->data.cmpl_f[y + a->ny*x]);
                        }
                    }
                }
                    break;
                    
                default:
                    fprintf(stderr, "Direction (%s) not supported in im_mult_line\n", dir2string(dir));
                    status->status = INVALID_TYPE;
                    break;
            }
        }
            break;
            
        default:
            fprintf(stderr, "Src image type (%s) not supported in im_mult_line\n", itype2string(a->image_type));
            status->status = INVALID_TYPE;
            break;
    }
    
    return(status);
}


// Gives the dot product of each element of two vector images.
// And outputs it as a double.
//
SPStatus* im_dot (SPImage *a, SPImage *b, SPImage *out, SPStatus *status)
{
    int64_t i;
    
    CHECK_STATUS(status);
    CHECK_PTR(out,status);
    CHECK_SPC(a,b,status);
    CHECK_IMINTEG(a,status);
    CHECK_IMINTEG(b,status);
    CHECK_IMSIZEIN(a,b,status);
    CHECK_TYPES_EQUAL(a, b, status);
    if ( (a->image_type != ITYPE_VECTOR) || (b->image_type != ITYPE_VECTOR) ) {
        status->status = INVALID_TYPE;
        fprintf(stderr, "Only vector type images can be used in im_dot\n");
        return (status);
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny == 0)) {
        im_create (out, ITYPE_DOUBLE, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_TYPE(out, ITYPE_DOUBLE, status);
    CHECK_IMSIZEOUT(a, out, status);
    
    CHECK_STATUS(status);
    
    for (i = 0; i < (a->nx*a->ny); i++) {
        out->data.d[i] = VECT_DOT(a->data.vect[i], b->data.vect[i]);
    }
    
    return(status);
}


// Gives the fmod (# man fmod) of each element of an image.
// And outputs it as a double.
//
SPStatus* im_fmod (SPImage *a, double div, SPImage *out, SPStatus *status)          
{
    int64_t i;
    
    CHECK_STATUS(status);                                                            
    CHECK_PTR(out,status);            
    CHECK_IMINTEG(a,status);
    
    if ( a->image_type != ITYPE_DOUBLE && a->image_type != ITYPE_FLOAT) {
        status->status = INVALID_TYPE;
        fprintf(stderr, "Only double or float type images can be used in im_fmod\n");
        return (status);
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0)) {
        im_create (out, a->image_type, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_IMSIZEOUT(a,out,status);
    
    CHECK_STATUS(status);
    
    if (a->image_type == ITYPE_DOUBLE) {
        for (i = 0; i < (a->nx*a->ny); i++) {
            out->data.d[i] = fmod(a->data.d[i], div);
        }
    } else {
        for (i = 0; i < (a->nx*a->ny); i++) {
            out->data.f[i] = fmodf(a->data.f[i], div);
        }
    }
    
    return(status);
}


// Gives the sin (# man sin) of each element of an image.
// And outputs it as a double.
//
SPStatus* im_sin (SPImage *in, SPImage *out, SPStatus *status)          
{
    int64_t i;
    
    CHECK_STATUS(status);                                                            
    CHECK_PTR(out,status);            
    CHECK_IMINTEG(in,status);
    
    CHECK_IMSIZEOUT(in,out,status);  
    if ( in->image_type != ITYPE_DOUBLE) {
        status->status = INVALID_TYPE;
        fprintf(stderr, "Only double type images can be used in im_sin\n");
        return (status);
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny ==0)) {
        im_create (out, ITYPE_DOUBLE, in->nx, in->ny, in->xspc, in->yspc, status);
    }
    
    CHECK_STATUS(status);
    
    TRIG(sin);
    
    return(status);
}


// Gives the cos (# man cos) of each element of an image.
// And outputs it as a double.
//
SPStatus* im_cos (SPImage *in, SPImage *out, SPStatus *status)          
{
    int64_t i;
    
    CHECK_STATUS(status);                                                            
    CHECK_PTR(out,status);            
    CHECK_IMINTEG(in,status);
    CHECK_IMSIZEOUT(in,out,status);  
    
    if ( in->image_type != ITYPE_DOUBLE) {
        status->status = INVALID_TYPE;
        fprintf(stderr, "Only double type images can be used in im_sin\n");
        return (status);
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny == 0)) {
        im_create (out, ITYPE_DOUBLE, in->nx, in->ny, in->xspc, in->yspc, status);
    }
    
    CHECK_STATUS(status);
    
    TRIG(cos);
    
    return(status);
}


// Gives the rint (# man rint) of each element of an image.
// And outputs it as an int64_t.
//
SPStatus* im_rint (SPImage *a, SPImage *out, SPStatus *status)          
{
    int64_t i;
    
    CHECK_STATUS(status);                                                            
    CHECK_PTR(out,status);            
    CHECK_IMINTEG(a,status);
    if ( a->image_type != ITYPE_DOUBLE) {
        status->status = INVALID_TYPE;
        fprintf(stderr, "Only double type images can be used in im_rint\n");
        return (status);
    }
    
    if ((out->data.v == NULL) && (out->nx == 0) && (out->ny == 0)) {
        im_create (out, ITYPE_INT64, a->nx, a->ny, a->xspc, a->yspc, status);
    }
    
    CHECK_IMSIZEOUT(a,out,status);  
    
    CHECK_STATUS(status);
    
    for (i = 0; i < (a->nx*a->ny); i++) {
        out->data.i64[i] = (int64_t) (rint(a->data.d[i]));
    }
    
    return(status);
}


