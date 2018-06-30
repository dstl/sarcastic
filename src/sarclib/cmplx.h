/** @file********************************************************************
 *
 *       Module:    cmplx.h
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 20/10/2004.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      a function for making a complex cartesian float image from an image which
 *      contains the real part and an image which contains the imaginary part.
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

#ifndef sarclib_CMPLX_H__
#define sarclib_CMPLX_H__

/// cartesian structure
///
typedef struct {
    float r;
    float i;
} SPCmplx;                                             

/// polar structure
///
typedef struct {
    float pabs;
    float parg;
} SPCmplxPol;                                          

/// Complex int64
///
typedef struct {
    int64_t r;
    int64_t i;
} SPCmplxInt64;                                          

/// Complex int32
///
typedef struct {
    int32_t r;
    int32_t i;
} SPCmplxInt32;                                         

/// Complex int short/int16 - be careful when using macros below
/// with this type - you'll need to think about overflows & rounding
///
typedef struct {
    int16_t r;
    int16_t i;
} SPCmplxInt16;

/// Complex bytes/int8 - be careful when using macros below
/// with this type - you'll need to think about overflows & rounding
///
typedef struct {
    int8_t r;
    int8_t i;
} SPCmplxInt8;

/// This is just reallt a place holder till I can think of a way to store
/// complex nibbles (4bits I, 4 bits Q) in a neat way. Will probabky wait until
/// there is a real requirement for them, not just cos CPHD has them as a
/// data type
///
typedef struct {
    int8_t data;
} SPCmplxNibb;

// Macros for cartesian complex numbers
//
/// Make a cartesian complex number
///
#define CMPLX_F_MAKE(a,b,out){out.r = (float)a; out.i = (float)b;}

/// Add two complex numbers together
///
#define CMPLX_ADD(a,b,out) {out.r = a.r + b.r; out .i = a.i + b.i;}

/// Subtract two complex numbers
///
#define CMPLX_SUB(a,b,out) {out.r = a.r - b.r; out .i = a.i - b.i;}

/// Multiply two complex numbers together
///
#define CMPLX_MULT(a,b,out) {out.r = a.r * b.r - a.i * b.i; out.i = a.i * b.r + a.r * b.i;}

/// Perform complex conjugate multiplication (A * Conjg(B) )
///
#define CMPLX_CONJG_MULT(a,b,out) {out.r = a.r * b.r + a.i * b.i; out.i = a.i * b.r - a.r * b.i;} 

/// Perform complex conjugate of a complex number
///
#define CMPLX_CONJG(a,out) {out.r = a.r; out.i = -a.i;}

/// Calculate complex magnitude
///
#define CMPLX_MAG(a) (sqrt(a.r * a.r + a.i * a.i))

/// Calculate the phase of a complex number
///
#define CMPLX_PHASE(a) (atan2(a.i, a.r))

/// Perform complex division
///
#define CMPLX_DIV(a,b,out) {SPCmplxPol tmp1, tmp2, tmp3; CMPLX_CART2POL(a,tmp1); CMPLX_CART2POL(b,tmp2); CMPLX_PDIV(tmp1,tmp2,tmp3); CMPLX_POL2CART(tmp3,out);}

/// Complex scalar multiply. multiply complex number cmp by scalar x
///
#define CMPLX_SCMULT(x,cmp,out) {out.r = x * cmp.r; out.i = x * cmp.i;}

/// convert cartesian complex number to polar complex number
///
#define CMPLX_CART2POL(a,out) {out.pabs = sqrt (a.r * a.r + a.i * a.i); out.parg = atan2 (a.i,a.r);}

// Macros for polar complex numbers
//
/// Convert polar complex number to cartesian complex number
///
#define CMPLX_POL2CART(a,out) {out.r = a.pabs * cos(a.parg); out.i = a.pabs * sin(a.parg);}

/// Perform addition of two polar complex numbers
///
#define CMPLX_PADD(a,b,out) {SPCmplx tempa, tempb, tempc; CMPLX_POL2CART(a,tempa); CMPLX_POL2CART(b,tempb); tempc.r = tempa.r + tempb.r; tempc.i = tempa.i + tempb.i; CMPLX_CART2POL(tempc,out);}

/// Perform subtraction of two polar complex numbers
///
#define CMPLX_PSUB(a,b,out) {SPCmplx tempa, tempb, tempc; CMPLX_POL2CART(a,tempa); CMPLX_POL2CART(b,tempb); tempc.r = tempa.r - tempb.r; tempc.i = tempa.i - tempb.i; CMPLX_CART2POL(tempc,out);}

/// Perform complex multiplication of two polar complex numbers
///
#define CMPLX_PMULT(a,b,out) {out.pabs = a.pabs * b.pabs; out.parg = a.parg + b.parg;}

/// Perform complex conjugate of a polar complex number
///
#define CMPLX_PCONJG(a,out) {out.pabs = a.pabs; out.parg = -a.parg;}

/// perform complex division of two polar complex numbers
///
#define CMPLX_PDIV(a,b,out) {out.pabs = a.pabs/b.pabs; out.parg = a.parg - b.parg;}

/// Perform scaler multiplication of a polar complex number
///
#define CMPLX_PSCMULT(x,cmp,out) {SPCmplx tmpa, tmpb; CMPLX_POL2CART(cmp,tmpa); CMPLX_SCMULT(x,tmpa,tmpb); CMPLX_CART2POL(tmpb,out);}


#endif
