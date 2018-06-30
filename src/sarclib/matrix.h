/** @file********************************************************************
 *
 *       Module:    matrix.h
 *      Program:    sarclib
 *   Created by:    Emma Griffiths on 26/01/2005.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      This file contains information for matrix functions
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

#ifndef sarclib_MATRIX_H__
#define sarclib_MATRIX_H__

#include "sarclib.h"

/// defines the structure for a matrix
///
typedef struct {
    size_t m;
    size_t n;
    double *cont;
} SPMatrix;                    

#define MATRIX_CREATE(out,mx,nx) {out.m = mx; out.n = nx; out.cont = calloc(((out.m)*(out.n)), sizeof(double));}
#define MATRIX_DESTROY(a) {free(a.cont); MATRIX_INIT(a);}
#define MATRIX_INIT(a) {a.m = 0; a.n = 0; a.cont = NULL;}
/// this macro is only used once in the change of coordinates
///
#define MATRIX_3VMULT(mat,vect,out,s) {CHECK_MAT_3(mat,s); out.x = mat.cont[0] * vect.x + mat.cont[1] * vect.y + mat.cont[2] * vect.z; out.y = mat.cont[3] * vect.x + mat.cont[4] * vect.y + mat.cont[5] * vect.z; out.z = mat.cont[6] * vect.x + mat.cont[7] * vect.y + mat.cont[8] * vect.z;}

#endif
