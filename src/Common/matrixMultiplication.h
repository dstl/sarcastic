/***************************************************************************
 * 
 *           Module :  matrixMultiplication.h
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 03/10/2014
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *      Description :  Matrix multiplication algorithms
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

#ifndef sarcastic_matrixMultiplication_h
#define sarcastic_matrixMultiplication_h
#include <sarclib/sarclib.h>

void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By) ;
void mat3by3inv(double *A, double *O) ;
void reduction(double a[][6],int size,int pivot ,int col) ;

void matmulCmplx(SPCmplx *A,SPCmplx *B, SPCmplx *O, int Ax, int Ay,int Bx, int By) ;
void matmulSCCmplx(double *A,SPCmplx *B, SPCmplx *O, int Ax, int Ay,int Bx, int By) ;
void mat3by3invCmplx(SPCmplx *A, SPCmplx *O) ;
void reductionCmplx(SPCmplx a[][6],int size,int pivot ,int col) ;


#endif
