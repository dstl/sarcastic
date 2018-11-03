/***************************************************************************
 * 
 *           Module :  matrixMultiplication.c
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

#include <stdio.h>
#include "matrixMultiplication.h"

// Matrix multiplication for a N size array
// Note : in order to remove the test for incompatable dimensions (for speed)
// it is essential that Ax = By
// Also the output matrix *O must have been allocated to the correct dimensions
//  double *A - matrix A
//  double *B - matrix B. By must equal Ax
//  double *O - output matrix. Dimensions must be Ay x Bx
//  int Ax  - Number of columns in A
//  int Ay  - Number of rows in A
//  int Bx  - number of columns in B
//  int By  - number of rows in B
//
void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By)
{
    // Ax must be equal to By !!!
    //
    //    if (Ax != By) {
    //        printf("Error : Multiplication by incompatable matrices\n");
    //        exit(1);
    //    }
    
    int i,j,k;
    
    for (i=0;i<Ay;i++){
        for(j=0;j<Bx;j++){
            
            O[i*Bx+j] = 0;
            for(k=0;k<Ax;k++)
                O[i*Bx+j] += A[i*Ax+k] * B[k*Bx+j];
        }
    }
    return ;
}

// 3x3 matrix inversion taken from
// http://www.c4learn.com/c-programs/c-program-to-find-inverse-of-3-x-3.html
// which has a very good graphical description
// double *A - pointer to an array of 9 doubles arranged {{c0r0, c1r0, c2r0},{c0r1,c1r1,c2r1},{c0r2,c1r2,c2r2}} c=column; r=row
// double *O - pointer to output array of 9 doubles arranged {{c0r0, c1r0, c2r0},{c0r1,c1r1,c2r1},{c0r2,c1r2,c2r2}} c=column; r=row
//
void mat3by3inv(double *A, double *O){
    double a[3][6];
    int y,x,i,j;
    for(y=0;y<3;y++){    // Append Unit Matrix
        for(x=0;x<6;x++){
            if(x<3){
                a[y][x] = A[y*3+x] ;
            }else if(x==y+3){
                a[y][x]=1;
            }else{
                a[y][x]=0;
            }
        }
    }
    for(i=0;i<3;i++){
        reduction(a,3,i,i);
    }
    
    for (i=0;i<3; i++){
        for(j=0; j<3; j++){
            O[i*3+j] = a[i][j+3] ;
        }
    }
    return ;
}

//                   a         3        i          i
void reduction(double a[][6],int size,int pivot ,int col)
{
    int i,j;
    double factor;
    
    factor=a[pivot][col];
    
    for(i=0;i<2*size;i++){
        a[pivot][i]/=factor;
    }
    
    for(i=0;i<size;i++){
        if(i!=pivot){
            factor=a[i][col];
            for(j=0;j<2*size;j++)
                a[i][j]=a[i][j]-a[pivot][j]*factor;
        }
    }
    return ;
}

// Complex Matrix multiplication for a N size array
// Note : in order to remove the test for incompatable dimensions (for speed)
// it is essential that Ax = By
// Also the output matrix *O must have been allocated to the correct dimensions
//  double *A - matrix A
//  double *B - matrix B. By must equal Ax
//  double *O - output matrix. Dimensions must be Ay x Bx
//  int Ax  - Number of columns in A
//  int Ay  - Number of rows in A
//  int Bx  - number of columns in B
//  int By  - number of rows in B
//
void matmulCmplx(SPCmplx *A,SPCmplx *B, SPCmplx *O, int Ax, int Ay,int Bx, int By)
{
    // Ax must be equal to By !!!
    //
    //    if (Ax != By) {
    //        printf("Error : Multiplication by incompatable matrices\n");
    //        exit(1);
    //    }
    
    int i,j,k;
    SPCmplx tmp1 ;

    
    for (i=0;i<Ay;i++){
        for(j=0;j<Bx;j++){
            
            O[i*Bx+j].r = O[i*Bx+j].r = 0;
            for(k=0;k<Ax;k++){
                CMPLX_MULT(A[i*Ax+k], B[k*Bx+j], tmp1) ;
                CMPLX_ADD(O[i*Bx+j], tmp1, O[i*Bx+j]) ;
            }
        }
    }
    return ;
}

void matmulSCCmplx(double *A,SPCmplx *B, SPCmplx *O, int Ax, int Ay,int Bx, int By)
{
    // Ax must be equal to By !!!
    //
    //    if (Ax != By) {
    //        printf("Error : Multiplication by incompatable matrices\n");
    //        exit(1);
    //    }
    
    int i,j,k;
    SPCmplx tmp1 ;
    
    
    for (i=0;i<Ay;i++){
        for(j=0;j<Bx;j++){
            
            O[i*Bx+j].r = O[i*Bx+j].r = 0;
            for(k=0;k<Ax;k++){
                CMPLX_SCMULT(A[i*Ax+k], B[k*Bx+j], tmp1) ;
                CMPLX_ADD(O[i*Bx+j], tmp1, O[i*Bx+j]) ;
            }
        }
    }
    return ;
}


// 3x3 complex matrix inversion taken from
// http://www.c4learn.com/c-programs/c-program-to-find-inverse-of-3-x-3.html
// which has a very good graphical description
// double *A - pointer to an array of 9 doubles arranged {{c0r0, c1r0, c2r0},{c0r1,c1r1,c2r1},{c0r2,c1r2,c2r2}} c=column; r=row
// double *O - pointer to output array of 9 doubles arranged {{c0r0, c1r0, c2r0},{c0r1,c1r1,c2r1},{c0r2,c1r2,c2r2}} c=column; r=row
//
void mat3by3invCmplx(SPCmplx *A, SPCmplx *O){
    SPCmplx a[3][6];
    int y,x,i,j;
    for(y=0;y<3;y++){    // Append Unit Matrix
        for(x=0;x<6;x++){
            if(x<3){
                a[y][x] = A[y*3+x] ;
            }else if(x==y+3){
                a[y][x].r = a[y][x].i = 1 ;
            }else{
                a[y][x].r = a[y][x].i = 0;
            }
        }
    }
    for(i=0;i<3;i++){
        reductionCmplx(a,3,i,i);
    }
    
    for (i=0;i<3; i++){
        for(j=0; j<3; j++){
            O[i*3+j] = a[i][j+3] ;
        }
    }
    return ;
}

//                       a         3        i          i
void reductionCmplx(SPCmplx a[][6],int size,int pivot ,int col)
{
    int i,j;
    SPCmplx factor, tmp1;
    
    factor=a[pivot][col];
    
    for(i=0;i<2*size;i++){
        CMPLX_DIV(a[pivot][i], factor, a[pivot][i]) ;
    }
    
    for(i=0;i<size;i++){
        if(i!=pivot){
            factor=a[i][col];
            for(j=0;j<2*size;j++){
                CMPLX_MULT(a[pivot][j], factor, tmp1) ;
                CMPLX_SUB(a[i][j], tmp1, a[i][j]) ;
            }
        }
    }
    return ;
}

