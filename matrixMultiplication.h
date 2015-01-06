//
//  matrixMultiplication.h
//  sarcastic
//
//  Created by Darren Muff on 03/10/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#ifndef sarcastic_matrixMultiplication_h
#define sarcastic_matrixMultiplication_h
#include <SIlib.h>

void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By) ;
void mat3by3inv(double *A, double *O) ;
void reduction(double a[][6],int size,int pivot ,int col) ;

void matmulCmplx(SPCmplx *A,SPCmplx *B, SPCmplx *O, int Ax, int Ay,int Bx, int By) ;
void matmulSCCmplx(double *A,SPCmplx *B, SPCmplx *O, int Ax, int Ay,int Bx, int By) ;
void mat3by3invCmplx(SPCmplx *A, SPCmplx *O) ;
void reductionCmplx(SPCmplx a[][6],int size,int pivot ,int col) ;


#endif
