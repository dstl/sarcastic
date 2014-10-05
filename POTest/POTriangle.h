//
//  POTriangle.h
//  sarcastic
//
//  Created by Darren Muff on 03/10/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#ifndef sarcastic_POTriangle_h
#define sarcastic_POTriangle_h
#include "materialProperties.h"

typedef struct triangle{
    int        id ;
    SPVector   AA ;
    SPVector   BB ;
    SPVector   CC ;
    SPVector   NN ;
    double     area ;
    double     globalToLocalMat[9];
    double     localToGlobalMat[9];
    double     Rs ;
    char       mat[MATBYTES] ;

} triangle;

typedef struct Ray {
    SPVector org ;               // Origin of ray
    SPVector dir ;               // unit vector direction of ray
    double   pow ;               // power of ray at origin
    double   len ;               // length or distance travelled by ray to get to origin
    SPVector pol ;               // unit vector of direction of E field of ray
} Ray ;

SPCmplx surfaceIntegral (double k, triangle tri, SPVector uvw_s) ;
void surfaceCurrent(triangle tri, Ray ray, double * Jx_dashdash, double * Jy_dashdash) ;

SPCmplx G0_func(double gamma);
SPCmplx G1_func(double gamma);
SPCmplx G2_func(double gamma);
SPCmplx G3_func(double gamma);
int  factorial(int n);
void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By);
void mat3by3inv(double *A, double *O);
void reduction(double a[][6],int size,int pivot ,int col);


#endif
