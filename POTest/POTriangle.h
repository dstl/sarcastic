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

typedef struct SPVectorCmplx{
    SPCmplx x ;
    SPCmplx y ;
    SPCmplx z ;
} SPVectorCmplx;


/// Multiply a vector by a constant
///
#define VECTCMPLX_SCMULT(inVect,inScal,outVect)  {  \
    outVect.x.r = (inVect.x.r)*inScal;              \
    outVect.x.i = (inVect.x.i)*inScal;              \
    outVect.y.r = (inVect.y.r)*inScal;              \
    outVect.y.i = (inVect.y.i)*inScal;              \
    outVect.z.r = (inVect.z.r)*inScal;              \
    outVect.z.i = (inVect.z.i)*inScal;              \
}

/// Find the dot product of two vectors
///
#define VECTCMPLX_DOT(aVect,bVect, out) {                                   \
    out.r = aVect.x.r+bVect.x.r+aVect.y.r+bVect.y.r+aVect.z.r+bVect.z.r ;   \
    out.i = aVect.x.i+bVect.x.i+aVect.y.i+bVect.y.i+aVect.z.i+bVect.z.i ;   \
}

/// find a cross product of two vectors
///
#define VECTCMPLX_CROSS(aVect,bVect,outVect) {                                                              \
    outVect.x.r=(aVect.y.r*bVect.z.r)-(aVect.y.i*bVect.z.i)-(aVect.z.r*bVect.y.r)+(aVect.z.i*bVect.y.i) ;   \
    outVect.x.i=(aVect.y.r*bVect.z.i)+(aVect.y.i*bVect.z.r)-(aVect.z.r*bVect.y.i)-(aVect.z.i*bVect.y.r) ;   \
    outVect.y.r=(aVect.z.r*bVect.x.r)-(aVect.z.i*bVect.x.i)-(aVect.x.r*bVect.z.r)+(aVect.x.i*bVect.z.i) ;   \
    outVect.y.i=(aVect.z.r*bVect.x.i)+(aVect.z.i*bVect.x.r)-(aVect.x.r*bVect.z.i)-(aVect.x.i*bVect.z.r) ;   \
    outVect.z.r=(aVect.x.r*bVect.y.r)-(aVect.x.i*bVect.y.i)-(aVect.y.r*bVect.x.r)+(aVect.y.i*bVect.x.i) ;   \
    outVect.z.i=(aVect.x.r*bVect.y.i)+(aVect.x.i*bVect.y.r)-(aVect.y.r*bVect.x.i)-(aVect.y.i*bVect.x.r) ;   \
}

/// Find the modulus (or magnitdue) of a complex vector
/// using De Moivre's formula
///
#define VECTCMPLX_MOD(aVect, outVect) {         \
    SPCmplx ___ax,___ay,___az, ___sum ;         \
    SPCmplxPol ___pol ;                         \
    float ___r , ___a;                          \
    CMPLX_MULT(aVect.x,aVect.x,ax);             \
    CMPLX_MULT(aVect.y,aVect.y,ay);             \
    CMPLX_MULT(aVect.z,aVect.z,az);             \
    CMPLX_ADD(ax, ay, ___sum);                  \
    CMPLX_ADD(az, ___sum, ___sum) ;             \
    CMPLX_CART2POL(___sum, ___pol) ;            \
    ___r = sqrt(___pol.pabs) ;                  \
    ___a = ___pol.parg / 2 ;                    \
    outVect.r = ___r * cos(___a) ;              \
    outVect.i = ___r * sin(___a) ;              \
}

#define VECTCMPLX_NORM(aVect,outVect) {                              \
double vect_unit_tmp = 1.0 / VECT_MOD(aVect);               \
outVect.x = aVect.x*vect_unit_tmp;                          \
outVect.y = aVect.y*vect_unit_tmp;                          \
outVect.z = aVect.z*vect_unit_tmp;                          \
}


void POTriangle(triangle tri, Ray ray, SPVector HitPoint, SPVector obsPnt, double lambda, SPCmplx *EsV, SPCmplx *EsH) ;

SPCmplx surfaceIntegral (double k, triangle tri, SPVector uvw_sg) ;
void EField(double k, double r, triangle tri, Ray ray, SPCmplx Ic, SPVector Es_parrdir, SPVector Es_perpdir, SPCmplx *Es_parr, SPCmplx *Es_perp );

SPCmplx G0_func(double gamma);
SPCmplx G1_func(double gamma);
SPCmplx G2_func(double gamma);
SPCmplx G3_func(double gamma);
int  factorial(int n);
void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By);
void mat3by3inv(double *A, double *O);
void reduction(double a[][6],int size,int pivot ,int col);


#endif
