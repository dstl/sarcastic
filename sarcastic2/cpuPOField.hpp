//
//  cpuPOField.hpp
//  sarcastic
//
//  Created by Darren Muff on 17/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef cpuPOField_hpp
#define cpuPOField_hpp

#include <stdio.h>
#include <SILib2/SILib2.h>
#include "TriangleMesh.hpp"
#include "threadCore.hpp"

#define RXPOL "V" // Can be V or H or - (both)
#define SURFMAXN 4
#define LT 0.04

typedef struct HitPoint {
    SPVector hit;       // Location of hitpoint in x,y,z
    int tri;            // index of triangle that this hit is on
} HitPoint ;

void cpuPOField(TriangleMesh       *mesh,
                Hit                 *hits,              // Array of hit locations to x-ref with triangles for material props
                int                 nRays,              // The number of reflected rays being considered
                Ray                 *rays,              // unit vector rays arriving at hitpoint
                Ray                 *shadowRays,        // Array of reflected rays - used for their origin as its the reflection point to Rx
                SPVector            RxPos,              // Location of Receiver in x,y,z
                double              k,                  // Wavenumber constant k = 2 * PI / Lambda
                double              *ranges,            // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                double              gainRx,             // Receiver gain used for power calculations
                int                 firstBounce,        // if 1 then PO calcs use origin for field calculations
                rangeAndPower       *rnp                // Output array of ray power at, and range to reciever
);

void POKernelCode(int ind, TriangleMesh * mesh, Ray *rays, int nrays, HitPoint *hitpoints, SPVector RxPnt, double k, SPVector Vdir, SPVector Hdir, int firstBounce, SPCmplx *EsVs, SPCmplx *EsHs) ;

// Forward declarations here
//
int factorial(int n) ;
SPCmplx G_func4(double gamma);
SPCmplx G_func3(double gamma);
SPCmplx G_func2(double gamma);
SPCmplx G_func1(double gamma);
SPCmplx G_func0(double gamma);
SPVector translateVector(SPVector v, double *matrix) ;
void matmul2(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By) ;
void surfaceCurrents(SPVector Eig, SPVector Ri_hat, Triangle tri, SPVector Vdir, SPVector Hdir, double *Jpar, double *Jper) ;
SPCmplx ludwigIntegral(double Dp,double Dq,double D0, double A) ;

static SPVector zz_hat = {0.0, 0.0, 1.0};
static double Z0 = 376.99111843077516; // Impedence of free space = 120 * PI

#endif /* cpuPOField_hpp */
