/***************************************************************************
 * 
 *           Module :  cpuPOField.hpp
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 25/01/2015
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Function to calculate the E field at the reciever
 *      using Physical Optics
 *
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

#ifndef cpuPOField_hpp
#define cpuPOField_hpp

#include <stdio.h>
#include <sarclib/sarclib.h>
#include "TriangleMesh.hpp"
#include "threadCore.hpp"

//#define RXPOL "V" // Can be V or H or - (both)
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
                int                 rxPol,              // Receive Polarisation. one of enums set in threadCore.hpp::polarisation
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
