//
//  main.c
//  POTest
//
//  Created by Darren on 21/05/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib.h>
#include "sarcastic.h"

// Create a flat plate
// Fire a ray in
// rotate around the object in Az and El
// For each position print out a point where R is teh Es field magnitude
//

int main(int argc, const char * argv[])
{
    
    // Set up a triangle
    //
    SPVector AA,BB,CC ;
    VECT_CREATE(10, 10, 0, AA);
    VECT_CREATE(5, 10, 0, BB);
    VECT_CREATE(10, 5, 0, CC);
    
    // Create a ray
    //
    Ray r ;
    SPVector hp,dir;
    VECT_CREATE(0, 0, 20, r.o);
    VECT_CREATE(7, 7, 0, hp);
    VECT_SUB(hp, r.o, dir);
    VECT_NORM(dir, r.d);
    r.p = 1.0;
    
    // Create observation point
    //
    SPVector op;
    float u,v,w ;
    float theta_s, phi_s ; // Scattering azimuth, elevation
    VECT_CREATE(-10, -10, -10, op);
    theta_s = atan2f(op.y, op.x);
    phi_s   = atan2f(op.z, sqrt(op.x*op.x+op.y*op.y));
    u     = cosf(phi_s) * cosf(theta_s);
    v     = cosf(phi_s) * sinf(theta_s);
    w     = sinf(phi_s);
    
    // Set up the problem parameters
    float Dp,Dq,D0;
    float Cp,Cq,C0;
    float k;
    float lambda = 0.04 ;
    float Z0 ;
    float r = 1000.0;
    float A ; // Area of triangle
    
    
    // First calculate Ic
    //
    k = 2 * SIPC_pi / lambda ;
    Cp = Cq = 0;
    C0 = 1;
    
    Dp = k * ( (AA.x - CC.x) * u + (AA.y - CC.y) * v + (AA.z - CC.z) * w) ;
    Dq = k * ( (BB.x - CC.x) * u + (BB.y - CC.y) * v + (BB.z - CC.z) * w) ;
    D0 = k * ( CC.x * u + CC.y * v + CC.z * w ) ;
    
    SPVector l1,l3,Av,triN;
    float a,modl1,modl3;
    VECT_SUB(BB, AA, l1) ;
    VECT_SUB(AA, CC, l3) ;
    VECT_CROSS(l1, l3, Av) ;
    A = 0.5 * VECT_MAG(Av) ;
    modl1 = VECT_MOD(l1);
    modl3 = VECT_MOD(l3);
    VECT_SCMULT(Av, (1.0/(modl1*modl3)), triN);
    
    float Lt = 0.05 ;  // Length of Taylor series region
    float magDp = fabs(Dp) ;
    float magDq = fabs(Dq) ;
    SPCmplx e_jDp, e_jDq, e_jD0 ;
    e_jD0.r = cosf(D0); e_jD0.i = sinf(D0) ;
    e_jDp.r = cosf(Dp); e_jDp.i = sinf(Dp) ;
    e_jDq.r = cosf(Dq); e_jDq.i = sinf(Dq) ;
    SPCmplx jDp, jD0, jDq;
    jDp.r = 0; jDp.i = Dp ;
    jD0.r = 0; jD0.i = D0 ;
    jDq.r = 0; jDq.i = Dq ;
    
    
    
    return 0;
}


