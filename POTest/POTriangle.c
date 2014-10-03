//
//  POTriangle.c
//  sarcastic
//
//  Created by Darren on 24/08/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib.h>
#include <math.h>
#include "POTriangle.h"
#include "matrixMultiplication.h"


int POTriangle(triangle tri, Ray ray, SPVector HitPoint, SPVector ObservationPoint, double lambda, SPCmplx *Es){


    SPVector l1,l3,Av,triN,tmp1,tmp2;
    double a,k,modl1,modl3, Area;
    double baryLambda_1,baryLambda_2,baryLambda_3;


    k = 2 * SIPC_pi / lambda ;
    
    // Sort out the direction cosines for this ray, hitpoint and observation point
    //
    
    SPVector uvw_i ; // Direction cosines for incident ray (ray origin to HitPoint)
    SPVector uvw_s ; // Direction cosines for scattering ray (HitPoint  to observation point)
    double u_i,v_i,w_i ;
    double u_s,v_s,w_s ;
    double theta_i, phi_i ; // incident elevation, azimuth
    double theta_s, phi_s ; // Scattering elevation, azimuth
    
    VECT_SUB(HitPoint, ray.org, tmp1);
    
    phi_i     = atan2(tmp1.y,tmp1.x);
    theta_i   = atan2(sqrt(tmp1.x*tmp1.x+tmp1.y*tmp1.y),tmp1.z);
    VECT_CREATE(sin(theta_i)*cos(phi_i), sin(theta_i)*sin(phi_i), cos(theta_i), uvw_i) ;
    
    phi_s   = atan2(ObservationPoint.y, ObservationPoint.x);
    theta_s = atan2(sqrt(ObservationPoint.x*ObservationPoint.x+ObservationPoint.y*ObservationPoint.y),ObservationPoint.z);
    VECT_CREATE(sin(theta_s) * cos(phi_s), sin(theta_s) * sin(phi_s), cos(theta_s), uvw_s) ;


    // Now need toperform two seperate calculations:
    // 1) The surface current over the triangle in terms of Jx & Jy
    // 2) The surface integral over triangle 'c' called Ic
    //
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    // 1) Find the surface current over the triangle in terms of Jx & Jy
    //
    
    double Jxdd, Jydd ;             // Jx'' & Jy''
    surfaceCurrent(&Jxdd, &Jydd) ;
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    // 2) Find the value of Ic
    //
    
    SPCmplx Ic = surfaceIntegral(k, tri, uvw_i, uvw_s ) ;
    
    return 0;
}


void surfaceCurrent(triangle tri, Ray ray, double * Jx_dashdash, double * Jy_dashdash){
    
    
    SPVector Eig = ray.pol ;
    SPVector dir = ray.dir ;
    SPVector Eil ;
    
    
    // Convert the incident ray into the coordinate system of the triangle
    // use _l for coordinates in local frame of facet and _g or nothing for global
    // world coordinates
    //
    
    double theta_il, theta_sl ; // angle to Z axis
    double phi_il, phi_sl ;
    double dir_ig[3], dir_il[3], dir_sg[3], dir_sl[3] ; // Direction cosines for incident and scattered Ray
    
    // Assuming that the direction of the ray has been normalised
    //
    dir_ig[0] = -dir.x;
    dir_ig[1] = -dir.y;
    dir_ig[2] = -dir.z;
    
    matmul(tri.globalToLocalMat, dir_ig, dir_il, 3, 3, 1, 3);
    
    theta_il = acos(dir_il[2]) ;
    phi_il   = atan2(dir_il[1], dir_il[0]);
    
    // Also convert the direction of the E field into the coordinate system
    // of the triangle
    //
    
    double E_ig[3], E_il[3] ;
    E_ig[0] = Eig.x ;
    E_ig[1] = Eig.y ;
    E_ig[2] = Eig.z ;
    
    matmul(tr.globalToLocalMat, E_ig, E_il, 3, 3, 1, 3) ;
    VECT_CREATE(E_il[0], E_il[1], E_il[2], Eil) ;
    SPVector theta_l_hat, phi_l_hat, z_l_hat ;
    VECT_CREATE(0, 0, 1, z_l_hat) ;
    VECT_CROSS(ray.dir, z_l_hat, phi_l_hat) ;
    VECT_CROSS(phi_l_hat, ray.dir, theta_l_hat) ;
    VECT_NORM(phi_l_hat, phi_l_hat);
    VECT_NORM(theta_l_hat, theta_l_hat) ;
    
    double Eiphi_l   = VECT_DOT(Eil, phi_l_hat) ;
    double Eitheta_l = VECT_DOT(Eil, theta_l_hat) ;
    
    double Rs, Z0, GamParr, GamPerp ;
    Rs = tri.Rs ;
    Z0 = 120 * SIPC_pi ; // Impedence of free space
    
    GamParr = -1.0 * Z0 * cos(theta_il) / (2*Rs + Z0*cos(theta_il));
    GamPerp = -1.0 * Z0 / ( 2.0*Rs*cos(theta_il) + Z0);
    
    
    
    // Need E_itheta_dashdash
    //
    

}


SPCmplx surfaceIntegral (double k, triangle tri, SPVector uvw_s){
    
    double Dp,Dq,D0;
    double Cp,Cq,C0;
    double u = uvw_s.x ;
    double v = uvw_s.y ;
    double w = uvw_s.z ;
    double x1 = tri.AA.x ;
    double y1 = tri.AA.y ;
    double z1 = tri.AA.z ;
    double x2 = tri.BB.x ;
    double y2 = tri.BB.y ;
    double z2 = tri.BB.z ;
    double x3 = tri.CC.x ;
    double y3 = tri.CC.y ;
    double z3 = tri.CC.z ;
    double A  = tri.area ;
    
    Cp = Cq = 0;
    C0 = 1;
    Dp = k * ( ((x1 - x3) * u) + ((y1 - y3) * v) + ((z1 - z3) * w) ) ;
    Dq = k * ( ((x2 - x3) * u) + ((y2 - y3) * v) + ((z2 - z3) * w) ) ;
    D0 = k * ( (x3 * u) + (y3 * v) + (z3 * w) ) ;

    double Lt      = 0.05 ;  // Length of Taylor series region
    double magDp   = fabs(Dp) ;
    double magDq   = fabs(Dq) ;
    int    TaylorN = 3;
    int    TaylorM = 3;
    
    SPCmplx e_jDp, e_jDq, e_jD0 ;
    SPCmplx jDp, jD0, jDq;
    SPCmplx Ic ;

    e_jD0.r = cosf(D0); e_jD0.i = sinf(D0) ;
    e_jDp.r = cosf(Dp); e_jDp.i = sinf(Dp) ;
    e_jDq.r = cosf(Dq); e_jDq.i = sinf(Dq) ;
      jDp.r = 0;          jDp.i = Dp ;
      jD0.r = 0;          jD0.i = D0 ;
      jDq.r = 0;          jDq.i = Dq ;
    
    SPCmplx G;
    
    if ( magDp < Lt && magDq >= Lt ){  // Case 1
        SPCmplx sum;
        sum.r = sum.i = 0;
        SPCmplx tmp, tmp1;
        SPCmplx C0_over_nplus1 ;
        SPCmplx braces ;
        SPCmplx mult_numer, mult_denom, multiplier ;
        
        // In order to avoid the recursive calculation for G(lambda) which is difficult to
        // implement on a GPU, we loop-unroll the following using order 4 taylor expansions
        //
        
        // n = 0
        //
        G = G0_func(-Dq);
        C0_over_nplus1.r = -C0 ;
        C0_over_nplus1.i = 0 ;
        CMPLX_MULT(e_jDq, G, tmp);
        CMPLX_ADD(C0_over_nplus1, tmp, sum);
        mult_numer.r = 1; mult_numer.i = 0;
        mult_denom.r = 1; mult_denom.i = 0;
        
        // n= 1
        //
        G = G1_func(-Dq);
        C0_over_nplus1.r = -C0 / 2;
        CMPLX_MULT(e_jDq, G, tmp);
        CMPLX_ADD(C0_over_nplus1, tmp, braces);
        tmp = mult_numer ;
        CMPLX_MULT(tmp, jDp, mult_numer);
        CMPLX_SCMULT(1, mult_denom, mult_denom);
        CMPLX_DIV(mult_numer, mult_denom, multiplier);
        CMPLX_MULT(multiplier, braces, tmp);
        CMPLX_ADD(tmp, sum, sum);
        
        // n = 2
        //
        G = G2_func(-Dq);
        C0_over_nplus1.r = -C0 / 3;
        CMPLX_MULT(e_jDq, G, tmp);
        CMPLX_ADD(C0_over_nplus1, tmp, braces);
        tmp = mult_numer ;
        CMPLX_MULT(tmp, jDp, mult_numer);
        CMPLX_SCMULT(2, mult_denom, mult_denom);
        CMPLX_DIV(mult_numer, mult_denom, multiplier);
        CMPLX_MULT(multiplier, braces, tmp);
        CMPLX_ADD(tmp, sum, sum);
        
        // n = 3
        //
        G = G3_func(-Dq);
        C0_over_nplus1.r = -C0 / 4;
        CMPLX_MULT(e_jDq, G, tmp);
        CMPLX_ADD(C0_over_nplus1, tmp, braces);
        tmp = mult_numer ;
        CMPLX_MULT(tmp, jDp, mult_numer);
        CMPLX_SCMULT(3, mult_denom, mult_denom);
        CMPLX_DIV(mult_numer, mult_denom, multiplier);
        CMPLX_MULT(multiplier, braces, tmp);
        CMPLX_ADD(tmp, sum, sum);
        
        // Finished summation
        //
        CMPLX_SCMULT(2*A, e_jD0, tmp);
        CMPLX_DIV(tmp, jDq, tmp1);
        CMPLX_MULT(tmp1, sum, Ic);
        
    }else if ( magDp < Lt && magDq < Lt ){  // Case 2
        
        SPCmplx jDp_a[4] ;
        SPCmplx jDq_a[4] ;
        SPCmplx sum;
        sum.r = sum.i = 0;
        
        jDp_a[0].r = 1 ; jDp_a[0].i = 0 ;
        CMPLX_MULT(jDp_a[0], jDp, jDp_a[1]);
        CMPLX_MULT(jDp_a[1], jDp, jDp_a[2]);
        CMPLX_MULT(jDp_a[2], jDp, jDp_a[3]);
        
        jDq_a[0].r = 1 ; jDq_a[0].i = 0 ;
        CMPLX_MULT(jDq_a[0], jDq, jDq_a[1]);
        CMPLX_MULT(jDq_a[1], jDq, jDq_a[2]);
        CMPLX_MULT(jDq_a[2], jDq, jDq_a[3]);
        
        for(int n=0; n<4; n++){
            for(int m=0; m<4; m++){
                CMPLX_MULT(jDp_a[n],jDq_a[m],tmp);
                CMPLX_SCMULT(1.f/((double)(factorial(m+n+2))), tmp, tmp1);
                CMPLX_ADD(tmp1, sum, sum);
            }
        }
        
        // Finished summation
        //
        CMPLX_SCMULT(2*A, e_jD0, tmp);
        CMPLX_MULT(tmp, sum, Ic);
        
    } else if ( magDp >= Lt && magDq < Lt ){    // Case 3
        
        SPCmplx tmp, tmp1, jDq_pow;
        
        // n = 0
        // C0 = 1, n+1 = 1, jDq^n = 1
        
        sum=G1_func(-Dp);
        
        // n=1
        
        G = G2_func(-Dp);
        CMPLX_SCMULT(0.5, jDq, tmp);
        CMPLX_MULT(tmp, G, tmp1);
        CMPLX_ADD(tmp1, sum, sum);
        
        // n=2
        
        G = G3_func(-Dp);
        CMPLX_MULT(jDq, jDq, jDq_pow);
        CMPLX_SCMULT(0.5, jDq_pow, tmp);
        CMPLX_SCMULT((1.0/3.0), tmp, tmp1);
        CMPLX_MULT(tmp1, G, tmp);
        CMPLX_ADD(tmp, sum, sum);
        
        // n=3
        
        G = G4_func(-Dp);
        tmp = jDq_pow ;
        CMPLX_MULT(jDq, tmp, jDq_pow);
        CMPLX_SCMULT((1.0/6.0), jDq_pow, tmp);
        CMPLX_SCMULT(0.25, tmp, tmp1);
        CMPLX_MULT(tmp1, G, tmp);
        CMPLX_ADD(tmp, sum, sum);
        
        // Finished summation
        //
        CMPLX_SCMULT(2*A, e_jD0, tmp);
        CMPLX_MULT(tmp, e_jDp, tmp1);
        CMPLX_MULT(tmp1, sum, Ic);
        
        
    } else if ( magDp >= Lt && magDq >= Lt && fabs(magDp - magDq) < Lt ){   // Case 4
        
        SPCmplx tmp, tmp1,tmp2, jDp_min_jDq, jD_subs_pow;
        
        CMPLX_SUB(jDp, jDq, jDp_min_jDq);
        jD_subs_pow = jDp_min_jDq ;
        
        // n=0
        G = G0_func(Dq);
        CMPLX_SCMULT(-1, G, G);
        CMPLX_ADD(G, e_jDq, sum);
        
        
        // n=1
        G = G1_func(Dq);
        CMPLX_SCMULT(-1, G, G);
        CMPLX_SCMULT(0.5, e_jDq, tmp);
        CMPLX_ADD(G, tmp, tmp1);
        CMPLX_MULT(jD_subs_pow, tmp1, tmp);
        CMPLX_ADD(tmp, sum, sum);
        
        // n=2
        G = G2_func(Dq);
        CMPLX_SCMULT(-1, G, G);
        CMPLX_SCMULT((1.0/3.0), e_jDq, tmp);
        CMPLX_ADD(G, tmp, tmp1);
        tmp = jD_subs_pow ;
        CMPLX_MULT(tmp, jDp_min_jDq, jD_subs_pow);
        CMPLX_SCMULT(0.5, jD_subs_pow, tmp);
        CMPLX_MULT(tmp, tmp1, tmp2);
        CMPLX_ADD(tmp2, sum, sum);
        
        // n=3
        G=G3_func(Dq);
        CMPLX_SCMULT(-1, G, G);
        CMPLX_SCMULT(0.25, e_jDq, tmp);
        CMPLX_ADD(G, tmp, tmp1);
        tmp = jD_subs_pow ;
        CMPLX_MULT(tmp, jDp_min_jDq, jD_subs_pow);
        CMPLX_SCMULT((1.0/6.0), jD_subs_pow, tmp);
        CMPLX_MULT(tmp, tmp1, tmp2);
        CMPLX_ADD(tmp2, sum, sum);
        
        // End of summation
        CMPLX_SCMULT(2*A, e_jD0, tmp);
        CMPLX_DIV(tmp, jDq, tmp1);
        CMPLX_MULT(tmp1, sum, Ic);
        
        
    } else {
        SPCmplx term1,term2,term3;
        double Dq_minus_Dp = Dq-Dp ;
        CMPLX_SCMULT((C0 / (Dp * Dq_minus_Dp)), e_jDp, term1);
        CMPLX_SCMULT((C0 / (Dq * Dq_minus_Dp)), e_jDq, term2);
        term3.r = C0 / (Dp * Dq) ; term3.i = 0;
        SPCmplx braces ;
        CMPLX_SUB(term1, term2,braces);
        CMPLX_SUB(braces, term3, braces) ;
        Ic.r = cosf(e_jD0) ; Ic.i = sinf(e_jD0) ;
        SPCmplx tmp ;
        CMPLX_SCMULT(2*A, Ic, tmp) ;
        CMPLX_MULT(tmp, braces, Ic);
    }
    
    return Ic ;
}

int factorial(int n)
{
    int c;
    int result = 1;
    
    for (c = 1; c <= n; c++)
        result = result * c;
    
    return result;
}

SPCmplx G4_func(double gamma){
    SPCmplx denominator;
    SPCmplx e_jgamma;
    SPCmplx tmp,ans;
    
    denominator.r = 0 ; denominator.i = gamma ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    
    SPCmplx G3 = G3_func(gamma);
    
    CMPLX_SCMULT(4, G3, G3);
    
    CMPLX_SUB(e_jgamma, G3, tmp);
    CMPLX_DIV(tmp, denominator, ans);
    return ans ;
    
}

SPCmplx G3_func(double gamma){
    SPCmplx denominator;
    SPCmplx e_jgamma;
    SPCmplx tmp,ans;
    
    denominator.r = 0 ; denominator.i = gamma ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    
    SPCmplx G2 = G2_func(gamma);
    
    CMPLX_SCMULT(3, G2, G2);
    
    CMPLX_SUB(e_jgamma, G2, tmp);
    CMPLX_DIV(tmp, denominator, ans);
    return ans ;
    
}

SPCmplx G2_func(double gamma){
    SPCmplx denominator;
    SPCmplx e_jgamma;
    SPCmplx tmp,ans;
    
    denominator.r = 0 ; denominator.i = gamma ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    
    SPCmplx G1 = G1_func(gamma);
    
    CMPLX_SCMULT(2, G1, G1);
    
    CMPLX_SUB(e_jgamma, G1, tmp);
    CMPLX_DIV(tmp, denominator, ans);
    return ans ;
    
}

SPCmplx G1_func(double gamma){
    SPCmplx denominator;
    SPCmplx e_jgamma;
    SPCmplx tmp,ans;
    
    denominator.r = 0 ; denominator.i = gamma ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    
    SPCmplx G0 = G0_func(gamma);
    
    CMPLX_SUB(e_jgamma, G0, tmp);
    CMPLX_DIV(tmp, denominator, ans);
    return ans ;
    
}

SPCmplx G0_func(double gamma){
    SPCmplx ans ;
    SPCmplx one ;
    one.r = 1;
    one.i = 0;
    
    SPCmplx e_jgamma;
    SPCmplx tmp ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    CMPLX_SUB(e_jgamma, one, tmp) ;
    CMPLX_DIV(tmp, e_jgamma, ans) ;
    return ans ;
}
