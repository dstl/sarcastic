//
//  POTriangle.c
//  sarcastic
//
//  Created by Darren on 24/08/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib.h>
#include "POTriangle.h"
#include "matrixMultiplication.h"


int POTriangle(triangle tri, Ray ray, SPVector HitPoint, SPVector obsDir, double lambda, SPCmplx *Es, SPVector *Epol){


    SPVector l1,l3,Av,triN,tmp1,tmp2;
    double a,k,modl1,modl3, Area;
    double baryLambda_1,baryLambda_2,baryLambda_3;

    k = 2 * SIPC_pi / lambda ;
    
    // Sort out the direction cosines for this ray, hitpoint and observation point
    //
    SPVector uvw_ig ; // Direction cosines for incident ray (ray origin to HitPoint)
    SPVector uvw_sg ; // Direction cosines for scattering ray (HitPoint  to observation point)
    double u_i,v_i,w_i ;
    double u_s,v_s,w_s ;
    double theta_i, phi_i ; // incident elevation, azimuth
    double theta_s, phi_s ; // Scattering elevation, azimuth
    
    // Assuming that the direction of the ray has been normalised then
    // the direction cosine is just the componet of direction
    //
    uvw_sg.x = obsDir.x ;
    uvw_sg.y = obsDir.y ;
    uvw_sg.z = obsDir.z ;
    
    
    // Now need to perform two seperate calculations:
    // 1) The surface integral over triangle 'c' called Ic
    // 2) The surface current over the triangle in terms of Jx & Jy
    //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    // 1) Find the value of Ic
    //
    SPCmplx Ic = surfaceIntegral(k, tri, obsPnt, uvw_sg) ;
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    // 2) Find the surface current over the triangle in terms of Jx & Jy. Multiply it
    // by the surfaceIntergral Ic and return the Electric field E (Cmplx amp & phase)
    // also return the polarisation direction of the Efield (in global coordinates)
    //
    Efield(k, tri, ray, Ic, Es, Epol) ;
    

    
    
    return 0;
}


void EField(triangle tri, Ray ray, SPCmplx Ic, SPCmplx Es, SPVector *Epol){
    
    SPVector Eig = ray.pol ;
    
    double theta_il ; // angle to Z axis
    double phi_il ;
    double uvw_ig[3], uvw_il[3]; // Direction cosines for incident Ray
    
    // Assuming that the direction of the ray has been normalised then
    // the direction cosine is just the componet of direction
    //
    uvw_ig[0] = -ray.dir.x;
    uvw_ig[1] = -ray.dir.y;
    uvw_ig[2] = -ray.dir.z;
    
    matmul(tri.globalToLocalMat, uvw_ig, uvw_il, 3, 3, 1, 3);
    
    double sin_theta_il, cos_theta_il, tan_phi_il, sin_phi_il, cos_phi_il ;
    sin_theta_il = sqrt(uvw_il[0] * uvw_il[0] ) ;
    cos_theta_il = sqrt(1 - sin_theta_il*sin_theta_il) ;
    tan_phi_il   = atan(uvw_il[1] / uvw_il[0]) ;
    sin_phi_il   = uvw_il[1] / sin_theta_il ;
    cos_phi_il   = uvw_il[0] / sin_theta_il ;
    
    // Also convert the direction of the E field into the coordinate system
    // of the triangle
    //
    double E_ig[3], E_il[3] ;
    SPVector Eil, theta_l_hat, phi_l_hat, z_l_hat ;

    E_ig[0] = Eig.x ;
    E_ig[1] = Eig.y ;
    E_ig[2] = Eig.z ;
    matmul(tri.globalToLocalMat, E_ig, E_il, 3, 3, 1, 3) ;
    VECT_CREATE(E_il[0], E_il[1], E_il[2], Eil) ;
    VECT_CREATE(0, 0, 1, z_l_hat) ;
    VECT_CROSS (ray.dir, z_l_hat, phi_l_hat) ;
    VECT_CROSS (phi_l_hat, ray.dir, theta_l_hat) ;
    VECT_NORM  (phi_l_hat, phi_l_hat);
    VECT_NORM  (theta_l_hat, theta_l_hat) ;
    double Eiphi_l   = VECT_DOT(Eil, phi_l_hat) ;
    double Eitheta_l = VECT_DOT(Eil, theta_l_hat) ;
    
    // Calculate Gamma_parallel and Gamma_perpendicular
    //
    double Rs, Z0, GamParr, GamPerp ;
    Rs      = tri.Rs ;
    Z0      = 120 * SIPC_pi ; // Impedence of free space
    GamParr = -1.0 * Z0 * cos_theta_il / (2*Rs + Z0*cos_theta_il) ;
    GamPerp = -1.0 * Z0 / ( 2.0*Rs*cos_theta_il + Z0);
    
    // Now calculate Jx_local and Jy_local
    //
    double Jx_l, Jy_l ;
    Jx_l = ((-1.0 * Eitheta_l * cos_phi_il * GamParr / Z0) + (Eiphi_l * sin_phi_il * GamPerp / Z0)) * cos_theta_il ;
    Jy_l = ((-1.0 * Eitheta_l * sin_phi_il * GamParr / Z0) - (Eiphi_l * cos_phi_il * GamPerp / Z0)) * cos_theta_il ;
    
    
    
    return ;
}


SPCmplx surfaceIntegral (double k, triangle tri, SPVector uvw_sg){
    
    double Dp,Dq,D0;
    double Cp,Cq,C0;
    double u = uvw_sg.x ;
    double v = uvw_sg.y ;
    double w = uvw_sg.z ;
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
