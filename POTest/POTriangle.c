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

SPCmplx G0_func(double gamma) ;
SPCmplx G1_func(double gamma) ;
SPCmplx G2_func(double gamma) ;
SPCmplx G3_func(double gamma) ;
SPCmplx G4_func(double gamma) ;
int factorial(int n) ;


void POTriangle(triangle tri, Ray ray, SPVector HitPoint, SPVector obsPnt, double lambda, SPCmplx *EsV, SPCmplx *EsH){

    double k;
    
    k = 2 * SIPC_pi / lambda ;
    
    // Sort out the direction cosines for this ray, hitpoint and observation point
    //
    SPVector uvw_sg ; // Direction cosines for scattering ray (HitPoint  to observation point)
    
    SPVector obsDir ;
    VECT_SUB(obsPnt, HitPoint, obsDir);
    double r = VECT_MAG(obsDir);
    VECT_SCMULT(obsDir, 1/r, obsDir) ;
    
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
    SPCmplx Ic = surfaceIntegral(k, tri, uvw_sg) ;
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //
    // 2) Find the surface current over the triangle in terms of Jx & Jy. Multiply it
    // by the surfaceIntergral Ic and return the Electric field E (Cmplx amp & phase)
    // also return the polarisation direction of the Efield (in global coordinates)
    //
    SPCmplx Es_parr, Es_perp ;
    // Define unit vectors for V & H fields
    //
    SPVector Es_parrdir, Es_perpdir, z_hat ;
    VECT_CREATE(0, 0, 1, z_hat) ;
    VECT_CROSS(z_hat, obsDir, Es_perpdir) ;
    VECT_CROSS(obsDir, Es_perpdir, Es_parrdir) ;
    VECT_NORM(Es_parrdir, Es_parrdir) ;
    VECT_NORM(Es_perpdir, Es_perpdir) ;
    
    EField(k, r, tri, ray, Ic, Es_parrdir, Es_perpdir, &Es_parr, &Es_perp) ;
    *EsV = Es_parr ;
    *EsH = Es_perp ;
    return ;
}


void EField(double k, double r, triangle tri, Ray ray, SPCmplx Ic, SPVector Es_parrdir, SPVector Es_perpdir, SPCmplx *Es_parr, SPCmplx *Es_perp ){
    
    SPVector Eig = ray.pol ;
    double uvw_ig[3], uvw_il[3]; // Direction cosines for incident Ray
    
    // Assuming that the direction of the ray has been normalised then
    // the direction cosine is just the component of direction
    //
    uvw_ig[0] = -ray.dir.x;
    uvw_ig[1] = -ray.dir.y;
    uvw_ig[2] = -ray.dir.z;
    
    matmul(tri.globalToLocalMat, uvw_ig, uvw_il, 3, 3, 1, 3);
    
    double sin_theta_il, cos_theta_il, tan_phi_il, sin_phi_il, cos_phi_il ;
    sin_theta_il = sqrt(uvw_il[0]*uvw_il[0] +  uvw_il[1] * uvw_il[1] ) ;
    cos_theta_il = sqrt(1 - sin_theta_il*sin_theta_il) ;
    tan_phi_il   = atan(uvw_il[1] / uvw_il[0]) ;
    sin_phi_il   = uvw_il[1] / sin_theta_il ;
    cos_phi_il   = uvw_il[0] / sin_theta_il ;
    
    // Also convert the direction of the E field into the coordinate system
    // of the triangle
    //
    double E_ig[3], E_il[3] ;
    SPVector Eil, theta_l_hat, phi_l_hat, z_hat ;

    E_ig[0] = Eig.x ;
    E_ig[1] = Eig.y ;
    E_ig[2] = Eig.z ;
    matmul(tri.globalToLocalMat, E_ig, E_il, 3, 3, 1, 3) ;
    VECT_CREATE(E_il[0], E_il[1], E_il[2], Eil) ;
    VECT_CREATE(0, 0, 1, z_hat) ;
    VECT_CROSS (ray.dir, z_hat, phi_l_hat) ;
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
    
    // Convert the current vector to global coordinates
    //
    double J_l[3], J_g[3] ;
    J_l[0] = Jx_l ;
    J_l[1] = Jy_l ;
    J_l[2] = 0.0  ;
    matmul(tri.localToGlobalMat, J_l, J_g, 3, 3, 1, 3) ;
    SPVector Jg;
    VECT_CREATE(J_g[0], J_g[1], J_g[2], Jg) ;
    double A = VECT_DOT(Es_parrdir, Jg);
    double B = VECT_DOT(Es_perpdir, Jg);
    
    // Work out scaler (complex) component of E field
    //
    SPCmplx jkZ0_o_4PIr, tmp1, scaler, e_jkr;
    CMPLX_F_MAKE(0, -k*Z0/(4*SIPC_pi*r), jkZ0_o_4PIr) ;
    e_jkr.r = cos(-k*r);
    e_jkr.i = sin(-k*r);
    CMPLX_MULT(jkZ0_o_4PIr, e_jkr, tmp1);
    CMPLX_MULT(tmp1, Ic, scaler);
    SPCmplx par,per ;
    CMPLX_SCMULT(A, scaler, par) ;
    CMPLX_SCMULT(B, scaler, per) ;
    Es_parr->r = par.r ; Es_parr->i = par.i ;
    Es_perp->r = per.r ; Es_perp->i = per.i ;
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
        SPCmplx t, t1;
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
        CMPLX_MULT(e_jDq, G, t);
        CMPLX_ADD(C0_over_nplus1, t, sum);
        mult_numer.r = 1; mult_numer.i = 0;
        mult_denom.r = 1; mult_denom.i = 0;
        
        // n= 1
        //
        G = G1_func(-Dq);
        C0_over_nplus1.r = -C0 / 2;
        CMPLX_MULT(e_jDq, G, t);
        CMPLX_ADD(C0_over_nplus1, t, braces);
        t = mult_numer ;
        CMPLX_MULT(t, jDp, mult_numer);
        CMPLX_SCMULT(1, mult_denom, mult_denom);
        CMPLX_DIV(mult_numer, mult_denom, multiplier);
        CMPLX_MULT(multiplier, braces, t);
        CMPLX_ADD(t, sum, sum);
        
        // n = 2
        //
        G = G2_func(-Dq);
        C0_over_nplus1.r = -C0 / 3;
        CMPLX_MULT(e_jDq, G, t);
        CMPLX_ADD(C0_over_nplus1, t, braces);
        t = mult_numer ;
        CMPLX_MULT(t, jDp, mult_numer);
        CMPLX_SCMULT(2, mult_denom, mult_denom);
        CMPLX_DIV(mult_numer, mult_denom, multiplier);
        CMPLX_MULT(multiplier, braces, t);
        CMPLX_ADD(t, sum, sum);
        
        // n = 3
        //
        G = G3_func(-Dq);
        C0_over_nplus1.r = -C0 / 4;
        CMPLX_MULT(e_jDq, G, t);
        CMPLX_ADD(C0_over_nplus1, t, braces);
        t = mult_numer ;
        CMPLX_MULT(t, jDp, mult_numer);
        CMPLX_SCMULT(3, mult_denom, mult_denom);
        CMPLX_DIV(mult_numer, mult_denom, multiplier);
        CMPLX_MULT(multiplier, braces, t);
        CMPLX_ADD(t, sum, sum);
        
        // Finished summation
        //
        CMPLX_SCMULT(2*A, e_jD0, t);
        CMPLX_DIV(t, jDq, t1);
        CMPLX_MULT(t1, sum, Ic);
        
    }else if ( magDp < Lt && magDq < Lt ){  // Case 2
        
        SPCmplx jDp_a[4] ;
        SPCmplx jDq_a[4] ;
        SPCmplx sum; sum.r = sum.i = 0;
        SPCmplx tmp, tmp1;
        
        
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
        SPCmplx sum; sum.r = sum.i = 0;
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
        
        
    } else if ( magDp >= Lt && magDq >= Lt && fabs(Dp - Dq) < Lt ){   // Case 4
        
        SPCmplx t, t1,t2, jDp_min_jDq, jD_subs_pow;
        SPCmplx sum; sum.r = sum.i = 0;

        CMPLX_SUB(jDp, jDq, jDp_min_jDq);
        jD_subs_pow = jDp_min_jDq ;
        
        // n=0
        G = G0_func(Dq);
        CMPLX_SCMULT(-1, G, G);
        CMPLX_ADD(G, e_jDq, sum);
        
        
        // n=1
        G = G1_func(Dq);
        CMPLX_SCMULT(-1, G, G);
        CMPLX_SCMULT(0.5, e_jDq, t);
        CMPLX_ADD(G, t, t1);
        CMPLX_MULT(jD_subs_pow, t1, t);
        CMPLX_ADD(t, sum, sum);
        
        // n=2
        G = G2_func(Dq);
        CMPLX_SCMULT(-1, G, G);
        CMPLX_SCMULT((1.0/3.0), e_jDq, t);
        CMPLX_ADD(G, t, t1);
        t = jD_subs_pow ;
        CMPLX_MULT(t, jDp_min_jDq, jD_subs_pow);
        CMPLX_SCMULT(0.5, jD_subs_pow, t);
        CMPLX_MULT(t, t1, t2);
        CMPLX_ADD(t2, sum, sum);
        
        // n=3
        G=G3_func(Dq);
        CMPLX_SCMULT(-1, G, G);
        CMPLX_SCMULT(0.25, e_jDq, t);
        CMPLX_ADD(G, t, t1);
        t = jD_subs_pow ;
        CMPLX_MULT(t, jDp_min_jDq, jD_subs_pow);
        CMPLX_SCMULT((1.0/6.0), jD_subs_pow, t);
        CMPLX_MULT(t, t1, t2);
        CMPLX_ADD(t2, sum, sum);
        
        // End of summation
        CMPLX_SCMULT(2*A, e_jD0, t);
        CMPLX_DIV(t, jDq, t1);
        CMPLX_MULT(t1, sum, Ic);
        
        
    } else {
        SPCmplx term1,term2,term3;
        double Dq_minus_Dp = Dq-Dp ;
        CMPLX_SCMULT((C0 / (Dp * Dq_minus_Dp)), e_jDp, term1);
        CMPLX_SCMULT((C0 / (Dq * Dq_minus_Dp)), e_jDq, term2);
        term3.r = C0 / (Dp * Dq) ; term3.i = 0;
        SPCmplx braces ;
        CMPLX_SUB(term1, term2,braces);
        CMPLX_SUB(braces, term3, braces) ;
        SPCmplx t ;
        CMPLX_SCMULT(2*A, e_jD0, t) ;
        CMPLX_MULT(t, braces, Ic);
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
    SPCmplx t,ans;
    
    denominator.r = 0 ; denominator.i = gamma ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    
    SPCmplx G3 = G3_func(gamma);
    
    CMPLX_SCMULT(4, G3, G3);
    
    CMPLX_SUB(e_jgamma, G3, t);
    CMPLX_DIV(t, denominator, ans);
    return ans ;
    
}

SPCmplx G3_func(double gamma){
    SPCmplx denominator;
    SPCmplx e_jgamma;
    SPCmplx t,ans;
    
    denominator.r = 0 ; denominator.i = gamma ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    
    SPCmplx G2 = G2_func(gamma);
    
    CMPLX_SCMULT(3, G2, G2);
    
    CMPLX_SUB(e_jgamma, G2, t);
    CMPLX_DIV(t, denominator, ans);
    return ans ;
    
}

SPCmplx G2_func(double gamma){
    SPCmplx denominator;
    SPCmplx e_jgamma;
    SPCmplx t,ans;
    
    denominator.r = 0 ; denominator.i = gamma ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    
    SPCmplx G1 = G1_func(gamma);
    
    CMPLX_SCMULT(2, G1, G1);
    
    CMPLX_SUB(e_jgamma, G1, t);
    CMPLX_DIV(t, denominator, ans);
    return ans ;
    
}

SPCmplx G1_func(double gamma){
    SPCmplx denominator;
    SPCmplx e_jgamma;
    SPCmplx t,ans;
    
    denominator.r = 0 ; denominator.i = gamma ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    
    SPCmplx G0 = G0_func(gamma);
    
    CMPLX_SUB(e_jgamma, G0, t);
    CMPLX_DIV(t, denominator, ans);
    return ans ;
    
}

SPCmplx G0_func(double gamma){
    SPCmplx ans ;
    SPCmplx one ;
    one.r = 1;
    one.i = 0;
    
    SPCmplx e_jgamma;
    SPCmplx t ;
    e_jgamma.r = cosf(gamma) ;
    e_jgamma.i = sinf(gamma) ;
    CMPLX_SUB(e_jgamma, one, t) ;
    CMPLX_DIV(t, e_jgamma, ans) ;
    return ans ;
}
