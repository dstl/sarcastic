/***************************************************************************
 *
 *       Module:    POTriangle.c
 *      Program:    SARCASTIC
 *   Created by:    Darren on 05/12/2014.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      Program to calculate the Electric field at a receive location caused by an incident ray scattering 
 *      from a triangular facet. The incident ray contains information about its origin and direction and its
 *      polarisation. The programme calculates the induced surface current on the facet and then calculates the
 *      surface integral over the facet. The field contributions are provided to the calling function
 *      in terms or parallel and perpendicular polarisations.
 *
 *  References:
 *
 *      [1] AC. Ludwig ”Computation of radiation patterns involving numerical double integration,”
 *          I€€€ Trans. Antennas Propagat.,vol. AP-16, pp. 767-769, NOV 1968.
 *      [2] R.J. Pogorzelski, ”The Ludwig integration algorithm for triangular subregions,”
 *          Proc. I€€€, vol. 23, pp. 837-838, Apr. 1985.
 *      [3] M.L.X. dos Santos and N.R.Rabelo, "On the Ludwig Integration Algorithm for Triangular Subregions,"
 *          Proc of the IEEE, 74 No.10, pp. 1455-1456, Oct. 1986
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  15/01/2014
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
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
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 *
 ***************************************************************************/

#include <stdio.h>
#include <SIlib/SIlib.h>
#include "POTriangle.h"
#include "matrixMultiplication.h"
#define SURFMAXN 4
#define LT 0.04

SPCmplx G_func(int n, double gamma) ;
int factorial(int n) ;
static const SPVector zz_hat = {0.0, 0.0, 1.0};

void POField(triangle tri, Ray ray, SPVector HitPnt, SPVector obsPnt, double k, SPVector Vdir, SPVector Hdir, SPCmplx *EsV, SPCmplx *EsH){
    
    
    // Some notes on variable terminology used in this function
    //
    // _il = incident ray in local coords
    // _ig = incident ray in global coords
    // _sl = scattered ray in local coords
    // _sg = scattered ray in global coords
    // uvw are direction cosines
    // E is Electric field
    // r is the magnitude of the vector defining the observation point in global coordinates (not the range)
    // J = surface current. J_l[3],J_g[3] are vectors in local and global coords
    // _par refers to parallel or V polarisation. _per refers to perpendicular polarisation
    // Jc_par, Jc_per are the complex surface current components in the parallel and perpendicular polarisation directions
    // E_par and E_per are the parallel and perpendicular E fields at the observation point in V/m
    //
    double Dp,Dq,D0;
    double Cp,Cq,C0;
    double u,v,w, x1,y1,z1,x2,y2,z2,x3,y3,z3,A;
    double r, Lt, magDp, magDq;
    double uvw_ig[3], uvw_il[3], sin_theta_il, cos_theta_il, tan_phi_il, sin_phi_il, cos_phi_il,E_ig[3], E_il[3] ;
    double Dq_minus_Dp, Eiphi_l, Eitheta_l, Rs, GamParr, GamPerp ;
    double J_l[3], J_g[3],J_par, J_per;
    double phs_ig, r_ig;

    SPVector uvw_sg, obsDir, rayDist, Raydir_l ;
    SPVector Eil, Eig, theta_l_hat, phi_l_hat, Jg;

    SPCmplx e_jDp, e_jDq, e_jD0, jDp, jD0, jDq, jDq_pow,jDp_min_jDq,jDp_a[SURFMAXN], jDq_a[SURFMAXN] ;
    SPCmplx C0_over_nplus1, jD_subs_pow, leftPart, G, Ic  ;
    SPCmplx Jc_par, Jc_per, Epar,Eper ;
    SPCmplx jkZ0_o_4PIr, e_jkr;
    SPCmplx sum, tmp, tmp1,t, t1, t2, t3,t4,t5, term1,term2,term3, braces ;

    // r is the magnitude of the vector defining the observation point in global coordinates (not the range)
    //
    VECT_SUB(obsPnt, HitPnt, obsDir);
    r = VECT_MAG(obsDir) ;
    VECT_NORM(obsDir, obsDir);
//    r  = VECT_MAG(obsPnt);
    VECT_SUB(HitPnt, ray.org, rayDist);
    r_ig   = VECT_MAG(rayDist);
    phs_ig = -(k * r_ig) - (SIPC_pi/2.0) ;
    
    // Sort out the direction cosines for this ray, hitpoint and observation point
    // Assuming that the direction of the ray has been normalised then
    // the direction cosine is just the component of direction
    //
    uvw_sg.x = obsDir.x -ray.dir.x ;
    uvw_sg.y = obsDir.y -ray.dir.y ;
    uvw_sg.z = obsDir.z -ray.dir.z ;
 
    // Calculate surface integral using Ludwig's integration algorithm [1], (modified for
    // triangular subregions by Pogorzelski [2] and then modified again for barycentric (simplex)
    // coordinates by Dos Santos [3].
    //
    u  = uvw_sg.x ;
    v  = uvw_sg.y ;
    w  = uvw_sg.z ;
    x1 = tri.AA.x ;
    y1 = tri.AA.y ;
    z1 = tri.AA.z ;
    x2 = tri.BB.x ;
    y2 = tri.BB.y ;
    z2 = tri.BB.z ;
    x3 = tri.CC.x ;
    y3 = tri.CC.y ;
    z3 = tri.CC.z ;
    A  = tri.area ;
    
    Cp = Cq = 0;
    C0 = 1;
    Dp = k * ( ((x1 - x3) * u) + ((y1 - y3) * v) + ((z1 - z3) * w) ) ;
    Dq = k * ( ((x2 - x3) * u) + ((y2 - y3) * v) + ((z2 - z3) * w) ) ;
    D0 = k * ( (x3 * u) + (y3 * v) + (z3 * w) ) ;
    
    Lt      = LT ;          // Length of Taylor series region
    magDp   = fabs(Dp) ;
    magDq   = fabs(Dq) ;
    e_jD0.r = cosf(D0) ; e_jD0.i = sinf(D0) ;
    e_jDp.r = cosf(Dp) ; e_jDp.i = sinf(Dp) ;
    e_jDq.r = cosf(Dq) ; e_jDq.i = sinf(Dq) ;
    jDp.r   = 0;           jDp.i = Dp ;
    jD0.r   = 0;           jD0.i = D0 ;
    jDq.r   = 0;           jDq.i = Dq ;
    
    if ( magDp < Lt && magDq >= Lt ){  // Case 1
        sum.r = sum.i = 0;
        
        // Set up zeroeth term of sum
        //
        jD_subs_pow.r = 1.0;
        jD_subs_pow.i = 0.0 ;
        t1 = G_func(0, -Dq) ;
        CMPLX_MULT(e_jDq, t1, t) ;
        C0_over_nplus1.r = -C0 ;
        C0_over_nplus1.i = 0 ;
        CMPLX_ADD(C0_over_nplus1, t, sum);
        
        // Now sum over n
        //
        for (int n=1; n<=SURFMAXN; n++) {
            CMPLX_MULT(jD_subs_pow, jDp, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func(n, -Dq) ;
            CMPLX_MULT(e_jDq, t1, t) ;
            C0_over_nplus1.r = -C0 / (n+1) ;
            C0_over_nplus1.i = 0.0 ;
            CMPLX_ADD(C0_over_nplus1, t, t1) ;
            CMPLX_MULT(leftPart, t1, t) ;
            CMPLX_ADD(sum, t, sum) ;
        }
        
        // Finished summation
        //
        CMPLX_SCMULT(2*A, e_jD0, t);
        CMPLX_DIV(t, jDq, t1);
        CMPLX_MULT(t1, sum, Ic);
        
    }else if ( magDp < Lt && magDq < Lt ){  // Case 2
        
        sum.r = sum.i = 0;
        jDp_a[0].r = 1 ; jDp_a[0].i = 0 ;
        jDq_a[0].r = 1 ; jDq_a[0].i = 0 ;
        
        for (int i=1; i<SURFMAXN; i++) {
            CMPLX_MULT(jDp_a[i-1], jDp, jDp_a[i]);
            CMPLX_MULT(jDq_a[i-1], jDq, jDq_a[i]);
        }
        for(int n=0; n<SURFMAXN; n++){
            for(int m=0; m<SURFMAXN; m++){
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
        
        sum.r = sum.i = 0;
        jDq_pow.r = 1.0; jDq_pow.i = 0.0 ;
        G = G_func(1, -Dp) ;
        sum = G ;
        
        for(int n=1; n<SURFMAXN; n++){
            G = G_func(n+1, -Dp) ;
            CMPLX_SCMULT(1.0/(double)(n+1), G, t);
            CMPLX_MULT(jDq_pow, jDq, t1);
            jDq_pow = t1 ;
            CMPLX_SCMULT(1/(double)factorial(n), jDq_pow, t1);
            CMPLX_MULT(t1, t, t2);
            CMPLX_ADD(sum, t2, sum);
        }
        
        // Finished summation
        //
        CMPLX_SCMULT(2*A, e_jD0, t);
        CMPLX_MULT(t, e_jDp, t1);
        CMPLX_MULT(t1, sum, Ic);
        
    } else if ( magDp >= Lt && magDq >= Lt && fabs(Dp - Dq) < Lt ){   // Case 4
       
        sum.r = sum.i = 0;
        // Set up zeroeth term of sum
        //
        CMPLX_SUB(jDp, jDq, jDp_min_jDq);
        jD_subs_pow.r = 1.0; jD_subs_pow.i = 0.0 ;
        // C0 = 1
        //
        t1 = G_func(0, Dq) ;
        CMPLX_SCMULT(-1, t1, t1);
        CMPLX_ADD(t1, e_jDp, sum);
        
        // Now sum over n
        //
        for (int n=1; n<=SURFMAXN; n++) {
            CMPLX_MULT(jD_subs_pow, jDp_min_jDq, t);
            jD_subs_pow = t ;
            CMPLX_SCMULT(1.0/(double)factorial(n), t, leftPart) ;
            t1 = G_func(n, Dq) ;
            CMPLX_SCMULT(-1, t1, t2);
            CMPLX_SCMULT(1.0/(n+1.0), e_jDp, t3);
            CMPLX_ADD(t2, t3, t4);
            CMPLX_MULT(leftPart, t4, t5);
            CMPLX_ADD(t5, sum, sum) ;
        }
        
        // End of summation
        CMPLX_SCMULT(2*A, e_jD0, t);
        CMPLX_DIV(t, jDq, t1);
        CMPLX_MULT(t1, sum, Ic);
        
    } else {
        Dq_minus_Dp = Dq-Dp ;
        CMPLX_SCMULT((C0 / (Dp * Dq_minus_Dp)), e_jDp, term1);
        CMPLX_SCMULT((C0 / (Dq * Dq_minus_Dp)), e_jDq, term2);
        term3.r = C0 / (Dp * Dq) ; term3.i = 0;
        CMPLX_SUB(term1, term2, braces);
        CMPLX_SUB(braces, term3, braces) ;
        CMPLX_SCMULT(2*A, e_jD0, t) ;
        CMPLX_MULT(t, braces, Ic);
    }
    
    // At this point we have integrated the complex field over the triangular facet and have the
    // value stored in the complex number Ic
    //
    VECT_SCMULT(ray.pol, sqrt(ray.pow), Eig) ;
    
    // Assuming that the direction of the ray has been normalised then
    // the direction cosine is just the component of direction
    //
    uvw_ig[0] = -ray.dir.x;
    uvw_ig[1] = -ray.dir.y;
    uvw_ig[2] = -ray.dir.z;
    
    matmul(tri.globalToLocalMat, uvw_ig, uvw_il, 3, 3, 1, 3);
    
    sin_theta_il = sqrt(uvw_il[0]*uvw_il[0] +  uvw_il[1] * uvw_il[1] ) ;
    cos_theta_il = sqrt(1 - sin_theta_il*sin_theta_il) ;
    if(sin_theta_il <= 0.0001){
        tan_phi_il = 0.0 ;
        sin_phi_il = 1.0 ;
        cos_phi_il = 0.0 ;
    }else{
        tan_phi_il   = atan2(uvw_il[1] , uvw_il[0]) ;
        sin_phi_il   = uvw_il[1] / sin_theta_il ;
        cos_phi_il   = uvw_il[0] / sin_theta_il ;
    }
    
    // Also convert the direction of the E field into the coordinate system
    // of the triangle
    //
    VECT_CREATE(uvw_il[0], uvw_il[1], uvw_il[2], Raydir_l) ;
    E_ig[0] = Eig.x ;
    E_ig[1] = Eig.y ;
    E_ig[2] = Eig.z ;
    matmul(tri.globalToLocalMat, E_ig, E_il, 3, 3, 1, 3) ;
    VECT_CREATE(E_il[0], E_il[1], E_il[2], Eil) ;
    if(fabs(VECT_DOT(zz_hat, Raydir_l)) >=0.99){
        VECT_CREATE(1, 0, 0, phi_l_hat) ;
    }else{
        VECT_CROSS (zz_hat, Raydir_l, phi_l_hat) ;
    }
    VECT_CROSS(Raydir_l, phi_l_hat, theta_l_hat) ;
    VECT_NORM(phi_l_hat, phi_l_hat) ;
    VECT_NORM(theta_l_hat, theta_l_hat) ;
    Eiphi_l   = VECT_DOT(Eil, phi_l_hat) ;
    Eitheta_l = VECT_DOT(Eil, theta_l_hat) ;
    
    // Calculate Gamma_parallel and Gamma_perpendicular
    //
    Rs      = materialProperties[tri.matId].resistivity ;
    GamParr = -1.0 * Z0 * cos_theta_il / (2*Rs + Z0*cos_theta_il) ;
    GamPerp = -1.0 * Z0 / ( 2.0*Rs*cos_theta_il + Z0);
    
    // Now calculate Jx_local and Jy_local and convert them to global
    // vectors.
    //
    J_l[0] = ((-1.0 * Eitheta_l * cos_phi_il * GamParr / Z0) + (Eiphi_l * sin_phi_il * GamPerp / Z0)) * cos_theta_il ;
    J_l[1] = ((-1.0 * Eitheta_l * sin_phi_il * GamParr / Z0) - (Eiphi_l * cos_phi_il * GamPerp / Z0)) * cos_theta_il ;
    J_l[2] = 0.0  ;
    matmul(tri.localToGlobalMat, J_l, J_g, 3, 3, 1, 3) ;
    VECT_CREATE(J_g[0], J_g[1], J_g[2], Jg) ;
    
    // Find the component of the surface current in the parallel polarisation
    // direction
    //
    J_par = VECT_DOT(Vdir, Jg);
    J_per = VECT_DOT(Hdir, Jg);
    CMPLX_F_MAKE(J_par*cos(phs_ig), J_par*sin(phs_ig), Jc_par);
    CMPLX_F_MAKE(J_per*cos(phs_ig), J_per*sin(phs_ig), Jc_per);
//    printf("Vdir: %f,%f,%f\n",Vdir.x,Vdir.y,Vdir.z);
//    printf("Jg : %f,%f,%f\n",Jg.x,Jg.y,Jg.z);
//    printf("J_par : %f\n",J_par);
//    printf("phs_ig : %f\n",phs_ig);
    
    // Work out scaler (complex) component of E field
    //
    CMPLX_F_MAKE(0, -k*Z0/(4*SIPC_pi*r), jkZ0_o_4PIr) ;
    e_jkr.r = cos(-k*r);
    e_jkr.i = sin(-k*r);
    
//    printf("phase(Jc_par)      : %8.5f (%6.1f deg) (%e, %e)\n",CMPLX_PHASE(Jc_par), RAD2DEG(CMPLX_PHASE(Jc_par)),Jc_par.r,Jc_par.i);
//    printf("phase(e_jkr)       : %8.5f (%6.1f deg) (total: %f)\n",CMPLX_PHASE(e_jkr), RAD2DEG(CMPLX_PHASE(e_jkr)), -k*r);
//    printf("phase(jkZ0_o_4PIr) : %8.5f (%6.1f deg) (%e, %e)\n",CMPLX_PHASE(jkZ0_o_4PIr), RAD2DEG(CMPLX_PHASE(jkZ0_o_4PIr)),jkZ0_o_4PIr.r,jkZ0_o_4PIr.i);
//    printf("phase(Ic)          : %8.5f (%6.1f deg) (%e, %e)\n",CMPLX_PHASE(Ic), RAD2DEG(CMPLX_PHASE(Ic)),Ic.r,Ic.i);

    CMPLX_MULT(jkZ0_o_4PIr, e_jkr, tmp1);
    CMPLX_MULT(tmp1, Ic, tmp);
    CMPLX_MULT(Jc_par, tmp, Epar);
    CMPLX_MULT(Jc_per, tmp, Eper);
//    printf("Ic:%f  J_par:%f, Jc_per:%f\n",CMPLX_PHASE(Ic),CMPLX_PHASE(Jc_par),CMPLX_PHASE(Jc_per));
    EsV->r = Epar.r ; EsV->i = Epar.i ;
    EsH->r = Eper.r ; EsH->i = Eper.i ;
    return ;

}

// Function to calculate n!
//
int factorial(int n)
{
    int c;
    int result = 1;
    
    for (c = 1; c <= n; c++)
        result = result * c;
    
    return result;
}

// Recursive function to calculate
//    G(n,gamma) = e^jgamma - nG(n-1,gamma)
//                -------------------------
//                        jgamma
//
SPCmplx G_func(int n, double gamma){
    SPCmplx e_jgamma, jgamma, G, t1, t2, ans ;
    e_jgamma.r = cos(gamma) ;
    e_jgamma.i = sin(gamma) ;
    jgamma.r = 0 ; jgamma.i = gamma ;
    
    if (n==0) {
        SPCmplx one ;
        one.r = 1;
        one.i = 0;
        CMPLX_SUB(e_jgamma, one, t1) ;
        CMPLX_DIV(t1, jgamma, ans) ;

    }else{
        G = G_func(n-1, gamma) ;
        CMPLX_SCMULT(n, G, t1);
        CMPLX_SUB(e_jgamma, t1, t2) ;
        CMPLX_DIV(t2, jgamma, ans) ;
    }
    return ans ;
}
