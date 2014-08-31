//
//  POTriangle.c
//  sarcastic
//
//  Created by Darren on 24/08/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib.h>
#include "sarcastic.h"

SPCmplx G0_func(double gamma);
SPCmplx G1_func(double gamma);
SPCmplx G2_func(double gamma);
SPCmplx G3_func(double gamma);
int factorial(int n);
void matmul(double *A,double *B, double **O, int Ay,int Ax,int By,int Bx);
void reduction(double a[][6],int size,int pivot ,int col) ;
void mat3by3inv(double A[3][3], double O[3][3]);

int POTriangle(TriCoords triCart, Triangle tri, Ray ray, SPVector HitPoint, SPVector ObservationPoint, double lambda, SPCmplx Es){

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // In this first section we'll calculate all the details associated with the triangle
    //
    
    
    // Find normal and area for triangle
    //
    SPVector l1,l3,Av,triN,tmp1,tmp2;
    double a,k,modl1,modl3, Area;
    double baryLambda_1,baryLambda_2,baryLambda_3;
    double Dp,Dq,D0;
    double Cp,Cq,C0;
    
    VECT_SUB(triCart.B, triCart.A, l1) ;
    VECT_SUB(triCart.A, triCart.Cc, l3) ;
    VECT_CROSS(l1, l3, Av) ;
    Area = 0.5 * VECT_MAG(Av) ;
    modl1 = VECT_MOD(l1);
    modl3 = VECT_MOD(l3);
    VECT_SCMULT(Av, (1.0/(modl1*modl3)), triN);

    k = 2 * SIPC_pi / lambda ;
    Cp = Cq = 0;
    C0 = 1;
    
    // find barycentric coords for triangle
    //
    
    
    Dp = k * ( (triCart.A.x - triCart.Cc.x) * baryLambda_1 + (triCart.A.y - triCart.Cc.y) * baryLambda_2 + (triCart.A.z - triCart.Cc.z) * baryLambda_3) ;
    Dq = k * ( (BB.x - CC.x) * baryLambda_1 + (BB.y - CC.y) * baryLambda_2 + (BB.z - Cc.z) * baryLambda_3) ;
    D0 = k * ( CC.x * baryLambda_1 + CC.y * baryLambda_2 + CC.z * baryLambda_3 ) ;
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    // Sort out the direction cosines for this ray, hitpoint and observation point
    //
  
    double u_i,v_i,w_i ;    // Direction cosines for incident ray (ray origin to HitPoint)
    double u_s,v_s,w_s ;    // Direction cosines for scattering ray (HitPoint  to observation point)
    double theta_i, phi_i ; // incident azimuth, elevation
    double theta_s, phi_s ; // Scattering azimuth, elevation
    VECT_SUB(HitPoint, ray.org, tmp1);
    phi_i     = atan2(tmp1.y,tmp1.x);
    theta_i   = atan2(sqrt(tmp1.x*tmp1.x+tmp1.y*tmp1.y),tmp1.z);
    u_i       = sin(theta_i)*cos(phi_i) ;
    v_i       = sin(theta_i)*sin(phi_i) ;
    w_i       = cos(theta_i) ;
    
    
    theta_s = atan2(ObservationPoint.y, ObservationPoint.x);
    phi_s   = atan2(ObservationPoint.z, sqrt(ObservationPoint.x*ObservationPoint.x+ObservationPoint.y*ObservationPoint.y));
    u_s     = cos(phi_s) * cos(theta_s);
    v_s     = cos(phi_s) * sin(theta_s);
    w_s     = sin(phi_s);

    
    // Set up the problem parameters
   
    double k;
    double Z0 ;
    double A ; // Area of triangle
    
    
    // First calculate Ic

    
    
    double Lt = 0.05 ;  // Length of Taylor series region
    double magDp = fabs(Dp) ;
    double magDq = fabs(Dq) ;
    SPCmplx e_jDp, e_jDq, e_jD0 ;
    e_jD0.r = cosf(D0); e_jD0.i = sinf(D0) ;
    e_jDp.r = cosf(Dp); e_jDp.i = sinf(Dp) ;
    e_jDq.r = cosf(Dq); e_jDq.i = sinf(Dq) ;
    SPCmplx jDp, jD0, jDq;
    jDp.r = 0; jDp.i = Dp ;
    jD0.r = 0; jD0.i = D0 ;
    jDq.r = 0; jDq.i = Dq ;
    
    SPCmplx Ic ;
    // The following is designed to calculate Ic - the scattered field
    //
    
    int TaylorN = 3;
    int TaylorM = 3;
    SPCmplx G;
    
    if ( magDp < Lt and magDq >= Lt ){  // Case 1
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
        CMPLX_DIV(mult_numer, mult_denom, multipler);
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
        CMPLX_DIV(mult_numer, mult_denom, multipler);
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
        CMPLX_DIV(mult_numer, mult_denom, multipler);
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
    
    // Now have Ic
    //
    
    // Calculate the surface current from the incident ray
    //
    
    // First calculate the rotation matrices from Global coordinates to local (this triangle) coordinates)
    //
    
    double alpha, beta;
    SPVector zhat;
    VECT_CREATE(0, 0, 1, zhat);
    
    alpha = atan2f(triN.y, triN.x);
    beta  = acos(VECT_DOT(zhat, triN));
    
    double T_dash[9];
    double T_dashdash[9];
    
    T_dash[0] = cos(alpha);
    T_dash[1] = sin(alpha);
    T_dash[2] = 0;
    T_dash[3] = -sin(alpha);
    T_dash[4] = cos(alpha);
    T_dash[5] = 0;
    T_dash[6] = 0;
    T_dash[7] = 0;
    T_dash[8] = 1;
    
    T_dashdash[0] = cos(beta);
    T_dashdash[1] = 0;
    T_dashdash[2] = -sin(beta);
    T_dashdash[3] = 0;
    T_dashdash[4] = 1;
    T_dashdash[5] = 0;
    T_dashdash[6] = sin(beta);
    T_dashdash[7] = 0;
    T_dashdash[8] = cos(beta);
    
    double globalToLocalMat[3][3];
    double localToGlobalMat[3][3];
    
    matmul(T_dashdash, T_dash, &globalToLocalMat, 3, 3, 3, 3);
    mat3by3inv(globalToLocalMat, localToGlobalMat);
    
    double u_dashdash,v_dashdash;
    double uvw_dashdash_arr[3][1] ;
    double uvw_arr[3][1] ;
    uvw_arr[0][0] = u;
    uvw_arr[1][0] = v;
    uvw_arr[2][0] = w;
    
    matmul(globalToLocalMat, uvw_arr, uvw_dashdash_arr, 3, 3, 3, 1);
    
    double u_dashdash, v_dashdash, w_dashdash ;
    u_dashdash = uvw_dashdash_arr[0][0] ;
    v_dashdash = uvw_dashdash_arr[1][0] ;
    w_dashdash = uvw_dashdash_arr[2][0] ;
    
    double sin_theta_dashdash, cos_theta_dash_dash, tan_phi_dashdash ;
    sin_theta_dashdash  = sqrt( u_dashdash*u_dashdash + v_dashdash*v_dashdash) ;
    cos_theta_dash_dash = sqrt(1-(u_dashdash*u_dashdash)+(v_dashdash*v_dashdash)) ;
    tan_phi_dashdash    = atan2f(v_dashdash, u_dashdash);
    
    double Rs, Z0, Gamma_TM, Gamma_TE ;
    Rs = 0 ; // Surface resistivity = 0 = perfect electric conductor
    Z0 = 120 * SIPC_pi ; // Impedence of free space
    
    Gamma_TM = -1.0 * Z0 * cos_theta_dash_dash / (2*Rs + Z0*cos_theta_dash_dash);
    Gamma_TE = -1.0 * Z0 / ( 2.0*Rs*cos_theta_dash_dash + Z0);
    
    // Need E_itheta_dashdash
    //
    
    
    return 0;
}

// 3x3 matrix inversion taken from
// http://www.c4learn.com/c-programs/c-program-to-find-inverse-of-3-x-3.html
// which has a very good graphical description
//
void mat3by3inv(double A[3][3], double O[3][3]){
    double a[3][6];
    int y,x;
    for(y=0;y<3;y++){    // Append Unit Matrix
        for(x=0;x<6;x++){
            if(x<3){
                a[y][x] = A[y][x] ;
            }else if(x==y+3){
                a[y][x]=1;
            }else{
                a[y][x]=0;
            }
        }
    }
    for(i=0;i<3;i++){
        reduction(a,3,i,i);
    }
    
    for (i=0;i<3; i++){
        for(j=0; j<3; j++){
            O[i][j] = a[i][j+3] ;
        }
    }
    return ;
}

//                   a         3        i          i
void reduction(double a[][6],int size,int pivot ,int col)
{
    int i,j;
    double factor;
    
    factor=a[pivot][col];
    
    for(i=0;i<2*size;i++){
        a[pivot][i]/=factor;
    }
    
    for(i=0;i<size;i++){
        if(i!=pivot){
            factor=a[i][col];
            for(j=0;j<2*size;j++)
                a[i][j]=a[i][j]-a[pivot][j]*factor;
        }
    }
    return ;
}

void matmul(double *A,double *B, double **O, int Ay,int Ax,int By,int Bx)
{
    // Ax must be equal to By !!!
    //
    if (Ax != By) {
        printf("Error : Multiplication by incompatable matrices\n");
        exit(1);
    }
    // O must be allocated
    //
    for (i=0;i<Ay;i++){
        for(j=0;j<Bx;j++){
            
            *O[i][j]=0;
            for(int k=0;k<Ax;k++)
                *O[i][j]+= A[i][k]*B[k][j];
        }
    }
    return ;
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
