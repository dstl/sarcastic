/***************************************************************************
 *
 *       Module:    POCalculation.c
 *      Program:    SARCASTIC
 *   Created by:    Darren on 03/08/2014.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  09/03/2014
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

#include "sarcastic.h"

/*
    GtoLMatrix - Global coordinates to local coordin
    ri_bar Vector incident Ray
    p   hitpoint on facet
    rs_bar   Vector to receive location
    Es - amplitude of E field at receiver
    Rs - resistivity
    a - ray width
    b - ray height

*/

void setupTriangle(TriCoords tri, SPVector *normal, double *matrix_gtol, double *matrix_ltog, double *Rs);
void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By);
void mat3by3inv(double *A, double *O);
void reduction(double a[][6],int size,int pivot ,int col);

void POCalculation(TriCoords tri,
                   Ray Ri_bar,      // Incident Ray
                   SPVector p,      // Global coordinates of scattered point
                   SPVector RxPnt,  // Observation point
                   SPCmplx Ei,      // incident E field at point P
                   double a,        // Ray size x-dimension
                   double b,        // Ray size y-dimension
                   SPCmplx *Es      // Scattered E field

){
    
    double Beta ;
    double lambda = 0.04 ;
    Beta = 2*SIPC_pi / lambda ;
    
    // setup triangle associated parameters such as the normal and rotation matrix
    // These can be performed when the triangles are built
    //
    SPVector triNormal ;
    double Glob2LocMat[9], Loc2GlobMat[9] ;
    double Rs ; // resistivity
    setupTriangle(tri, &triNormal, Glob2LocMat, Loc2GlobMat, &Rs);
    
    
    // Calculate the Observation direction and range
    //
    SPVector Rscat, RscatDir ; // Scattered ray, Scattered ray direction (normalised)
    VECT_SUB(RxPnt, p, Rscat);
    double r = VECT_MAG(Rscat);
    VECT_NORM(Rscat, RscatDir);
    

    // Convert the incident ray into the coordinate system of the triangle
    // use _l for coordinates in local frame of facet and _g or nothing for global
    // world coordinates
    //
    
    double theta_il, theta_sl ; // angle to Z axis
    double phi_il, phi_sl ;
    double uvw_ig[3], uvw_il[3], uvw_sg[3], uvw_sl[3] ; // Direction cosines for incident and scattered Ray
    
    // Assuming that the direction of the ray has been normalised
    //
    uvw_ig[0] = -Ri_bar.dir.x;
    uvw_ig[1] = -Ri_bar.dir.y;
    uvw_ig[2] = -Ri_bar.dir.z;
    
    matmul(Glob2LocMat, uvw_ig, uvw_il, 3, 3, 1, 3);
    
    theta_il = acos(uvw_il[2]) ;
    phi_il   = atan2(uvw_il[1], uvw_il[0]);
    
    // Now calculate the direction to the viewpoint (along the scattered ray)
    //
    uvw_sg[0] = RscatDir.x ;
    uvw_sg[1] = RscatDir.y ;
    uvw_sg[2] = RscatDir.z ;
    
    matmul(Glob2LocMat, uvw_sg, uvw_sl, 3, 3, 1, 3);
    
    theta_sl = acos(uvw_sl[2]) ;
    phi_sl   = atan2(uvw_sl[1],uvw_sl[0]);
    
    // Calculate X and Y, SINCX and SINCY
    //
    double X = 0.5 * Beta*(uvw_sl[0]+uvw_il[0])*a ;
    double Y = 0.5 * Beta*(uvw_sl[1]+uvw_il[1])*b ;
    double sincX = (X != 0 ? sin(X)/X : 1);
    double sincY = (Y != 0 ? sin(Y)/Y : 1);
    
    // Calculate C0
    //
    double Z0 = 120 * SIPC_pi ; // Impedence of free space
    SPCmplx tmp,tmp1 ;
    SPCmplx j ;
    SPCmplx C0 ;
    SPCmplxPol e_jbeta_r;
    SPCmplx e_jBeta_r_cart ;
    
    CMPLX_F_MAKE(0, 1, j);
    e_jbeta_r.pabs = 1.0f;
    e_jbeta_r.parg = -1*Beta*r ;
    CMPLX_POL2CART(e_jbeta_r, e_jBeta_r_cart);

    CMPLX_SCMULT((Z0*Beta*a*b*sincX*sincY / (2*SIPC_pi*r)), j, tmp);
    CMPLX_MULT(tmp, e_jBeta_r_cart, C0);
    
    // Now calculate Es
    //
    double rtmp;
    double costheta_il = uvw_il[2] ;
    double costheta_sl = uvw_sl[2] ;
//    rtmp = (-1 * costheta_il / ((Z0 * costheta_il) + (2*Rs))) * ((costheta_sl * sin(phi_sl)) + cos(phi_sl)) ;
//    rtmp = (-1 * costheta_il / ((Z0 * costheta_il) + (2*Rs))) * ((costheta_sl * cos(phi_sl-phi_il)) + sin(phi_il-phi_sl)) ;
    SPCmplx E_theta_theta;
    rtmp = (-1 * costheta_il / ((Z0 * costheta_il) + (2*Rs))) * (costheta_sl*cos(phi_sl-phi_il)) ; // E_theta_theta
    // Assume Ei is theta polarised (ie in the plane of incidence)
    //
    CMPLX_SCMULT(rtmp, C0, tmp);
    CMPLX_MULT(tmp, Ei, E_theta_theta);
    
    SPCmplx E_phi_theta;
    rtmp = (-1 * costheta_il / ((Z0 * costheta_il) + (2*Rs))) * (sin(phi_il-phi_sl)) ; // E_phi_theta
    // Assume Ei is theta polarised (ie in the plane of incidence)
    //
    CMPLX_SCMULT(rtmp, C0, tmp);
    CMPLX_MULT(tmp, Ei, E_phi_theta);

    // To find Es we need to find the magnitude of the Es vector
    // from the parallel E field (E_theta) and the perpendicular E_s
    // field (E_phi)
    ///
    
    SPCmplx Es_bar_phi, Es_bar_theta; // r,phi,theta
    
//    SPVector V_E_theta_theta_l, V_E_phi_theta_l ;
//    VECT_CROSS(RscatDir, Z_DIR, V_E_phi_theta_l);
//    VECT_NORM(V_E_phi_theta_l, V_E_phi_theta_l);
//    VECT_CROSS(V_E_phi_theta_l, RscatDir, V_E_theta_theta_l);
//    VECT_NORM(V_E_theta_theta_l, V_E_theta_theta_l);
    
    // Now rotate the local vectors to global vectors
    //
    
//    Es->r = E_phi_theta.r ;
//    Es->i = E_phi_theta.i ;
    Es->r = E_theta_theta.r ;
    Es->i = E_theta_theta.i ;
    
    return ;
}

void setupTriangle(TriCoords tri, SPVector *normal, double *matrix_gtol, double *matrix_ltog, double *Rs){
    
    // Calculate triangle Normal
    //
    SPVector l1,l3,Av,triN;
    double Area;
    VECT_SUB(tri.B, tri.A, l1) ;
    VECT_SUB(tri.A, tri.Cc, l3) ;
    VECT_CROSS(l3, l1, Av) ;
    Area = 0.5 * VECT_MAG(Av) ;
    VECT_NORM(Av, triN);
    normal->x = triN.x ; normal->y = triN.y; normal->z = triN.z ;
    
    // Calculate rotation matrices
    //
    double T_dash[9];
    double T_dashdash[9];
    double alpha,beta ;
    SPVector zhat;
    VECT_CREATE(0, 0, 1, zhat);
    
    alpha = atan2(triN.y, triN.x);
    beta  = acos(VECT_DOT(zhat, triN));
    
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
    
    // Calculate global to local coordinate transformation matrix
    //
    matmul(T_dashdash, T_dash, matrix_gtol, 3, 3, 3, 3);
    mat3by3inv(matrix_gtol, matrix_ltog);
    
    *Rs = 0; // Perfect Electrical Conductor
    
    return ;


}

// Matrix multiplication for a N size array
// Note : in order to remove the test for incompatable dimensions (for speed)
// it is essential that Ax = By
// Also the output matrix *O must have been allocated to the correct dimensions
//  double *A - matrix A
//  double *B - matrix B. By must equal Ax
//  double *O - output matrix. Dimensions must be Ay x Bx
//  int Ax  - Number of columns in A
//  int Ay  - Number of rows in A
//  int Bx  - number of columns in B
//  int By  - number of rows in B
//
void matmul(double *A,double *B, double *O, int Ax, int Ay,int Bx, int By)
{
    // Ax must be equal to By !!!
    //
    //    if (Ax != By) {
    //        printf("Error : Multiplication by incompatable matrices\n");
    //        exit(1);
    //    }
   
    int i,j,k;
    
    for (i=0;i<Ay;i++){
        for(j=0;j<Bx;j++){
            
            O[i*Bx+j] = 0;
            for(k=0;k<Ax;k++)
                O[i*Bx+j] += A[i*Ax+k] * B[k*Bx+j];
        }
    }
    return ;
}

// 3x3 matrix inversion taken from
// http://www.c4learn.com/c-programs/c-program-to-find-inverse-of-3-x-3.html
// which has a very good graphical description
// double *A - pointer to an array of 9 doubles arranged {{c0r0, c1r0, c2r0},{c0r1,c1r1,c2r1},{c0r2,c1r2,c2r2}} c=column; r=row
// double *O - pointer to output array of 9 doubles arranged {{c0r0, c1r0, c2r0},{c0r1,c1r1,c2r1},{c0r2,c1r2,c2r2}} c=column; r=row
//
void mat3by3inv(double *A, double *O){
    double a[3][6];
    int y,x,i,j;
    for(y=0;y<3;y++){    // Append Unit Matrix
        for(x=0;x<6;x++){
            if(x<3){
                a[y][x] = A[y*3+x] ;
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
            O[i*3+j] = a[i][j+3] ;
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

