//
//  main.c
//  sarcastic
//
//  Created by Darren Muff on 03/10/2014.
//  Copyright (c) 2014 Dstl. All rights reserved.
//

#include <stdio.h>
#include <SIlib.h>
#include "POTriangle.h"
#include "matrixMultiplication.h"
#define LAMBDA 0.04
#define TXPOL "H"
#define RXPOL "V"

void buildTriangle(SPVector AA, SPVector BB, SPVector CC, triangle * tri) ;

int main(int argc, const char * argv[])
{
    
    triangle tri1, tri2 ;
    SPVector AA,BB,CC ;
    VECT_CREATE(  0.1, 0.0, 0.0, AA);
    VECT_CREATE(  0.0, 0.1, 0.0, BB);
    VECT_CREATE( -0.1, 0.0, 0.0, CC);
    buildTriangle(AA, BB, CC, &tri1) ;
    
    VECT_CREATE( -0.1,  0.0, 0.0, AA);
    VECT_CREATE(  0.0, -0.1, 0.0, BB);
    VECT_CREATE(  0.1,  0.0, 0.0, CC);
    buildTriangle(AA, BB, CC, &tri2) ;
    
    // Create a ray
    //
    Ray r1, r2 ;
    SPVector hp, dir;
    VECT_CREATE(0.0, 0.0, 200.0, r1.org);
    VECT_CREATE(0.0, 0.0, 0.0, hp);
    VECT_SUB(hp, r1.org, dir);
    r1.len = 0.0;
    r1.pow = 1.0e12 ;
    r1.pow = r1.pow / (SIPC_pi*(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z));
    VECT_NORM(dir, r1.dir);
    SPVector zhat, Hpol, Vpol ;
    VECT_CREATE(0, 0, 1, zhat);
    if(fabs(VECT_DOT(r1.dir, zhat)) == 1.0){
        VECT_CREATE(1, 0, 0, Hpol);
    }else{
        VECT_CROSS(r1.dir, zhat, Hpol) ;
        VECT_NORM(Hpol, Hpol) ;
    }
    VECT_CROSS(Hpol, r1.dir, Vpol) ;
    if(TXPOL == "V"){
        r1.pol = Vpol ;
    }else{
        r1.pol = Hpol ;
    }
    
    VECT_CREATE(0.0, 0.0, 200.0, r2.org);
    VECT_CREATE(0.0, 0.0, 0.0, hp);
    VECT_SUB(hp, r2.org, dir);
    r2.len = 0.0;
    r2.pow = 1.0e12 ;
    r2.pow = r2.pow / (SIPC_pi*(dir.x*dir.x+dir.y*dir.y+dir.z*dir.z));
    VECT_NORM(dir, r2.dir);
    VECT_CREATE(0, 0, 1, zhat);
    if(fabs(VECT_DOT(r2.dir, zhat)) == 1.0){
        VECT_CREATE(1, 0, 0, Hpol);
    }else{
        VECT_CROSS(r2.dir, zhat, Hpol) ;
        VECT_NORM(Hpol, Hpol) ;
    }
    VECT_CROSS(Hpol, r2.dir, Vpol) ;
    r2.pol = Vpol ;

    
    printf("Transmit location %f,%f,%f\n", r1.org.x,r1.org.y,r1.org.z);
    printf("Hipoint location  %f,%f,%f\n", hp.x,hp.y,hp.z);
    printf(" Facet:\n");
    printf("    %8.2f %8.2f %8.2f\n",AA.x,AA.y,AA.z);
    printf("    %8.2f %8.2f %8.2f\n",BB.x,BB.y,BB.z);
    printf("    %8.2f %8.2f %8.2f\n",CC.x,CC.y,CC.z);
    printf("    %8.2f %8.2f %8.2f\n",AA.x,AA.y,AA.z);
    
    int iphi, niphis;
    int itheta, nitheta;
    double startPhi,endPhi,startTheta,endTheta;
    niphis     = 360 ;
    nitheta    = 180 ;
    startPhi   = 0;
    endPhi     = DEG2RAD(360.0) ;
    startTheta = DEG2RAD(0) ;
    endTheta   = SIPC_pi ;
    
    double deltaiphi, deltaitheta;
    deltaiphi = (endPhi-startPhi) / niphis ;
    deltaitheta = (endTheta-startTheta) / (2*nitheta) ;
    double phi_s, theta_s;
    SPVector RxPnt, obsDir ;
    double obsDist = 10000 ;
    
    for(itheta=0; itheta<nitheta; itheta++){
        theta_s = startTheta + itheta * deltaitheta ;
        
        for(iphi=0; iphi < niphis; iphi++){
            phi_s = startPhi + iphi * deltaiphi ;
            
            obsDir.x = sin(theta_s) * cos(phi_s);
            obsDir.y = sin(theta_s) * sin(phi_s);
            obsDir.z = cos(theta_s);
            VECT_SCMULT(obsDir, obsDist, RxPnt) ;
            
            // Define unit vectors for V & H fields
            // We do this as the definition of H and V from the sensor may not be truly horizontal
            // or vertical from the sensor. By allowing us to define the V & H directions we can
            // accurately model the V & H at the sensor.
            // We also do it here as we dont want POTriangle() calculating this each time it is called
            //
            SPVector RXVdir, RXHdir, z_hat ;
            VECT_CREATE(0, 0, 1, z_hat) ;
            if(VECT_DOT(z_hat, obsDir) == 1.0){
                // Ie looking from above - for diagnostic purposes - remove from SAR
                // simulation as this never occurs.
                //
                VECT_CREATE(r1.pol.x, r1.pol.y, 0, RXVdir) ;
                VECT_CROSS(RXVdir, z_hat, RXHdir) ;
            }else{
                VECT_CROSS(z_hat, obsDir, RXHdir) ;
                VECT_CROSS(obsDir, RXHdir, RXVdir) ;
            }
            VECT_NORM(RXVdir, RXVdir) ;
            VECT_NORM(RXHdir, RXHdir) ;
            
            SPCmplx EsV1, EsH1, EsV2, EsH2, EsV, EsH ;
            POTriangle(tri1, r1, RxPnt, LAMBDA, RXVdir, RXHdir, &EsV1, &EsH1) ;
            POTriangle(tri2, r2, RxPnt, LAMBDA, RXVdir, RXHdir, &EsV2, &EsH2) ;
            
            if (RXPOL == "V") {
                CMPLX_ADD(EsV1, EsV2, EsV) ;
                if(CMPLX_MAG(EsV) < 1.0e-4){
                    printf("%f, %f, %f \n",0.0,phi_s,theta_s);
                }else{
                    printf("%f, %f, %f \n",10*log10(CMPLX_MAG(EsV)),phi_s,theta_s);
                }
            }else{
                CMPLX_ADD(EsH1, EsH2, EsH) ;
                if(CMPLX_MAG(EsH) < 1.0e-4){
                    printf("%f, %f, %f \n",0.0,phi_s,theta_s);
                }else{
                    printf("%f, %f, %f \n",10*log10(CMPLX_MAG(EsV)),phi_s,theta_s);
                }

            }
            
            
            
        }
    }

    return 0;
}


void buildTriangle(SPVector AA, SPVector BB, SPVector CC, triangle * tri){
    

    SPVector NN, l1,l2, zhat ;
    double alpha, beta;
    double T_dash[9];
    double T_dashdash[9];
    
    tri->id = 0;
    tri->AA = AA ;
    tri->BB = BB ;
    tri->CC = CC ;
    VECT_SUB(BB, AA, l1) ;
    VECT_SUB(CC, AA, l2) ;
    VECT_CROSS(l1, l2, NN) ;
    tri->area = VECT_MAG(NN) * 0.5 ;
    VECT_NORM(NN, tri->NN) ;
    strcpy (tri->mat, "METAL") ;
    tri->Rs = -66.0;
    for(int imat=0; imat < NMATERIALS; imat++){
        scatProps m = materialProperties[imat] ;
        if( !strcmp(tri->mat, m.matname)){
            tri->Rs = m.resistivity ;
        }
    }
    if(tri->Rs < 0){
        printf("ERROR : Triangle material %s not found\n",tri->mat);
        exit (-1);
    }
    tri->MP.x = (AA.x+BB.x+CC.x) / 3.0 ;
    tri->MP.y = (AA.y+BB.y+CC.y) / 3.0 ;
    tri->MP.z = (AA.z+BB.z+CC.z) / 3.0 ;
    
    VECT_CREATE(0, 0, 1, zhat);
    alpha = atan2(tri->NN.y, tri->NN.x);
    beta  = acos(VECT_DOT(zhat, tri->NN));
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
    
    
    
    matmul(T_dashdash, T_dash, tri->globalToLocalMat, 3, 3, 3, 3);
    mat3by3inv(tri->globalToLocalMat, tri->localToGlobalMat);

    return ;
}