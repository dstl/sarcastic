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

void POCalculation(TriCoords tri,
                   Ray Ri_bar,      // Incident Ray
                   SPVector p,      // Global coordinates of scattered point
                   SPVector RxPnt,  // Observation point
                   SPCmplx Ei,      // incident E field at point P
                   double a,        // Ray size x-dimension
                   double b,        // Ray size y-dimension
                   SPCmplx *Es      // Scattered E field
);
void test1();
void test2();
void test3();

int main(int argc, const char * argv[])
{
    test1() ; // Print E field over full spehrical coords
//    test2() ; // print phase over full spherical coords
//    test3() ; // print phase as function of scattering point across plate
  
    return 0;
}

void test1(){
    // Set up a triangle
    //
    SPVector AA,BB,CC ;
    VECT_CREATE(10, 0, 0, AA);
    VECT_CREATE(0, 10, 0, BB);
    VECT_CREATE(-10, -10, 0, CC);
    TriCoords triCoords ;
    triCoords.A=AA ;
    triCoords.B=BB ;
    triCoords.Cc=CC ;
    
    // Create a ray
    //
    Ray r ;
    SPVector hp, dir;
    VECT_CREATE(20.0, 20.0, 20.0, r.org);
    VECT_CREATE(0.0, 0.0, 0.0, hp);
    VECT_SUB(hp, r.org, dir);
    r.len = VECT_MAG(dir);
    VECT_NORM(dir, r.dir);
    r.pow = 1000000.0;
    
    printf("Transmit location %f,%f,%f\n", r.org.x,r.org.y,r.org.z);
    printf("Hipoint location  %f,%f,%f\n", hp.x,hp.y,hp.z);
    printf(" Facet:\n");
    printf("    %f%f%f\n",AA.x,AA.y,AA.z);
    printf("    %f%f%f\n",BB.x,BB.y,BB.z);
    printf("    %f%f%f\n",CC.x,CC.y,CC.z);
    printf("    %f%f%f\n",AA.x,AA.y,AA.z);
    
    
    int iphi, niphis;
    int itheta, nitheta;
    niphis = 360 ;
    nitheta = 100 ;
    
    double deltaiphi, deltaitheta;
    deltaiphi = 2*SIPC_pi / niphis ;
    deltaitheta = SIPC_pi / (2*nitheta) ;
    
    double phi_s, theta_s;
    SPVector RxPnt, op ;
    double obsDist = 10000 ;
    SPCmplx Ei, Es;
    
    for(itheta=0; itheta<nitheta; itheta++){
        theta_s = itheta * deltaitheta ;
        
        //    theta_s = DEG2RAD(45) ;
        //    printf("Observation Incidence Angle : %f deg\n", RAD2DEG(theta_s));
        
        for(iphi=0; iphi < niphis; iphi++){
            phi_s = iphi * deltaiphi ;
            
            //    phi_s = DEG2RAD(270);
            //    printf("Observation  Azimuth Angle : %f deg\n", RAD2DEG(phi_s));
            
            
            RxPnt.x = obsDist * sin(theta_s) * cos(phi_s);
            RxPnt.y = obsDist * sin(theta_s) * sin(phi_s);
            RxPnt.z = obsDist * cos(theta_s);
            CMPLX_F_MAKE(r.pow, 0.0, Ei);
            
            POCalculation(triCoords, r, hp, RxPnt, Ei, 0.1, 0.1, &Es);
            
            op.x = CMPLX_MAG(Es) * sin(theta_s) * cos(phi_s) ;
            op.y = CMPLX_MAG(Es) * sin(theta_s) * sin(phi_s) ;
            op.z = CMPLX_MAG(Es) * cos(theta_s) ;
            
            if(CMPLX_MAG(Es) == 0){
                printf("%f, %f, %f \n",(CMPLX_MAG(Es)),phi_s,theta_s);
                
            }else{
                printf("%f, %f, %f \n",10*log(CMPLX_MAG(Es)),phi_s,theta_s);
            }
            
        }
    }
    return ;
    
}

void test2(){
    // Set up a triangle
    //
    SPVector AA,BB,CC ;
    VECT_CREATE(10, 0, 0, AA);
    VECT_CREATE(0, 10, 0, BB);
    VECT_CREATE(-10, -10, 0, CC);
    TriCoords triCoords ;
    triCoords.A=AA ;
    triCoords.B=BB ;
    triCoords.Cc=CC ;
    
    // Create a ray
    //
    Ray r ;
    SPVector hp, dir;
    VECT_CREATE(20.0, 20.0, 20.0, r.org);
    VECT_CREATE(0.0, 0.0, 0.0, hp);
    VECT_SUB(hp, r.org, dir);
    r.len = VECT_MAG(dir);
    VECT_NORM(dir, r.dir);
    r.pow = 1000000.0;
    
    printf("Transmit location %f,%f,%f\n", r.org.x,r.org.y,r.org.z);
    printf("Hipoint location  %f,%f,%f\n", hp.x,hp.y,hp.z);
    printf(" Facet:\n");
    printf("    %f%f%f\n",AA.x,AA.y,AA.z);
    printf("    %f%f%f\n",BB.x,BB.y,BB.z);
    printf("    %f%f%f\n",CC.x,CC.y,CC.z);
    printf("    %f%f%f\n",AA.x,AA.y,AA.z);
    
    
    int iphi, niphis;
    int itheta, nitheta;
    niphis = 360 ;
    nitheta = 100 ;
    
    double deltaiphi, deltaitheta;
    deltaiphi = 2*SIPC_pi / niphis ;
    deltaitheta = SIPC_pi / (2*nitheta) ;
    
    double phi_s, theta_s;
    SPVector RxPnt, op ;
    double obsDist = 10000 ;
    SPCmplx Ei, Es;
    
    for(itheta=0; itheta<nitheta; itheta++){
        theta_s = itheta * deltaitheta ;
        
        //    theta_s = DEG2RAD(45) ;
        //    printf("Observation Incidence Angle : %f deg\n", RAD2DEG(theta_s));
        
        for(iphi=0; iphi < niphis; iphi++){
            phi_s = iphi * deltaiphi ;
            
            //    phi_s = DEG2RAD(270);
            //    printf("Observation  Azimuth Angle : %f deg\n", RAD2DEG(phi_s));
            
            
            RxPnt.x = obsDist * sin(theta_s) * cos(phi_s);
            RxPnt.y = obsDist * sin(theta_s) * sin(phi_s);
            RxPnt.z = obsDist * cos(theta_s);
            CMPLX_F_MAKE(r.pow, 0.0, Ei);
            
            POCalculation(triCoords, r, hp, RxPnt, Ei, 0.1, 0.1, &Es);
            
            op.x = CMPLX_MAG(Es) * sin(theta_s) * cos(phi_s) ;
            op.y = CMPLX_MAG(Es) * sin(theta_s) * sin(phi_s) ;
            op.z = CMPLX_MAG(Es) * cos(theta_s) ;
            
            if(CMPLX_MAG(Es) == 0){
                printf("%f, %f, %f \n",-10*(CMPLX_PHASE(Es)),phi_s,theta_s);
                
            }else{
                printf("%f, %f, %f \n",-10*(CMPLX_PHASE(Es)),phi_s,theta_s);
            }
            
        }
    }
    return ;
    
}

void test3(){
    
    double xsize = 0.4;
    double ysize = 0.4;
    double lambda = 0.04;
    double xpos,ypos;
    int rn;
    int npoints = 100;
    SPCmplx Ei, Es;
    float phs,real,imag;
    
    // Set up a triangle
    //
    SPVector AA,BB,CC ;
    VECT_CREATE(10, 0, 0, AA);
    VECT_CREATE(0, 10, 0, BB);
    VECT_CREATE(-10, -10, 0, CC);
    TriCoords triCoords ;
    triCoords.A=AA ;
    triCoords.B=BB ;
    triCoords.Cc=CC ;
    
    // Create a ray
    //
    Ray r ;
    SPVector hp, dir, obsPnt;
    VECT_CREATE(0, 0.5, 11.002, r.org);
    VECT_CREATE(0, 0.5, 11.002, obsPnt);
    

    srand(time(NULL));

    for (int n=0; n<npoints; n++){
        rn = rand() % 1000 ;
        xpos = (((double)rn / 1000) -0.5) * xsize ;
        
        rn = rand() % 1000 ;
        ypos = (((double)rn / 1000) -0.5) * ysize ;
        
        VECT_CREATE(xpos, ypos, 0.0, hp);
        VECT_SUB(hp, r.org, dir);
        r.len = VECT_MAG(dir);
        VECT_NORM(dir, r.dir);
        r.pow = 1000000.0;
        phs = (-2*SIPC_pi/lambda) * r.len ;
        real = r.pow * cos(phs) ;
        imag = r.pow * sin(phs) ;
        CMPLX_F_MAKE(real, imag, Ei);
        
        phs = fmod((-2*SIPC_pi*r.len/lambda), 2*SIPC_pi);
        
        POCalculation(triCoords, r, hp, obsPnt, Ei, 1, 1, &Es);
        
        printf("%f, %f, %f\n",xpos,ypos,CMPLX_PHASE(Es));
        phs = fmod((-4*SIPC_pi*r.len/lambda), 2*SIPC_pi);


    }
    
    
    
    
    return ;
    
}


