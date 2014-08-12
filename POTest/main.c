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

int main(int argc, const char * argv[])
{
    
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
    VECT_CREATE(0.0, 20.0, 20.0, r.org);
    VECT_CREATE(0.0, 0.0, 0.0, hp);
    VECT_SUB(hp, r.org, dir);
    r.len = VECT_MAG(dir);
    VECT_NORM(dir, r.dir);
    r.pow = 1000.0;
    
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
    double obsDist = 1000 ;
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
    
    POCalculation(triCoords, r, hp, RxPnt, Ei, 0.8, 0.8, &Es);
    
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

    return 0;
}


