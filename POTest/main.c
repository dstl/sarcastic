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
                   SPCmplx *Es      // Scattered E field

);

int main(int argc, const char * argv[])
{
    
    // Set up a triangle
    //
    SPVector AA,BB,CC ;
    VECT_CREATE(10, 10, 0, AA);
    VECT_CREATE(5, 10, 0, BB);
    VECT_CREATE(10, 5, 0, CC);
    TriCoords triCoords ;
    triCoords.A=AA ;
    triCoords.B=BB ;
    triCoords.Cc=CC ;

    // Create a ray
    //
    Ray r ;
    SPVector hp, dir;
    VECT_CREATE(0.0, 0.0, 20.0, r.org);
    VECT_CREATE(7.0, 7.0, 0.0, hp);
    VECT_SUB(hp, r.org, dir);
    r.len = VECT_MAG(dir);
    VECT_NORM(dir, r.dir);
    r.pow = 1.0;
    
    int iphi, niphis;
    int itheta, nitheta;
    niphis = 360 ;
    nitheta = 180 ;
    
    float deltaiphi, deltaitheta;
    deltaiphi = 2*SIPC_pi / niphis ;
    deltaitheta = SIPC_pi / nitheta ;
    
    float phi_s, theta_s;
    SPVector RxPnt, op ;
    float obsDist = 1000 ;
    SPCmplx Ei, Es;
    
    for(iphi=0; iphi < niphis; iphi++){
        phi_s = iphi * deltaiphi ;
        
        for(itheta=0; itheta< nitheta; itheta++){
            theta_s = nitheta * deltaitheta ;
            
            RxPnt.x = obsDist * sinf(theta_s) * cosf(phi_s);
            RxPnt.y = obsDist * sinf(theta_s) * sinf(phi_s);
            RxPnt.z = obsDist * cosf(theta_s);
            CMPLX_F_MAKE(r.pow, 0.0, Ei);
            
            POCalculation(triCoords, r, hp, RxPnt, Ei, &Es);
            
            op.x = CMPLX_MAG(Es) * sinf(theta_s) * cosf(phi_s) ;
            op.y = CMPLX_MAG(Es) * sinf(theta_s) * sinf(phi_s) ;
            op.x = CMPLX_MAG(Es) * cosf(theta_s) ;
            
            printf("%f, %f, %f \n",op.x,op.y,op.z);
            
        }
    }
    
    return 0;
}


