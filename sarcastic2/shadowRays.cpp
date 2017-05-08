//
//  shadowRays.cpp
//  sarcastic
//
//  Created by Darren Muff on 08/05/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#include "shadowRays.hpp"
#include <SILib2/SIlib2.h>

void buildShadowRays(int                 nRays,              // The number of reflected rays being considered
                     SPVector            RxPos,              // The Receiver location in x,y,z
                     Ray                 *reflectedRays,     // Array of reflected rays - used for their origin as its the reflection point to Rx
                     Ray                 *shadowRays,        // Output - shadow rays to be tested for occlusion later
                     double              *ranges             // Output - the range to the receiver for each shadowray. Calculated here to save calcs later
)
{
    
    SPVector dir, origin, rxp ;
    double unitvect_tmp, rng ;
    
    for(int ind=0; ind<nRays; ++ind){
        origin              = reflectedRays[ind].org ;
        rxp                 = RxPos ;
        
        VECT_SUB(rxp, origin, dir);
        
        rng                 = VECT_MOD(dir);
        unitvect_tmp        = 1.0 / rng;
        
        VECT_SCMULT(dir, unitvect_tmp, dir) ;
        
        ranges[ind]         = rng ;
        shadowRays[ind].org = origin ;
        shadowRays[ind].dir = dir ;
        shadowRays[ind].pow = reflectedRays[ind].pow ;
        shadowRays[ind].len = reflectedRays[ind].len ;
        shadowRays[ind].pol = reflectedRays[ind].pol ;

    }
    return ;
}
