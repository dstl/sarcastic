//
//  shadowRays.hpp
//  sarcastic
//
//  Created by Darren Muff on 08/05/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef shadowRays_hpp
#define shadowRays_hpp

#include <stdio.h>
#include "rayTrace.hpp"

void buildShadowRays(int                 nRays,              // The number of reflected rays being considered
                SPVector            RxPos,              // The Receiver location in x,y,z
                Ray                 *reflectedRays,     // Array of reflected rays - used for their origin as its the reflection point to Rx
                Ray                 *shadowRays,        // Output - shadow rays to be tested for occlusion later
                double              *ranges             // Output - the range to the receiver for each shadowray. Calculated here to save calcs later
);

#endif /* shadowRays_hpp */
