//
//  reflect.hpp
//  sarcastic
//
//  Created by Darren Muff on 08/05/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef reflect_hpp
#define reflect_hpp

#include <stdio.h>
#include "rayTrace.hpp"

void reflect(int    nRays,              // Number of rays to reflect
             Ray    *rays,              // Array of rays to consider
             Hit    *hits,              // Array of hit points for each ray
             ATS    *accelTriangles,    // Array of triangles in 'accelerated' format
             Ray    *reflectedRays      // output array of reflected rays
);

#endif /* reflect_hpp */
