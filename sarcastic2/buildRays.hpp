//
//  buildRays.hpp
//  sarcastic
//
//  Created by Darren Muff on 08/05/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef buildRays_hpp
#define buildRays_hpp

#include <stdio.h>
#include "AABB.hpp"
#include "rayTrace.hpp"
void buildRays(Ray **rayArray, int *nRays, int nAzRays, int nElRays, TriangleMesh *mesh, SPVector TxPos,
               double PowPerRay, AABB SceneBoundingBox,SPVector **rayAimPoints);
#endif /* buildRays_hpp */
