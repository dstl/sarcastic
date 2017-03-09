//
//  clipToAABB.hpp
//  sarcastic
//
//  Created by Darren Muff on 09/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef clipToAABB_hpp
#define clipToAABB_hpp

#include <stdio.h>
#include <SIlib2/SIlib2.h>
#include "AABB.hpp"

typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000
#define NOINTERSECTION -1


int     clipToAABB    (AABB boundingBox, SPVector *lineStart, SPVector *lineEnd);
void    ClipToBox     (SPVector *p0,     SPVector *p1,        SPVector min, SPVector max, int *status);
OutCode ComputeOutCode(SPVector p,       SPVector min,        SPVector max);

#endif /* clipToAABB_hpp */
