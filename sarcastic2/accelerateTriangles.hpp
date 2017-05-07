//
//  accelerateTriangles.hpp
//  sarcastic
//
//  Created by Darren Muff on 17/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef accelerateTriangles_hpp
#define accelerateTriangles_hpp

#include <stdio.h>
#include "Trianglemesh.hpp"

typedef struct ATS {
    int  triNum;    // Triangle ID
    double d;         // Constant of plane equation
    double nd_u;      // Normal.u / normal.k
    double nd_v;      // normal.v / normal.k
    int k;          // projection dimension
    double kbu;
    double kbv;
    double kbd;
    double kcu;
    double kcv;
    double kcd;
    int textureInd;
} ATS;

void accelerateTriangles(TriangleMesh *mesh, ATS **accelTriangles) ;

#endif /* accelerateTriangles_hpp */
