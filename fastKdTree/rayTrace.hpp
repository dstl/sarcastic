//
//  rayTrace.hpp
//  sarcastic
//
//  Created by Darren Muff on 05/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef rayTrace_hpp
#define rayTrace_hpp

#include <stdio.h>
#include "buildTree.hpp"
#define NOINTERSECTION -1
#define MAXTRAVERSAL 1000           // Maximum kdTree traversal steps before we blow up
#define EPSILON          ((double) 0.0001)

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
typedef struct Hit {
    double dist;
    int trinum;
    double u;
    double v;
} Hit;
typedef struct Ray {
    SPVector org;    // Origin
    SPVector dir;    // Direction
    double   pow;    // Power for this ray
    double   len;    // Distance travelled to this ray's origin from transmission
    SPVector pol ;   // unit vector of direction of E field of ray
} Ray;

void  rayTrace(TriangleMesh *mesh, kdTree::KdData **kdTree, int *numNodesInTree) ;
int clipToAABB(AABB boundingBox, SPVector *lineStart, SPVector *lineEnd) ;
void Intersect(ATS *tri, Ray *ray, Hit *hit) ;

#endif /* rayTrace_hpp */
