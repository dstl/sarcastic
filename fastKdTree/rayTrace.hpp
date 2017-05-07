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
#include "clipToAABB.hpp"
extern "C" {
#include "boxMullerRandom.h"
}
#include "accelerateTriangles.hpp"

#define NOINTERSECTION -1
#define MAXTRAVERSAL 1000           // Maximum kdTree traversal steps before we blow up
#define EPSILON          ((double) 0.0001)

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

static const unsigned int quickmodulo[] = {0,1,2,0,1};

void rayTrace(TriangleMesh *mesh, kdTree::KdData *kdTree, int *numNodesInTree, int nAzRays, int nElRays) ;
void Intersect(ATS *tri, Ray *ray, Hit *hit) ;
void BuildRopesAndBoxes(kdTree::KdData * Node, int *RS, AABB aabb, kdTree::KdData * KdTree) ;
void shootRay(kdTree::KdData * KdTree,ATS * accelTriangles,const int nRays, Ray * rays, Hit *hits);
void stacklessTraverse(const int ind, kdTree::KdData * KdTree, ATS * accelTriangles, const int nRays, Ray * rays, Hit *hits );

#endif /* rayTrace_hpp */
