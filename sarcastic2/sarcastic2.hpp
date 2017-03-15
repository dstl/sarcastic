//
//  sarcastic2.hpp
//  sarcastic
//
//  Created by Darren Muff on 15/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef sarcastic2_h
#define sarcastic2_h

#include <SIlib2/SIlib2.h>
#include "buildTree.hpp"

typedef struct HitPoint {
    SPVector hit;       // Location of hitpoint in x,y,z
    int tri;            // index of triangle that this hit is on
} HitPoint ;

typedef struct Ray {
    SPVector org;    // Origin
    SPVector dir;    // Direction
    double   pow;    // Power for this ray
    double   len;    // Distance travelled to this ray's origin from transmission
    SPVector pol ;   // unit vector of direction of E field of ray
} Ray;

typedef struct SPCmplxD{
    double r ;
    double i ;
} SPCmplxD;

typedef struct rangeAndPower {
    double  range ;
    SPCmplxD Es ;
} rangeAndPower ;


void rayTrace(kdTree::KdData *tree, int treeSize, TriangleMesh *newMesh, HitPoint **hitPoints, Ray **incidentRays, Ray **observationRays, int *nHits) ;
void POFields(const TriangleMesh *mesh, const HitPoint *hitPoints, const Ray *incidentRays, const Ray *observationRays, const int nHits, rangeAndPower **rnp) ;

#endif /* sarcastic2_h */
