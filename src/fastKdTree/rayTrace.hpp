/***************************************************************************
 * 
 *           Module :  rayTrace.hpp
 *          Program :  fastKdTree
 *       Created by :  Darren Muff on 05/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  04-Nov-2018
 *      Description :  Program to cast an array of rays through a polygon
 *                     hierarchy. The program uses a stackless traversal 
 *                     through a Kd-Tree [1]
 *
 *        [1] Popov, S., et al. "Stackless kd‐tree traversal for high 
 *            performance GPU ray tracing." Computer Graphics Forum. 
 *            Vol. 26. No. 3. Oxford, UK: Blackwell Publishing Ltd, 2007.
 *
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 * 
 ***************************************************************************/

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
    int      id;     // unique identifier for this ray
} Ray;

static const unsigned int quickmodulo[] = {0,1,2,0,1};

void rayTrace(TriangleMesh *mesh, kdTree::KdData *kdTree, int *numNodesInTree, int nAzRays, int nElRays) ;
void Intersect(ATS *tri, Ray *ray, Hit *hit) ;
void BuildRopesAndBoxes(kdTree::KdData * Node, int *RS, AABB aabb, kdTree::KdData * KdTree) ;
void shootRay(kdTree::KdData * KdTree,ATS * accelTriangles,const int nRays, Ray * rays, Hit *hits);
void stacklessTraverse(const int ind, kdTree::KdData * KdTree, ATS * accelTriangles, const int nRays, Ray * rays, Hit *hits );

#endif /* rayTrace_hpp */
