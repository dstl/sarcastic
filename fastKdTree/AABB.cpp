//
//  AABB.cpp
//  sarcastic
//
//  Created by Darren on 23/10/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#include "AABB.hpp"

static const int fastmod[5] = {0,1,2,0,1} ;
SPVector crossingPoint(int dim, double linePos, SPVector start, SPVector end) ;

void AABB::clipToTriangle(SPVector vertA, SPVector vertB, SPVector vertC, double splitPos, int splitDim, AABB &leftBox, AABB &rightBox)
{
    int k = splitDim ;
    std::vector<double> kvals;
    
    kvals.push_back(vertA.cell[k]);
    kvals.push_back(vertB.cell[k]);
    kvals.push_back(vertC.cell[k]);
    
    sort(kvals.begin(),kvals.end());
    
    // make sure min of triangle is left of splitPos and max is right of splitPos else return
    //
    if(kvals[0] > splitPos || kvals[2] < splitPos){
        return;
    }
    
    std::vector<SPVector> points;
    points.push_back(vertA);
    points.push_back(vertB);
    points.push_back(vertC);
    
    // Find when each side of the triangle crosses the splitplane
    //
    SPVector ABxp = crossingPoint(k, splitPos, vertA, vertB) ;
    SPVector BCxp = crossingPoint(k, splitPos, vertB, vertC) ;
    SPVector CAxp = crossingPoint(k, splitPos, vertC, vertA) ;
    
    points.push_back(ABxp);
    points.push_back(BCxp);
    points.push_back(CAxp);
    
    leftBox  = AABB(points, splitDim, -9e9, splitPos);
    rightBox = AABB(points, splitDim, splitPos, 9e9) ;
    
    return;
}

AABB::AABB(std::vector<SPVector> points) {
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> zs;
    for(auto it=points.begin(); it!=points.end(); ++it){
        xs.push_back(it->x);
        ys.push_back(it->y);
        zs.push_back(it->z);
    }
    sort(xs.begin(), xs.end()) ;
    sort(ys.begin(), ys.end()) ;
    sort(zs.begin(), zs.end()) ;
    
    VECT_CREATE(xs[0], ys[0], zs[0], AA);
    VECT_CREATE(*xs.end(), *ys.end(), *zs.end(), BB) ;
}

AABB::AABB(std::vector<SPVector> points, int maxDim, float boundMin, float boundMax){
    
    AABB aabb = AABB(points) ;
    
    AA.cell[maxDim] = (aabb.AA.cell[maxDim] < boundMin) ? boundMin : aabb.AA.cell[maxDim] ;
    BB.cell[maxDim] = (aabb.BB.cell[maxDim] > boundMax) ? boundMax : aabb.BB.cell[maxDim] ;
}

float AABB::surfaceArea(){
    float xlen,ylen,zlen;
    xlen = BB.x - AA.x;
    ylen = BB.y - AA.y;
    zlen = BB.z - AA.z;
    return (xlen*ylen + xlen*zlen + ylen*zlen) * 2 ;
}

SPVector crossingPoint(int dim, double linePos, SPVector start, SPVector end)
// Returns the crossing point of a line through a plane defined in the dim
// dimension at position linepos
//
{
    SPVector ans;
    double sx,sy,ex,ey,grad,y ;
    int k;
    
    ans.cell[dim] = linePos ;
    
    for(int i=1; i<3; ++i){
        k = fastmod[dim+i];
        sx = start.cell[dim] ;
        sy = end.cell[k] ;
        ex = end.cell[dim] ;
        ey = end.cell[k] ;
        if(sx==ex){     // line doesn't change in x
            return end ;
        }
        grad = (ey-sy)/(ex-sx) ;
        y = ((linePos-sx) * grad) + sy ;
        ans.cell[k] = y ;
    }
    
    return ans;
}
