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
// This function takes three vertices as an input (representing the 3 corners of a triangle) together with a split
// dimension and a split plane position. The Axis-Aligned Bounding Box (AABB) for the triangle is calculated and returned
// in either leftBox or rightBox depending on whether it is before or after the split plane.
// If the triangle straddles the split plane in the indicated dimension then two bounding boxes are calculated with
// each one limited by the split plane and bounding only the parts of teh triangle on the boxes corresponding
// side of the plane.
//
{
    int k = splitDim ;
    std::vector<double> kvals;
    
    // Initialise the retruned bounding boxes so that the calling routine can test to see if they are empty
    //
    SPVector origin = {0.0,0.0,0.0} ;
    leftBox.AA  = origin ;
    leftBox.BB  = origin ;
    rightBox.AA = origin ;
    rightBox.BB = origin ;
      
    kvals.push_back(vertA.cell[k]);
    kvals.push_back(vertB.cell[k]);
    kvals.push_back(vertC.cell[k]);
    
    sort(kvals.begin(),kvals.end());
    
    // make sure min of triangle is left of splitPos and max is right of splitPos else return
    //
    if(kvals[0] > splitPos || kvals[2] < splitPos){
        return;
    }
    
    std::vector<SPVector> leftBoxPoints;
    std::vector<SPVector> rghtBoxPoints;
    
    SPVector * verts = new SPVector [3] ;
    verts[0] = vertA ;
    verts[1] = vertB ;
    verts[2] = vertC ;
    SPVector va,vb;
    
    // loop through each side of the triangle
    //
    for(int i=0; i<3; ++i){
        va = verts[i] ;
        vb = verts[fastmod[i+1]] ;
    
        // Find when each side of the triangle crosses the splitplane
        //
        SPVector ABxp = crossingPoint(k, splitPos, va, vb) ;
        if ( (ABxp.x != vb.x) && (ABxp.y != vb.y) && (ABxp.z != vb.z) ) // ie line AB does pass through plane
        {
            leftBoxPoints.push_back(ABxp)  ;
            rghtBoxPoints.push_back(ABxp)  ;
            if(va.cell[k] <= splitPos) {
                leftBoxPoints.push_back(va) ;
                rghtBoxPoints.push_back(vb) ;
            }else{
                rghtBoxPoints.push_back(va) ;
                leftBoxPoints.push_back(vb) ;
            }
        }else // line does not cross the split plane
        {
            if(va.cell[k] <= splitPos) {
                leftBoxPoints.push_back(va) ;
                leftBoxPoints.push_back(vb) ;
            }
            if(va.cell[k] >= splitPos) {
                rghtBoxPoints.push_back(va) ;
                rghtBoxPoints.push_back(vb) ;
            }
        }
    }
    
    if(leftBoxPoints.size() > 1){
    	leftBox  = AABB(leftBoxPoints);
	}
    if(rghtBoxPoints.size() > 1){
    	rightBox = AABB(rghtBoxPoints);
    }
    delete verts ;
    
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

bool AABB::isEmpty(){
    if (surfaceArea() == 0) {
        return true ;
    }else{
        return false ;
    }
}
SPVector crossingPoint(int dim, double linePos, SPVector start, SPVector end)
// Returns the crossing point of a line through a plane defined in the dim
// dimension at position linepos.
// If line doesn't cross the plane then the end point of the line is returned
//
{
    SPVector ans;
    double s0,s1,s2,e0,e1,e2,grad1,grad2,ans1,ans2 ;
    int dim1,dim2;
    
    dim1 = fastmod[dim+1];
    dim2 = fastmod[dim+2];
    s0 = start.cell[dim] ;
    s1 = start.cell[dim1] ;
    s2 = start.cell[dim2] ;
    e0 = end.cell[dim] ;
    e1 = end.cell[dim1] ;
    e2 = end.cell[dim2] ;

    // Test to see if line is only on one side of teh split plain or does not vary in x
    //
    if( (linePos > e0 && linePos > s0) || (linePos < e0 && linePos < s0 ) || (s0 == e0) ){
        ans = end ;
    }else{
    
        grad1 = (e1-s1)/(e0-s0) ;
        grad2 = (e2-s2)/(e0-s0) ;
        ans1  = ((linePos-s0) * grad1) + s1 ;
        ans2  = ((linePos-s0) * grad2) + s2 ;
        
        ans.cell[dim]  = linePos ;
        ans.cell[dim1] = ans1 ;
        ans.cell[dim2] = ans2 ;
   
    }

    return ans;
}
