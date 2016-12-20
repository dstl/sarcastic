//
//  AABB.hpp
//  sarcastic
//
//  Created by Darren on 23/10/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#ifndef AABB_hpp
#define AABB_hpp

#include <stdio.h>
#include <SIlib2/SIlib2.h>
#include <vector>

typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000

OutCode ComputeOutCode(SPVector * p, SPVector * min, SPVector * max);

class AABB {
    
public:
    SPVector AA ;   // Lower coord of box
    SPVector BB ;   // Upper coord of box
    
    AABB(){} ;
    AABB(SPVector min, SPVector max): AA(min), BB(max) {}
    AABB(std::vector<SPVector> points) ;
    AABB(std::vector<SPVector> points, int maxDim, float boundMin, float boundMax);
    
    void clipToTriangle(SPVector vertA, SPVector vertB, SPVector vertC, double splitPos, int splitDim, AABB &left, AABB &right);
    float surfaceArea();

};

/*
// Axis Aligned Bounding Box
@interface AABB : NSObject <NSCopying>
@property (copy, nonatomic) MVector * AA;   // Lower coord of box
@property (copy, nonatomic) MVector * BB;   // Upper coord of box
- (id) initWithTriangle: (Triangle *) triangle;
- (id) initWithTriangleArray: (NSArray *) triangles;
- (NSUInteger) maxDim; // largest dimension of AABB
- (void) clipToTriangle: (Triangle *) t ;
- (BOOL) isPlanarInDimension: (NSUInteger) dim ;
- (AABB *) CopyAndSplitLeftDim: (NSUInteger) dim Pos: (float) splitPosition ;
- (AABB *) CopyAndSplitRightDim: (NSUInteger) dim Pos: (float) splitPosition ;
- (float) SurfaceArea ;
- (BOOL) clipLineStart: (MVector *) start End: (MVector *) end ;
@end
 */


#endif /* AABB_hpp */
