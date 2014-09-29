/****************************************************************************
 *  KdTree.h
 *  Sadilac
 *
 *  Created by Darren on 23/05/2012.
 *  Copyright (c) 2012 [dstl]. All rights reserved.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   23/05/2012
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
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
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 ***************************************************************************/

#import <Foundation/Foundation.h>
#import "COLScene.h"
#import "AABB.h"
#import "radarMaterial.h"

#define EVENT_END    0
#define EVENT_PLANAR 1
#define EVENT_START  2
#define LEFT         0
#define RIGHT        1
#define BOTH         2
#define TRAVERSALCOST ((float)(15.0))
#define INTERSECTIONCOST ((float)(20.0))
#define MAXDEPTH 32

void eventCheck(NSArray * events);
float cost(float ProbLeft, float ProbRight, NSUInteger NLeft, NSUInteger NRight, BOOL maxSpace);
NSMutableArray * calculateEventsForTris( NSArray * tris, AABB * V) ;

typedef union {
    struct KdTreeLeaf {
        unsigned int flagDimAndOffset;
        // bits 0..1        : splitting dimension
        // bits 2..30       : offset bits
        // bits 31 (sign)   : flag whether node is a leaf
        float splitPosition;
    } leaf;
    
    struct KdTreeBranch {
        unsigned int flagDimAndOffset;
        // bits 0..30       : offset to first child
        // bit 31 (sign)    : flag whether node is a leaf
        float splitPosition;
    } branch;
} KdData ;


@interface KdTreeEvent : NSObject <NSCopying>
@property unsigned int ptype;     // 0 ending event; 1 planar event; 2 starting event
@property unsigned int dim;       // dimension for this event
@property (assign, nonatomic) Triangle * tri; // Pointer to triangle associated with this event
@property float pos;    // position for this event

- (void) setType: (unsigned int) type Dim: (unsigned int) d Tri: (Triangle *) t Pos: (float) p ;
- (NSComparisonResult) compareEvents: (KdTreeEvent *) e;
- (BOOL) isEqual:(id)object ;
- (NSUInteger) hash ;
@end

@interface KdTreeNode : NSObject <NSCopying> {
    AABB * aabb;
    NSArray * sortDescriptors;
}
@property BOOL isLeaf;
@property float splitPosition;
@property (assign, nonatomic) KdTreeNode * leftNode;
@property (assign, nonatomic) KdTreeNode * rightNode;
@property int axis;
@property int depth; // depth of this node
@property (retain, nonatomic) NSArray * triangles;
@property (retain, nonatomic) NSArray * Events;
@property int cheapSide;
- (void) setAabb: (AABB *) v;
- (AABB *) Aabb ; 
- (id) init ;
- (id) initWithTriangleArray: (NSArray *) triangles andAABB: (AABB *) v ;
- (void) subDivide;
- (void) splitTreeAndEvents;
- (void) findBestSplitPlane;
- (float) SAHDim: (NSUInteger) dim Pos: (float) p Nleft: (NSUInteger) Nl Nright: (NSUInteger) Nr Nplane: (NSUInteger) Np ;
- (KdTreeNode *) splitLeftWithEvents: (NSArray *) LeftEvents ;
- (KdTreeNode *) splitRightWithEvents: (NSArray *) RightEvents ;
- (void) PackNode: (KdTreeNode *) node inArray: (NSMutableArray *) array ;
@end

@interface KdTree : NSObject {
    KdTreeNode * _root;
    triDataStruct * _triangleData;
}
@property (copy, nonatomic) NSMutableArray * packArray;
@property (copy, nonatomic) NSArray * textures;
- (id) initWithScene: (COLScene *) scene ;
- (void) Build ;
- (void) setRoot: (KdTreeNode *) n ;
- (void) PacktoFile: (NSString *) filename ;
- (void) printBoxes ;
- (void) printSummary ;
@end


//
//#define KDTREE_ISLEAF(n)    (n->flagDimAndOffset & (unsigned int)1<<31)
//#define KDTREE_DIMENSION(n) (n->flagDimAndOffset & 0x3)
//#define KDTREE_OFFSET(n)    (n->flagDimAndOffset & 0x7FFFFFFC)
