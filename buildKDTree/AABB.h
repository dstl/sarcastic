/***************************************************************************
 *  AABB.h
 *  Sadilac
 *
 *  Created by Muff Darren on 03/06/2012.
 *  Copyright (c) 2012 [dstl]. All rights reserved.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   03/06/2012
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
#import "MVectorObjc.h"
#import "Triangle.h"

typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000

OutCode ComputeOutCode(MVector * p, MVector * min, MVector * max);

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
