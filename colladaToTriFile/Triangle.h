/***************************************************************************
 *  Triangle.h
 *  Sadilac
 *
 *  Created by Muff Darren on 20/05/2012.
 *  Copyright (c) 2014 [Dstl]. All rights reserved.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   20/05/2012
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

@interface Triangle : NSObject
@property int texture ;
@property int triId ;
- (id) initWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c;
- (id) initWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c andId: (int) val;
- (id) initWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c andNormal:(MVector *)normal;
+ (id) TriangleWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c;
+ (id) TriangleWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c andId: (int) val;
+ (id) TriangleWithVerticesA:(MVector *)a B:(MVector *)b C:(MVector *)c andNormal:(MVector *)normal ;
- (void) setVertex: (int) index toValue: (MVector *) coordinate;
- (BOOL) isPlanarInDimension: (unsigned int) k ;
- (MVector *) vertexAt: (int) index;
- (MVector *) normal ;
- (NSComparisonResult) compareTriangle: (Triangle *) t ;
- (double) minInDim: (unsigned int) k ;
- (double) maxInDim: (unsigned int) k ;
- (void) setTexture: (int) texInd ;
- (BOOL) isEqual:(id)object ;
- (NSUInteger) hash ;

@end
