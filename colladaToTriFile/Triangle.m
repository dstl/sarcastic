/***************************************************************************
 *  Triangle.m
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

#import "Triangle.h"
#import "MVectorObjc.h"

@implementation Triangle {
    MVector * _vertices[3];
    MVector * _normal;
}
@synthesize texture=_texture;
@synthesize triId=_triId;

- (id) initWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c{
    self = [super init];
    if (self) {
        _vertices[0] = [[MVector alloc] initWithMVector:a];
        _vertices[1] = [[MVector alloc] initWithMVector:b];
        _vertices[2] = [[MVector alloc] initWithMVector:c];
    }
    return self;
}
- (id) initWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c andId:(int)val{
    self = [super init];
    if (self) {
        _vertices[0] = [[MVector alloc] initWithMVector:a];
        _vertices[1] = [[MVector alloc] initWithMVector:b];
        _vertices[2] = [[MVector alloc] initWithMVector:c];
        [self setTriId:val];
    }
    return self;
}
- (id) initWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c andNormal:(MVector *)normal{
    self = [super init];
    if (self) {
        _vertices[0] = [[MVector alloc] initWithMVector:a];
        _vertices[1] = [[MVector alloc] initWithMVector:b];
        _vertices[2] = [[MVector alloc] initWithMVector:c];
        _normal      = [[MVector alloc] initWithMVector:normal];
    }
    return self;
}

+ (id) TriangleWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c{
    return [[Triangle alloc] initWithVerticesA:a B:b C:c];
}

+ (id) TriangleWithVerticesA:(MVector *)a B:(MVector *)b C:(MVector *)c andId:(int)val {
    return [[Triangle alloc] initWithVerticesA:a B:b C:c andId:val];
}
+ (id) TriangleWithVerticesA:(MVector *)a B:(MVector *)b C:(MVector *)c andNormal:(MVector *)normal {
    return [[Triangle alloc] initWithVerticesA:a B:b C:c andNormal:normal];
}
- (MVector *) vertexAt: (int) index {
    return _vertices[index];
}
- (void) setVertex: (int) index toValue: (MVector *) coordinate{
    _vertices[index] = [coordinate copy];
}
- (MVector *) normal {
    return _normal ;
}
- (BOOL) isPlanarInDimension: (unsigned int) k {
    if ( ([_vertices[0] cell:k] == [_vertices[1] cell:k]) 
        && ([_vertices[1] cell:k] == [_vertices[2] cell:k]) ){
        return YES;
    } else return NO ;
}

- (NSString *) description {
    return [NSString stringWithFormat: @"<%p>[%d]A:%@,B:%@,C:%@ material:%d",self, _triId , [self vertexAt:0],[self vertexAt:1],[self vertexAt:2],[self texture]];
}

- (id) copyWithZone: (NSZone *) zone {
    Triangle * t = [[Triangle allocWithZone:zone] init];
    for (int k=0; k<3; k++)
        [t setVertex:k toValue:[self vertexAt:k]];
    return t;
}
- (NSComparisonResult) compareTriangle: (Triangle *) t {
    unsigned int a = (unsigned int)self;
    unsigned int b = (unsigned int)t;

    if ( a < b )return NSOrderedAscending;
    else if ( a > b) return NSOrderedDescending;
    else return NSOrderedSame;
}
- (double) minInDim: (unsigned int) k {
    double min = 9.0e99;
    for ( int i=0; i<3; i++){
        if( [_vertices[i] cell:k] < min )min=[_vertices[i] cell:k];
    }
    return min;
}
- (double) maxInDim: (unsigned int) k {
    double max = -9.0e99;
    for ( int i=0; i<3; i++){
        if( [_vertices[i] cell:k] > max )max=[_vertices[i] cell:k];
    }
    return max;
}
- (BOOL) isEqual:(id)other{
    if (other == self)      // Pointer test
        return YES;
    return NO;              // Tris are equal only if pointer is the same
}
- (NSUInteger) hash{
    int prime = 31;
    unsigned long result =1;
    for(int k=0; k<3; k++) result = prime * result + [_vertices[k] hash];
    return result;
}
@end




