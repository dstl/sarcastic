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

@implementation Triangle {
    MVector * _vertices[3];
//    MVector * _normal;
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

+ (id) TriangleWithVerticesA: (MVector *) a B: (MVector *) b C: (MVector *) c{
    return [[Triangle alloc] initWithVerticesA:a B:b C:c];
}

- (id) TriangleWithVerticesA:(MVector *)a B:(MVector *)b C:(MVector *)c andId:(int)val {
    return [[Triangle alloc] initWithVerticesA:a B:b C:c andId:val];
}
- (MVector *) vertexAt: (int) index {
    return _vertices[index];
}
- (void) setVertex: (int) index toValue: (MVector *) coordinate{
    _vertices[index] = [coordinate copy];
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
- (REAL) minInDim: (unsigned int) k {
    REAL min = 9.0e99;
    for ( int i=0; i<3; i++){
        if( [_vertices[i] cell:k] < min )min=[_vertices[i] cell:k];
    }
    return min;
}
- (REAL) maxInDim: (unsigned int) k {
    REAL max = -9.0e99;
    for ( int i=0; i<3; i++){
        if( [_vertices[i] cell:k] > max )max=[_vertices[i] cell:k];
    }
    return max;
}
- (triDataStruct) returnAsData {
    TriAccel * t = [[TriAccel alloc] initWithTriangle:self];
    triDataStruct data;
    data.d = [t d];
    data.nd_u = [t nd_u];
    data.nd_v = [t nd_v];
    data.k = [t k];
    data.kbu = [t kbu];
    data.kbv = [t kbv];
    data.kbd = [t kbd];
    data.kcu = [t kcu];
    data.kcv = [t kcv];
    data.kcd = [t kcd];
    data.texture = [t texture];
    return data;
}
- (BOOL) isEqual:(id)other{
    if (other == self)      // Pointer test
        return YES;
    return NO;              // Tris are equal only if pointer is the same
//    return ([[self vertexAt:0] isEqual:[other vertexAt:0]]
//            && [[self vertexAt:1] isEqual:[other vertexAt:1]]
//            && [[self vertexAt:2] isEqual:[other vertexAt:2]]);
}
- (NSUInteger) hash{
    int prime = 31;
    unsigned long result =1;
    for(int k=0; k<3; k++) result = prime * result + [_vertices[k] hash];
    return result;
}
@end

@implementation TriAccel
@synthesize d=_d;
@synthesize nd_u=_nd_u;
@synthesize nd_v=_nd_v;
@synthesize k=_k;
@synthesize kbu=_kbu;
@synthesize kbv=_kbv;
@synthesize kbd=_kbd;
@synthesize kcu=_kcu;
@synthesize kcv=_kcv;
@synthesize kcd=_kcd;
@synthesize texture=_texture;
- (id) initWithTriangle: (Triangle *) triangle {
    
    self = [super init];
    if (self){
        
        MVector * A = [triangle vertexAt:0];
        MVector * B = [triangle vertexAt:1];
        MVector * C = [triangle vertexAt:2];
        
        // calculate the projection dimension
//        MVector * N = [triangle normal];
        MVector * N = [[B subtract:A] cross:[C subtract:A]];
        REAL nx = [N X];
        REAL ny = [N Y];
        REAL nz = [N Z];
        REAL anx = ABS(nx);
        REAL any = ABS(ny);
        REAL anz = ABS(nz);
        if( anx > any )
            if (anx > anz) _k = 0; /* X */ else _k=2; /* Z */
        else
            if ( any > anz) _k=1; /* Y */ else _k=2; /* Z */
        int u = quickmodulo[_k+1];
        int v = quickmodulo[_k+2];
        
        _nd_u = [N cell:u] / [N cell:_k];
        _nd_v = [N cell:v] / [N cell:_k];
        _d = [A dot:[N divide:[N cell:_k]]];
        
        MVector * b = [C subtract:A];
        MVector * c = [B subtract:A];
        
        REAL denom = 1./(([b cell:u]*[c cell:v]) - ([b cell:v]*[c cell:u]));
                
        _kbu = -[b cell:v] * denom;
        _kbv =  [b cell:u] * denom;
        _kbd =  (([b cell:v] * [A cell:u]) - ([b cell:u] * [A cell:v])) * denom;
        _kcu =  [c cell:v] * denom;
        _kcv = -[c cell:u] * denom;
        _kcd =  (([c cell:u] * [A cell:v]) - ([c cell:v] * [A cell:u])) * denom;
        
        _texture = [triangle texture];
        
    }
    return self;
}

@end






