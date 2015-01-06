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
#include "matrixMultiplication.h"
#include "materialProperties.h"

@implementation Triangle {
    NSString *_materialName ;
    MVector * _vertices[3]; // Three triangle vertices
    MVector * _normal;      // The triangle normal
    MVector * _midpoint;    // The triangle midpoint
    int       _matId ;       // Id (index) of material this triangle is made of
    int       _triId ;       // Index of this triangle
    double    _area;         // The triangle's area
    double    _g2lMTX[9] ;   // Rotation matrix to convert from global coords into the triangle's local coords
    double    _l2gMTX[9] ;   // Rotation matrix to convert from the triangle's local coords into global coords
    double    _Rs ;          // The triangle's electrical resistivity
 }
@synthesize triId ;

- (id) initWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c{
    self = [super init];
    if (self) {
        double T_dash[9], T_dashdash[9], l, alpha, beta ;
        MVector * x, *zhat ;
        
        _vertices[0]  = [[MVector alloc] initWithMVector:a];
        _vertices[1]  = [[MVector alloc] initWithMVector:b];
        _vertices[2]  = [[MVector alloc] initWithMVector:c];
        _triId        = 0 ;
        _materialName = [NSString stringWithUTF8String: materialProperties[0].matname ] ;
        _matId        = 0 ;
        _Rs           = materialProperties[_matId].resistivity ;
        x             = [[MVector alloc] initWithMVector:[[b subtract:a] cross:[c subtract:a]]] ;
        l             = [x mag] ;
        _area         = l * 0.5 ;
        _normal       = [x multiply:1/l] ;
        _midpoint     = [[MVector alloc] initWithValuesX:([_vertices[0] X]+[_vertices[1] X]+[_vertices[2] X])/3.0
                                                       Y:([_vertices[0] Y]+[_vertices[1] Y]+[_vertices[2] Y])/3.0
                                                       Z:([_vertices[0] Z]+[_vertices[1] Z]+[_vertices[2] Z])/3.0] ;
        zhat          = [[MVector alloc] initWithValuesX:0 Y:0 Z:1] ;
        alpha         = atan2([_normal Y], [_normal X]);
        beta          = acos([zhat dot:_normal]) ;
        T_dash[0]     = cos(alpha);
        T_dash[1]     = sin(alpha);
        T_dash[2]     = 0;
        T_dash[3]     = -sin(alpha);
        T_dash[4]     = cos(alpha);
        T_dash[5]     = 0;
        T_dash[6]     = 0;
        T_dash[7]     = 0;
        T_dash[8]     = 1;
        T_dashdash[0] = cos(beta);
        T_dashdash[1] = 0;
        T_dashdash[2] = -sin(beta);
        T_dashdash[3] = 0;
        T_dashdash[4] = 1;
        T_dashdash[5] = 0;
        T_dashdash[6] = sin(beta);
        T_dashdash[7] = 0;
        T_dashdash[8] = cos(beta);
        
        matmul(T_dashdash, T_dash, _g2lMTX, 3, 3, 3, 3);
        mat3by3inv(_g2lMTX, _l2gMTX);

    }
    return self;
}
- (id) initWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c andId: (int) val{
    self = [super init];
    if (self) {
        double T_dash[9], T_dashdash[9], l, alpha, beta ;
        MVector * x, *zhat ;
        _vertices[0]  = [[MVector alloc] initWithMVector:a];
        _vertices[1]  = [[MVector alloc] initWithMVector:b];
        _vertices[2]  = [[MVector alloc] initWithMVector:c];
        _triId        = val ;
        _materialName = [NSString stringWithUTF8String: materialProperties[0].matname ] ;
        _matId        = 0 ;
        _Rs           = materialProperties[_matId].resistivity ;
        x             = [[MVector alloc] initWithMVector:[[b subtract:a] cross:[c subtract:a]]] ;
        l             = [x mag] ;
        _area         = l * 0.5 ;
        _normal       = [x multiply:1/l] ;
        _midpoint     = [[MVector alloc] initWithValuesX:([_vertices[0] X]+[_vertices[1] X]+[_vertices[2] X])/3.0
                                                       Y:([_vertices[0] Y]+[_vertices[1] Y]+[_vertices[2] Y])/3.0
                                                       Z:([_vertices[0] Z]+[_vertices[1] Z]+[_vertices[2] Z])/3.0] ;
        zhat          = [[MVector alloc] initWithValuesX:0 Y:0 Z:1] ;
        alpha         = atan2([_normal Y], [_normal X]);
        beta          = acos([zhat dot:_normal]) ;
        T_dash[0]     = cos(alpha);
        T_dash[1]     = sin(alpha);
        T_dash[2]     = 0;
        T_dash[3]     = -sin(alpha);
        T_dash[4]     = cos(alpha);
        T_dash[5]     = 0;
        T_dash[6]     = 0;
        T_dash[7]     = 0;
        T_dash[8]     = 1;
        T_dashdash[0] = cos(beta);
        T_dashdash[1] = 0;
        T_dashdash[2] = -sin(beta);
        T_dashdash[3] = 0;
        T_dashdash[4] = 1;
        T_dashdash[5] = 0;
        T_dashdash[6] = sin(beta);
        T_dashdash[7] = 0;
        T_dashdash[8] = cos(beta);
        
        matmul(T_dashdash, T_dash, _g2lMTX, 3, 3, 3, 3);
        mat3by3inv(_g2lMTX, _l2gMTX);

    }
    return self;
}
- (id) initWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c andNormal: (MVector *) N{
    self = [super init];
    if (self) {
        double T_dash[9], T_dashdash[9], l, alpha, beta ;
        MVector * x, *zhat ;
        
        _vertices[0]  = [[MVector alloc] initWithMVector:a];
        _vertices[1]  = [[MVector alloc] initWithMVector:b];
        _vertices[2]  = [[MVector alloc] initWithMVector:c];
        _triId        = 0 ;
        _materialName = [NSString stringWithUTF8String: materialProperties[0].matname ] ;
        _matId        = 0 ;
        _Rs           = materialProperties[_matId].resistivity ;
        x             = [[MVector alloc] initWithMVector:[[b subtract:a] cross:[c subtract:a]]] ;
        l             = [x mag] ;
        _area         = l * 0.5 ;
        _normal       = [[MVector alloc] initWithMVector:N ];
        _midpoint     = [[MVector alloc] initWithValuesX:([_vertices[0] X]+[_vertices[1] X]+[_vertices[2] X])/3.0
                                                       Y:([_vertices[0] Y]+[_vertices[1] Y]+[_vertices[2] Y])/3.0
                                                       Z:([_vertices[0] Z]+[_vertices[1] Z]+[_vertices[2] Z])/3.0] ;
        zhat          = [[MVector alloc] initWithValuesX:0 Y:0 Z:1] ;
        alpha         = atan2([_normal Y], [_normal X]);
        beta          = acos([zhat dot:_normal]) ;
        T_dash[0]     = cos(alpha);
        T_dash[1]     = sin(alpha);
        T_dash[2]     = 0;
        T_dash[3]     = -sin(alpha);
        T_dash[4]     = cos(alpha);
        T_dash[5]     = 0;
        T_dash[6]     = 0;
        T_dash[7]     = 0;
        T_dash[8]     = 1;
        T_dashdash[0] = cos(beta);
        T_dashdash[1] = 0;
        T_dashdash[2] = -sin(beta);
        T_dashdash[3] = 0;
        T_dashdash[4] = 1;
        T_dashdash[5] = 0;
        T_dashdash[6] = sin(beta);
        T_dashdash[7] = 0;
        T_dashdash[8] = cos(beta);
        
        matmul(T_dashdash, T_dash, _g2lMTX, 3, 3, 3, 3);
        mat3by3inv(_g2lMTX, _l2gMTX);
        
    }
    return self;
}

//- (void) setTriId: (int) val { _triId=val; }
- (void) setMaterialName: (NSString *) name {
    _materialName=name;
    _Rs = -66.0;

    for(int imat=0; imat < NMATERIALS; imat++){
        scatProps m = materialProperties[imat] ;
        if( !strcmp([name UTF8String], m.matname)){
            _Rs = m.resistivity ;
        }
    }
    if(_Rs < 0){
        printf("ERROR : Triangle material %s not found\n",[name UTF8String]);
        exit (-1);
    }

}
- (NSString *) materialName { return _materialName; }
- (MVector *)  normal { return _normal ; }
- (MVector *)  midpoint { return _midpoint ; }
- (double)     area { return _area ; }
- (double *)   globalToLocalMatrix { return _g2lMTX; }
- (double *)   localToGlobalMatrix { return _l2gMTX; }

+ (id) TriangleWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c{
    return [[Triangle alloc] initWithVerticesAa:a Bb:b Cc:c];
}

+ (id) TriangleWithVerticesAa:(MVector *)a Bb:(MVector *)b Cc:(MVector *)c andId:(int)val {
    return [[Triangle alloc] initWithVerticesAa:a Bb:b Cc:c andId:val];
}
+ (id) TriangleWithVerticesAa:(MVector *)a Bb:(MVector *)b Cc:(MVector *)c andNormal:(MVector *) N {
    return [[Triangle alloc] initWithVerticesAa:a Bb:b Cc:c andNormal: N];
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
    return [NSString stringWithFormat: @"<%p>[%d]A:%@,B:%@,C:%@ material:%@",self, _triId , [self vertexAt:0],[self vertexAt:1],[self vertexAt:2],[self materialName]];
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




