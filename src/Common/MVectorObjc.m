/***************************************************************************
 * 
 *           Module :  MVectorObjc.m
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 19/05/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *      Description :  Mathematical vector classes in Objective C
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

#import "MVectorObjc.h"

@implementation MVector {
    double _cell[3];
}
- (id) init {
    self = [super init];
    if(self){
        _cell[0] = _cell[1] = _cell[2] = 0.f;
    }
    return self;
}
- (id) initWithValuesX: (double) x Y: (double) y Z: (double) z {
    self = [super init];
    if (self) {
        _cell[0] = x;
        _cell[1] = y;
        _cell[2] = z;
    }
    return self;
}
- (id) initWithMVector : (MVector *) v {
    self = [super init];
    if (self) {
        _cell[0] = [v X];
        _cell[1] = [v Y];
        _cell[2] = [v Z];
    }
    return self;
}
+ (id) MVectorWithMVector: (MVector *) v {return [[MVector alloc] initWithMVector:v];}
+ (id) MVectorWithValuesX: (double) x Y: (double) y Z: (double) z {return [[MVector alloc] initWithValuesX:x Y:y Z:z]; }
- (void) setX: (double) x { _cell[0] = x; }
- (void) setY: (double) y { _cell[1] = y; }
- (void) setZ: (double) z { _cell[2] = z; }
- (void) setX: (double) x withY: (double) y andZ: (double) z {
    _cell[0] = x; _cell[1] = y; _cell[2] = z;
}
- (void) setCell:(NSUInteger)index with:(double)value {
    _cell[index] = value;
}
- (double) X { return _cell[0]; }
- (double) Y { return _cell[1]; }
- (double) Z { return _cell[2]; }
- (double) cell: (NSUInteger) index { return _cell[index]; }

- (double) length { return (double)sqrt( _cell[0] * _cell[0] + _cell[1] * _cell[1] + _cell[2] * _cell[2] ); }
- (double) mag {return [self length]; }
- (void) normalise { double l = 1.0f / [self length]; _cell[0] *= l; _cell[1] *= l; _cell[2] *= l; }
- (MVector *) norm { double l = 1.0f / [self length]; return [[MVector alloc] initWithValuesX:_cell[0]*l Y:_cell[1]*l Z:_cell[2]*l]; }
- (double) sqrLength { return _cell[0] * _cell[0] + _cell[1] * _cell[1] + _cell[2] * _cell[2]; }
- (double) dot: (MVector *) v { return [v X]*_cell[0] + [v Y]*_cell[1] + [v Z]*_cell[2]; }
- (MVector *) cross : (MVector *) v {return [[MVector alloc] initWithValuesX:_cell[1]*[v Z]-_cell[2]*[v Y] Y:_cell[2]*[v X]-_cell[0]*[v Z] Z:_cell[0]*[v Y]-_cell[1]*[v X]]; }
- (void) negate { _cell[0] *= -1.f; _cell[1] *= -1.f; _cell[2] *= -1.f; }
- (void) minus  { _cell[0] *= -1.f; _cell[1] *= -1.f; _cell[2] *= -1.f; }
- (void) rotateX: (double) a { 
    double c = cos(a); 
    double s = sin(a); 
    MVector *t = [[MVector alloc] initWithValuesX:_cell[0] Y:_cell[1]*c-_cell[2]*s Z:_cell[1]*s+_cell[2]*c ]; 
    _cell[0] = [t X]; _cell[1] = [t Y]; _cell[2]=[t Z];

}
- (void) rotateY: (double) a { 
    double c = cos(a); 
    double s = sin(a);
    MVector *t = [[MVector alloc] initWithValuesX:_cell[0]*c+_cell[2]*s Y:_cell[1] Z:_cell[2]*c-_cell[0]*s ];
    _cell[0] = [t X]; _cell[1]=[t Y]; _cell[2]=[t Z];
}
- (void) rotateZ: (double) a {
    double c = cos(a); 
    double s = sin(a); 
    MVector *t = [[MVector alloc] initWithValuesX:_cell[0]*c-_cell[1]*s Y:_cell[0]*s+_cell[1]*c Z:_cell[2] ];
    _cell[0] = [t X]; _cell[1]=[t Y]; _cell[2]=[t Z];
}
- (void) rotateAroundAxis: (MVector *) axis byAngle: (double) angle {
    if([axis X] == 0 && [axis Y] == 0 ){
        [self rotateZ:angle];
    }else if([axis X] == 0 && [axis Z] == 0){
        [self rotateY:angle];
    }else if([axis Y] == 0 && [axis Z] == 0){
        [self rotateX:angle];
    }else{
        double thetaz = atan([axis Y]/[axis X]);
        double thetay = atan( sqrt( [axis X] * [axis X] + [axis Y] * [axis Y] ) / [axis Z] );
        
        // First rotate around the z axis
        [self rotateZ: -thetaz];
        
        // Now rotate around y axis
        [self rotateY: -thetay];
        
        // now rotate around z by the required angle theta
        [self rotateZ: angle];
        
        // Now add on the rotation axis around y and z to get back to the original reference frame
        [self rotateY:thetay];
        [self rotateZ:thetaz];
    }
    return;
}

- (NSString *) description {
    return [NSString stringWithFormat:@"[%1.5f,%1.5f,%1.5f]",_cell[0],_cell[1],_cell[2]];
}

- (id) copyWithZone:(NSZone *)zone { 
    MVector * v = [[MVector allocWithZone:zone] initWithMVector:self];
    return v; 
}

- (MVector *) add:(MVector *)vector {
    return [[MVector alloc]initWithValuesX:_cell[0]+[vector X] Y:_cell[1]+[vector Y] Z:_cell[2]+[vector Z]];
}

- (MVector *) subtract:(MVector *)vector {
    return [[MVector alloc]initWithValuesX:_cell[0]-[vector X] Y:_cell[1]-[vector Y] Z:_cell[2]-[vector Z]];
}

- (MVector *) multiply: (double) n {
    return [[MVector alloc] initWithValuesX:_cell[0]*n Y:_cell[1]*n Z:_cell[2]*n ];
}

- (MVector *) divide: (double) n {
    return [[MVector alloc] initWithValuesX:_cell[0]/n Y:_cell[1]/n Z:_cell[2]/n];
}
- (BOOL) isEqual:(id)other{
    if (other == self)              // Pointer test
        return YES;
    return ( ([self X]==[other X]) && ([self Y]==[other Y]) && ([self Z]==[other Z]) );
}
- (NSUInteger) hash{
    int prime = 31;
    unsigned long result =1;
    for(int k=0; k<3; k++)
        result = prime * result + (int) ((unsigned long)[self cell:k] ^ ((unsigned long)[self cell:k] >> 32));
    return result;
}

@end
