/***************************************************************************
 *  MVector.h
 *  Sadilac
 *
 *  Created by Muff Darren on 19/05/2012.
 *  Copyright (c) 2014 [dstl]. All rights reserved.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   19/05/2012
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
#ifndef REAL
#define REAL double
#endif

@interface MVector : NSObject <NSCopying>
- (id) init;
- (id) initWithValuesX : (REAL) x Y: (REAL) y Z: (REAL) z ;
- (id) initWithMVector : (MVector *) v;
+ (id) MVectorWithMVector: (MVector *) v ;
+ (id) MVectorWithValuesX: (REAL) x Y: (REAL) y Z: (REAL) z ;
- (void) setX: (REAL) x;
- (void) setY: (REAL) y;
- (void) setZ: (REAL) z;
- (REAL) X;
- (REAL) Y;
- (REAL) Z;
- (REAL) cell: (NSUInteger) index;
- (void) setX: (REAL) x withY: (REAL) y andZ: (REAL) z ;
- (void) setCell: (NSUInteger) index with: (REAL) value ;
- (REAL) length;
- (REAL) mag;
- (void) normalise;
- (MVector *) norm;
- (REAL) sqrLength;
- (REAL) dot: (MVector *) v;
- (MVector *) cross : (MVector *) v;
- (void) negate;
- (void) minus;
- (void) rotateX: (REAL) angle;
- (void) rotateY: (REAL) angle;
- (void) rotateZ: (REAL) angle;
- (void) rotateAroundAxis: (MVector *) axis byAngle: (REAL) angle ;
- (id) copyWithZone:(NSZone *)zone ;
- (MVector *) add: (MVector *) vector;
- (MVector *) subtract: (MVector *) vector;
- (MVector *) multiply: (REAL) n;
- (MVector *) divide: (REAL) n;
- (BOOL) isEqual:(id)object ;
- (NSUInteger) hash ;
@end

