/***************************************************************************
 * 
 *           Module :  MVectorObjc.h
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 19/05/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *	    Description :  Mathematical vector classes in Objective C
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

#import <Foundation/Foundation.h>

@interface MVector : NSObject <NSCopying>
- (id) init;
- (id) initWithValuesX : (double) x Y: (double) y Z: (double) z ;
- (id) initWithMVector : (MVector *) v;
+ (id) MVectorWithMVector: (MVector *) v ;
+ (id) MVectorWithValuesX: (double) x Y: (double) y Z: (double) z ;
- (void) setX: (double) x;
- (void) setY: (double) y;
- (void) setZ: (double) z;
- (double) X;
- (double) Y;
- (double) Z;
- (double) cell: (NSUInteger) index;
- (void) setX: (double) x withY: (double) y andZ: (double) z ;
- (void) setCell: (NSUInteger) index with: (double) value ;
- (double) length;
- (double) mag;
- (void) normalise;
- (MVector *) norm;
- (double) sqrLength;
- (double) dot: (MVector *) v;
- (MVector *) cross : (MVector *) v;
- (void) negate;
- (void) minus;
- (void) rotateX: (double) angle;
- (void) rotateY: (double) angle;
- (void) rotateZ: (double) angle;
- (void) rotateAroundAxis: (MVector *) axis byAngle: (double) angle ;
- (id) copyWithZone:(NSZone *)zone ;
- (MVector *) add: (MVector *) vector;
- (MVector *) subtract: (MVector *) vector;
- (MVector *) multiply: (double) n;
- (MVector *) divide: (double) n;
- (BOOL) isEqual:(id)object ;
- (NSUInteger) hash ;
@end

