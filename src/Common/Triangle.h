/***************************************************************************
 * 
 *           Module :  Triangle.h
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 20/05/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018

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
#import "MVectorObjc.h"

@interface Triangle : NSObject
@property NSString *materialName ;
@property int matId ;
@property int triId ;
- (id) initWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c;
- (id) initWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c andId: (int) val;
- (id) initWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c andNormal:(MVector *)N;
+ (id) TriangleWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c;
+ (id) TriangleWithVerticesAa: (MVector *) a Bb: (MVector *) b Cc: (MVector *) c andId: (int) val;
+ (id) TriangleWithVerticesAa:(MVector *)a Bb:(MVector *)b Cc:(MVector *)c andNormal:(MVector *) N;
- (void) setVertex: (int) index toValue: (MVector *) coordinate;
- (BOOL) isPlanarInDimension: (unsigned int) k ;
- (MVector *) vertexAt: (int) index;
- (NSComparisonResult) compareTriangle: (Triangle *) t ;
- (double) minInDim: (unsigned int) k ;
- (double) maxInDim: (unsigned int) k ;
- (BOOL) isEqual:(id)object ;
- (NSUInteger) hash ;
- (void) setTriId: (int) val ;
- (void) setMaterialName: (NSString *) name ;
- (NSString *) materialName ;
- (MVector *) normal ;
- (MVector *) midpoint ;
- (double) area ;
- (double *) globalToLocalMatrix ;
- (double *) localToGlobalMatrix ;
@end

