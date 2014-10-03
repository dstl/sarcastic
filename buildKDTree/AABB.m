/***************************************************************************
 *  AABB.m
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


#import "AABB.h"
#import "PolyClip.h"

int isClockwise(Poly *triangle) ;

@implementation AABB
@synthesize AA;
@synthesize BB;

- (id) initWithTriangle: (Triangle *) triangle{
    self = [super init];
    if(self){
        float minx, miny, minz, maxx, maxy, maxz;
        float f;
        minx = miny = minz = 9e99;
        maxx = maxy = maxz = -9e99;
        for (int i=0; i<3; i++){
            minx = ( (f=[[triangle vertexAt:i] X]) < minx) ? f : minx;
            miny = ( (f=[[triangle vertexAt:i] Y]) < miny) ? f : miny;
            minz = ( (f=[[triangle vertexAt:i] Z]) < minz) ? f : minz;
            
            maxx = ( (f=[[triangle vertexAt:i] X]) > maxx) ? f : maxx;
            maxy = ( (f=[[triangle vertexAt:i] Y]) > maxy) ? f : maxy;
            maxz = ( (f=[[triangle vertexAt:i] Z]) > maxz) ? f : maxz;
        }
        AA = [[MVector alloc] initWithValuesX:minx Y:miny Z:minz];
        BB = [[MVector alloc] initWithValuesX:maxx Y:maxy Z:maxz];
    }
    return self;
}

- (id) initWithTriangleArray: (NSArray *) triangles{
    self = [super init];
    if (self) {
        NSUInteger n = [triangles count];
        float minx, miny, minz, maxx, maxy, maxz;
        float f;
        minx = miny = minz = 9e99;
        maxx = maxy = maxz = -9e99;
        Triangle * tri;
        for (int i=0; i<n; i++){
            tri = [triangles objectAtIndex:i];
            for (int i=0; i<3; i++){
                minx = ( (f=[[tri vertexAt:i] X]) < minx) ? f : minx;
                miny = ( (f=[[tri vertexAt:i] Y]) < miny) ? f : miny;
                minz = ( (f=[[tri vertexAt:i] Z]) < minz) ? f : minz;
                
                maxx = ( (f=[[tri vertexAt:i] X]) > maxx) ? f : maxx;
                maxy = ( (f=[[tri vertexAt:i] Y]) > maxy) ? f : maxy;
                maxz = ( (f=[[tri vertexAt:i] Z]) > maxz) ? f : maxz;
            }
        }
        AA = [[MVector alloc] initWithValuesX:minx Y:miny Z:minz];
        BB = [[MVector alloc] initWithValuesX:maxx Y:maxy Z:maxz];
    }
    return self;
}

- (NSUInteger) maxDim {
    float xs = fabs([BB X] - [AA X]);
    float ys = fabs([BB Y] - [AA Y]);
    float zs = fabs([BB Z] - [AA Z]);
    
    if (xs > ys) {
        if( xs > zs ) return 0;     // x dim
        else return 2;              // z dim
    }else if( ys > zs) return 1;    // y dim
    else return 2;                  // z dim
}

// for polygon intersection to work, the triangle has to be 'clockwise' (mathematically positive).
// The rectangle is positive by construction
// This function takes the triangle and puts it into the correct sense
//
int isClockwise(Poly *triangle){
    Coordinate A,B,C;
    A=triangle->coordinates[0] ;
    B=triangle->coordinates[1] ;
    C=triangle->coordinates[2] ;
    
    float acx = C.x-A.x;
    float acy = C.y-A.y;
    float cbx = B.x-C.x;
    float cby = B.y-C.y;
    
    if ( acx*cby > acy*cbx ){
        // Triangle is clockwise
        return 1 ;
        
    }else{
        // Triangle is counterclockwise
        return 0 ;
    }
}

- (void) clipToTriangle:(Triangle *)t {
    Poly xyRectangle, triangle, xzRectangle ;
    Coordinate tA,tB,tC,rA,rB,rC,rD ;
    int i;
    NSMutableArray * xPoints = [NSMutableArray array];
    NSMutableArray * yPoints = [NSMutableArray array];
    NSMutableArray * zPoints = [NSMutableArray array];
    
    tA.x = [[t vertexAt:0] X];
    tA.y = [[t vertexAt:0] Y];
    tB.x = [[t vertexAt:1] X];
    tB.y = [[t vertexAt:1] Y];
    tC.x = [[t vertexAt:2] X];
    tC.y = [[t vertexAt:2] Y];
    
//    NSLog(@"Triangle is ");
//    printf("%f,%f,%f\n",[[t vertexAt:0] X], [[t vertexAt:0] Y], [[t vertexAt:0] Z]);
//    printf("%f,%f,%f\n",[[t vertexAt:1] X], [[t vertexAt:1] Y], [[t vertexAt:1] Z]);
//    printf("%f,%f,%f\n",[[t vertexAt:2] X], [[t vertexAt:2] Y], [[t vertexAt:2] Z]);
//    printf("%f,%f,%f\n",[[t vertexAt:0] X], [[t vertexAt:0] Y], [[t vertexAt:0] Z]);
//
//    NSLog(@"AABB is ");
//    printf("%f,%f,%f\n", [AA X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [BB X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [BB X], [BB Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [BB Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [BB Z]);
//    printf("%f,%f,%f\n", [BB X], [AA Y], [BB Z]);
//    printf("%f,%f,%f\n", [BB X], [BB Y], [BB Z]);
//    printf("%f,%f,%f\n", [AA X], [BB Y], [BB Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [BB Z]);
    
    triangle.nVerts  = 0;
    addToPoly(&triangle, tA);
    addToPoly(&triangle, tB);
    addToPoly(&triangle, tC);
    if( ! isClockwise(&triangle) ){
        Coordinate tmp = triangle.coordinates[1] ;
        triangle.coordinates[1] = triangle.coordinates[2] ;
        triangle.coordinates[2] = tmp ;
    }
    
    rA.x = [AA X] ;
    rA.y = [AA Y] ;
    rB.x = [AA X] ;
    rB.y = [BB Y] ;
    rC.x = [BB X] ;
    rC.y = [BB Y] ;
    rD.x = [BB X] ;
    rD.y = [AA Y] ;

    xyRectangle.nVerts = 0;
    addToPoly(&xyRectangle, rA) ;
    addToPoly(&xyRectangle, rB) ;
    addToPoly(&xyRectangle, rC) ;
    addToPoly(&xyRectangle, rD) ;
    
    Poly xyIntersect;
    xyIntersect.nVerts=0;
    // if AABB is planar then by definition it cannot contain any triangles
    // however it may not be planar in another dimension so to make sure we
    // correctly calculate AABB intersect Triangle then use the triangle to cookie-cut
    // the rectangle.
    if ([self isPlanarInDimension:0] || [self isPlanarInDimension:1] ){
        xyIntersect = intersection(triangle, xyRectangle );
    }else{
        xyIntersect = intersection(xyRectangle, triangle );
    }
    
    for(i=0; i<xyIntersect.nVerts; i++){
        [xPoints addObject:[NSNumber numberWithFloat: xyIntersect.coordinates[i].x ]];
        [yPoints addObject:[NSNumber numberWithFloat: xyIntersect.coordinates[i].y ]];
    }
    
    tA.x = [[t vertexAt:0] X];
    tA.y = [[t vertexAt:0] Z];
    tB.x = [[t vertexAt:1] X];
    tB.y = [[t vertexAt:1] Z];
    tC.x = [[t vertexAt:2] X];
    tC.y = [[t vertexAt:2] Z];
    
    triangle.nVerts=0;
    addToPoly(&triangle, tA);
    addToPoly(&triangle, tB);
    addToPoly(&triangle, tC);
    if( ! isClockwise(&triangle) ){
        Coordinate tmp = triangle.coordinates[1] ;
        triangle.coordinates[1] = triangle.coordinates[2] ;
        triangle.coordinates[2] = tmp ;
    }
    
    rA.x = [AA X] ;
    rA.y = [AA Z] ;
    rB.x = [AA X] ;
    rB.y = [BB Z] ;
    rC.x = [BB X] ;
    rC.y = [BB Z] ;
    rD.x = [BB X] ;
    rD.y = [AA Z] ;
   
    xzRectangle.nVerts = 0 ;
    addToPoly(&xzRectangle, rA) ;
    addToPoly(&xzRectangle, rB) ;
    addToPoly(&xzRectangle, rC) ;
    addToPoly(&xzRectangle, rD) ;
    Poly xzIntersect;
    xzIntersect.nVerts = 0;
    // if AABB is planar then by definition it cannot contain any triangles
    // however it may not be planar in another dimension so to make sure we
    // correctly calculate AABB intersect Triangle then use the triangle to cookie-cut
    // the rectangle.
    if ([self isPlanarInDimension:2]){
        xzIntersect = intersection(triangle, xzRectangle );
    }else{
        xzIntersect = intersection(xzRectangle, triangle );
    }
    
    for(i=0; i<xzIntersect.nVerts; i++){
        [zPoints addObject:[NSNumber numberWithFloat: xzIntersect.coordinates[i].y ]];
    }
    
    if( [xPoints count] > 0) {
        [xPoints sortUsingSelector:@selector(compare:)] ;
        [AA setX:[[xPoints objectAtIndex:0] doubleValue]] ;
        [BB setX: [[xPoints lastObject] doubleValue]] ;
    }
    
    if ([yPoints count] > 0) {
        [yPoints sortUsingSelector:@selector(compare:)] ;
        [AA setY:[[yPoints objectAtIndex:0] doubleValue]] ;
        [BB setY: [[yPoints lastObject] doubleValue]] ;
    }
    
    if ([zPoints count] > 0) {
        [zPoints sortUsingSelector:@selector(compare:)] ;
        [AA setZ:[[zPoints objectAtIndex:0] doubleValue]] ;
        [BB setZ: [[zPoints lastObject] doubleValue]] ;
    }
    
//    NSLog(@"AABB after clip : ");
//    printf("%f,%f,%f\n", [AA X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [BB X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [BB X], [BB Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [BB Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [BB Z]);
//    printf("%f,%f,%f\n", [BB X], [AA Y], [BB Z]);
//    printf("%f,%f,%f\n", [BB X], [BB Y], [BB Z]);
//    printf("%f,%f,%f\n", [AA X], [BB Y], [BB Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [BB Z]);
    
    return ;

}

/*

 // code commented out for revision completeness 
 // As described below there is a bug in this section that is fixed above.
 //

- (void) clipToTriangle:(Triangle *)t {
    int fastmod[] = {0, 1, 2, 0, 1};
    bool status;
    
    NSLog(@"Triangle is %@", t);
    NSLog(@"AABB is ");
    printf("%f,%f,%f\n", [AA X], [AA Y], [AA Z]);
    printf("%f,%f,%f\n", [BB X], [AA Y], [AA Z]);
    printf("%f,%f,%f\n", [BB X], [BB Y], [AA Z]);
    printf("%f,%f,%f\n", [AA X], [BB Y], [AA Z]);
    printf("%f,%f,%f\n", [AA X], [AA Y], [AA Z]);
    printf("%f,%f,%f\n", [AA X], [AA Y], [BB Z]);
    printf("%f,%f,%f\n", [BB X], [AA Y], [BB Z]);
    printf("%f,%f,%f\n", [BB X], [BB Y], [BB Z]);
    printf("%f,%f,%f\n", [AA X], [BB Y], [BB Z]);
    printf("%f,%f,%f\n", [AA X], [AA Y], [BB Z]);
   

    NSMutableArray * xPoints = [NSMutableArray array];
    NSMutableArray * yPoints = [NSMutableArray array];
    NSMutableArray * zPoints = [NSMutableArray array];
    MVector * start, * end;
    
    // Bug ??
    // consider line of triangle beyond volume V in Y direction. [ClipLineStart End] returns status=0
    // as the line does not intersect volume. the points from the line therefore to not get added to yPoints
    // list. If another side ofthe triangle does clip the volume then its Y-coordinate will get added
    // to the list. When the volume is clipped (below) using the largest y value in the list it will get clipped
    // to the lower value of Y. It should get clipped to the boundary of the volume in y.
    
    for (int v=0; v<3; v++){
        start = [[t vertexAt:v] copy];
        end   = [[t vertexAt:fastmod[v+1]] copy];
        status = [self clipLineStart:start End:end];
        if (status){
            [xPoints addObject:[NSNumber numberWithFloat:[start X]]];
            [xPoints addObject:[NSNumber numberWithFloat:[end   X]]];
            [yPoints addObject:[NSNumber numberWithFloat:[start Y]]];
            [yPoints addObject:[NSNumber numberWithFloat:[end   Y]]];
            [zPoints addObject:[NSNumber numberWithFloat:[start Z]]];
            [zPoints addObject:[NSNumber numberWithFloat:[end   Z]]];
        }
    }

    if ([xPoints count] > 0 ) {
        // reduce AABB to tightly enclose triangle
        [xPoints sortUsingSelector:@selector(compare:)];
        [yPoints sortUsingSelector:@selector(compare:)];
        [zPoints sortUsingSelector:@selector(compare:)];
        [AA setX:[[xPoints objectAtIndex:0] doubleValue]
           withY:[[yPoints objectAtIndex:0] doubleValue]
            andZ:[[zPoints objectAtIndex:0] doubleValue] ];
        [BB setX:[[xPoints lastObject] doubleValue]
           withY:[[yPoints lastObject] doubleValue]
            andZ:[[zPoints lastObject] doubleValue] ];

    }    
    
//    NSLog(@"AABB after clip : ");
//    printf("%f,%f,%f\n", [AA X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [BB X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [BB X], [BB Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [BB Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [AA Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [BB Z]);
//    printf("%f,%f,%f\n", [BB X], [AA Y], [BB Z]);
//    printf("%f,%f,%f\n", [BB X], [BB Y], [BB Z]);
//    printf("%f,%f,%f\n", [AA X], [BB Y], [BB Z]);
//    printf("%f,%f,%f\n", [AA X], [AA Y], [BB Z]);
    
}
*/

- (BOOL) isPlanarInDimension: (NSUInteger) dim {
    if ( [BB cell:dim] == [AA cell:dim] ){
        return YES;
    }
    return  NO;
}
- (AABB *) CopyAndSplitLeftDim: (NSUInteger) dim Pos: (float) splitPosition {
    AABB * sub = [self copy];
    [[sub BB] setCell:dim with:splitPosition];
    return sub;
}
- (AABB *) CopyAndSplitRightDim: (NSUInteger) dim Pos: (float) splitPosition {
    AABB * sub = [self copy];
    [[sub AA] setCell:dim with:splitPosition];
    return sub;
}
- (float) SurfaceArea {
    float xlen,ylen,zlen;
    xlen = [BB X] - [AA X];
    ylen = [BB Y] - [AA Y];
    zlen = [BB Z] - [AA Z];
    return (xlen*ylen + xlen*zlen + ylen*zlen) * 2 ;
}
- (id) copyWithZone:(NSZone *)zone {
    AABB * v = [[AABB allocWithZone:zone] init ];
    [v setAA:AA];
    [v setBB:BB];
    return v;
}
- (NSString *) description {
    NSMutableString * s = [[NSMutableString alloc] init];
    [s appendFormat:@"[%1.5f,%1.5f,%1.5f\\%1.5f,%1.5f,%1.5f]"
     ,[AA X],[AA Y],[AA Z],[BB X],[BB Y],[BB Z]];
    return s ;
}
- (BOOL) clipLineStart:(MVector *)start End:(MVector *)end {
    
    // Cohen–Sutherland clipping algorithm clips a line from
    // P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with 
    // diagonal from (xmin, ymin) to (xmax, ymax).
    
    MVector * min = AA;
    MVector * max = BB;
    
    // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
    OutCode outcode0 = ComputeOutCode(start, min, max);
    OutCode outcode1 = ComputeOutCode(end, min, max);
    BOOL accept = false;
    
    while (true) {
        if (!(outcode0 | outcode1)) { // Bitwise OR is 0. Both ends must be inside box
            accept = true;
            break;
        } else if (outcode0 & outcode1) { // Bitwise AND is not 0. Reject as both ends must be in same region outside box
            break;
        } else {
            // failed both tests, so calculate the line segment to clip
            // from an outside point to an intersection with clip edge
            double x, y, z, m;
            
            // At least one endpoint is outside the clip rectangle; pick it.
            OutCode outcodeOut = outcode0 ? outcode0 : outcode1;
            
            // Now find the intersection point;
            
            if (outcodeOut & FARAABB) {             // Point is beyond the box (y-axis)
                m = ([max Y] - [start Y]) / ([end Y] - [start Y]);
                x = [start X] + m * ([end X] - [start X]) ;
                y = [max Y];
                z = [start Z] + m * ([end Z] - [start Z]) ;
            }else if (outcodeOut & NEARAABB){     // Point is before the box (y-axis)
                m = ([min Y] - [start Y]) / ([end Y] - [start Y]);
                x = [start X] + m * ([end X] - [start X]) ;
                y = [min Y];
                z = [start Z] + m * ([end Z] - [start Z]) ;
            }else if (outcodeOut & LEFTAABB){       // Point is left of box (x-axis)
                m = ([min X] - [start X]) / ([end X] - [start X]);
                x = [min X];
                y = [start Y] + m * ([end Y] - [start Y]);
                z = [start Z] + m * ([end Z] - [start Z]);
            }else if (outcodeOut & RIGHTAABB){      // Point is right of box (x-axis)
                m = ([max X] - [start X]) / ([end X] - [start X]);
                x = [max X];
                y = [start Y] + m * ([end Y] - [start Y]);
                z = [start Z] + m * ([end Z] - [start Z]);
            }else if (outcodeOut & TOPAABB){        // Point is above box (z-axis)
                m = ([max Z] - [start Z]) / ([end Z] - [start Z]);
                x = [start X] + m * ([end X] - [start X]);
                y = [start Y] + m * ([end Y] - [start Y]);
                z = [max Z];
            }else{                              // Point is below box (z-axis)
                m = ([min Z] - [start Z]) / ([end Z] - [start Z]);
                x = [start X] + m * ([end X] - [start X]);
                y = [start Y] + m * ([end Y] - [start Y]);
                z = [min Z];
            }
            
            // Now we move outside point to intersection point to clip
            // and get ready for next pass.
            if (outcodeOut == outcode0) {
                [start setX:x];
                [start setY:y];
                [start setZ:z];
                outcode0 = ComputeOutCode(start, min, max);
            } else {
                [end setX:x];
                [end setY:y];
                [end setZ:z];
                outcode1 = ComputeOutCode(end, min, max);
            }
        }
    }
    if (accept) {
        return true;
    }else{
        return false;
    }
}

OutCode ComputeOutCode(MVector * p, MVector * min, MVector * max)
{
    OutCode code;
    
    code = INSIDEAABB;          // initialised as being inside of clip window
    
    if ([p X] < [min X])        // to the left of box
        code |= LEFTAABB;
    else if ([p X] > [max X])   // to the right of box
        code |= RIGHTAABB;
    if ([p Y] < [min Y])        // nearer that the box
        code |= NEARAABB;
    else if ([p Y] > [max Y])   // farther than the box
        code |= FARAABB;
    if ([p Z] < [min Z])        // Below bottom of box
        code |= BOTTOMAABB;     
    else if([p Z] > [max Z])    // Above top of box
        code |= TOPAABB;        
    
    return code;
}

@end
