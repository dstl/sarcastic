/** @file********************************************************************
 *
 *       Module:    polygonIntersection.h
 *      Program:    sarclib
 *   Created by:    Darren Muff on 26/07/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *      Functions to clip polygons using another polygon as a 'cookie-cutter'
 *      All these functions are based on 2-dimensional polygons.
 *      For efficency reasons a polygon is limited to having 8 vertices
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  14/03/2013
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
 *
 ***************************************************************************/

#ifndef sarclib_polygonIntersection_h
#define sarclib_polygonIntersection_h

#include "sarclib.h"

#define NOINTERSECTION NULL

enum side {centre=0,rightinside= -1,leftoutside=1};

typedef struct Coordinate {
    float x ;
    float y ;
} Coordinate;

typedef struct Line {
    Coordinate point1;
    Coordinate point2;
} Line;

typedef struct Polygon {
    int nVerts;                 ///< Number of vertices (max is 8)
    Coordinate coordinates[8];  ///< Max points for two intersecting quadrilaterals
} Polygon;

/// Function that determines whether a coordinate c is to the left
/// or right of a line l. It works by using the sign of the determinant
/// of Vectors (AB,AM) where M(X,Y) is the query point.
/// The Position, left or right, is given by:
///     leftorright = sign( (Bx-Ax)*(Y-Ay) - (By-Ay)*(X-Ax) )
///
enum side leftright(Line l, Coordinate c);

/// Function to find the intersection point where a line defined by a start point and an end point
/// intersect with a line.
///
/// It works like this:
///
/// The equations of the lines are
/// Pa = P1 + ua ( P2 - P1 )
/// Pb = P3 + ub ( P4 - P3 )
///
/// Solving for the point where Pa = Pb gives the following two
/// equations in two unknowns (ua and ub)
///
/// x1 + ua (x2 - x1) = x3 + ub (x4 - x3)
/// and
/// y1 + ua (y2 - y1) = y3 + ub (y4 - y3)
///
/// Solving gives the following expressions for ua and ub
///
/// d  = (y4-y3)(x2-x1) - (x4-x3)(y2-y1)
/// n1 = (x4-x3)(y1-y3) - (y4-y3)(x1-x3)
/// n2 = (x2-x1)(y1-y3) - (y2-y1)(x1-x3)
/// ua = n1/d
/// ub = n2/d
///
/// Substituting either of these into the corresponding equation for
/// the line gives the intersection point. For example the intersection point (x,y) is
///
/// x = x1 + ua (x2 - x1)
/// y = y1 + ua (y2 - y1)
///
Coordinate intersectionPoint(Line line, Coordinate start, Coordinate end);

/// Function to add a 2 dimensional coordinate 'c' to a polygon 'p'
///
void addToPoly(Polygon *p, Coordinate c);

/// Function to calculate the intersection polygon when a polygon is used
/// as a cookie-cutter for another polygon
///
Polygon intersection(Polygon p1, Polygon p2);

#endif
