/***************************************************************************
 * 
 *           Module :  polygonIntersection.c
 *          Program :  sarclib
 *       Created by :  Darren Muff on 26/07/2013
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Functions to clip polygons using another polygon as a 'cookie-cutter'
 *      All these functions are based on 2-dimensional polygons.
 *      For efficency reasons a polygon is limited to having 8 vertices
 *
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

#include "polygonIntersection.h"

/// Function that determines whether a coordinate c is to the left
/// or right of a line l. It works by using the sign of the determinant
/// of Vectors (AB,AM) where M(X,Y) is the query point.
/// The Position, left or right, is given by:
///     leftorright = sign( (Bx-Ax)*(Y-Ay) - (By-Ay)*(X-Ax) )
///
enum side leftright(Line l, Coordinate c){
    
    Coordinate point1 = l.point1;
    Coordinate point2 = l.point2;
    float sign = ((point2.x-point1.x)*(c.y-point1.y)) - ((point2.y-point1.y)*(c.x-point1.x));
    if(sign == 0){
        return(centre);
    }else if (sign < 0){
        return(rightinside);
    }else{
        return(leftoutside);
    }
}

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
Coordinate intersectionPoint(Line line, Coordinate start, Coordinate end){
    Line l;
    l.point1 = start ;
    l.point2 = end ;
    
    Coordinate p1,p2,p3,p4;
    float x1,x2,x3,x4,y1,y2,y3,y4;
    float d,n1,n2,ua;
    float xans,yans;
    p1 = line.point1;
    p2 = line.point2;
    p3 = l.point1;
    p4 = l.point2;
    x1 = p1.x;
    x2 = p2.x;
    x3 = p3.x;
    x4 = p4.x;
    y1 = p1.y;
    y2 = p2.y;
    y3 = p3.y;
    y4 = p4.y;
    d  = ((y4-y3)*(x2-x1)) - ((x4-x3)*(y2-y1));
    n1 = ((x4-x3)*(y1-y3)) - ((y4-y3)*(x1-x3));
    n2 = ((x2-x1)*(y1-y3)) - ((y2-y1)*(x1-x3));
    ua = n1/d;
    //ub = n2/d;
    
    xans = x1 + ua*(x2-x1);
    yans = y1 + ua*(y2-y1);
    Coordinate ans;
    ans.x = xans;
    ans.y = yans;
    return(ans);
    
}

/// Function to add a 2 dimensional coordinate 'c' to a polygon 'p'
///
void addToPoly(Polygon *p, Coordinate c){
    if(p->nVerts == 8){
        printf("Error: request for extra vertices: %d\n",p->nVerts);
        return ;
    }
    p->coordinates[p->nVerts].x = c.x;
    p->coordinates[p->nVerts].y = c.y;
    p->nVerts++;
    return ;
}

/// Function to calculate the intersection polygon when a polygon is used
/// as a cookie-cutter for another polygon
/// 
Polygon intersection(Polygon p1, Polygon p2){
    
    Polygon outputPoly = p2;
    int nVertices = p1.nVerts;
    
    for(int i=0; i<nVertices; i++){
        Coordinate c1 = p1.coordinates[i];
        Coordinate c2 = p1.coordinates[(i+1)%nVertices];
        Line clipEdge;
        clipEdge.point1 = c1;
        clipEdge.point2 = c2;
        Polygon inputPoly = outputPoly;
        outputPoly.nVerts=0;
        if(inputPoly.nVerts == 0){
            return inputPoly;
        }
        Coordinate S;
        S = inputPoly.coordinates[inputPoly.nVerts-1];
        for(int j=0; j<inputPoly.nVerts; j++){
            Coordinate E ;
            E = inputPoly.coordinates[j];
            if ( leftright(clipEdge,E) == rightinside){
                if ( leftright(clipEdge,S) != rightinside){
                    Coordinate c3 = intersectionPoint(clipEdge,S, E);
                    addToPoly(&outputPoly, c3);
                }
                addToPoly(&outputPoly, E);
            }else if ( leftright(clipEdge, S) == rightinside){
                Coordinate c4 = intersectionPoint(clipEdge, S, E);
                addToPoly(&outputPoly, c4);
            }
            S = E;
        }
    }
    
    return(outputPoly);
}

