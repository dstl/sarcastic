/***************************************************************************
 *
 *       Module:    PolyClip.c
 *      Program:    Sadilac
 *   Created by:    Darren on 15/05/2013.
 *                  Copyright (c) 2013 [dstl]. All rights reserved.
 *
 *   Description:
 *   Function to clip a polygon using another polygon as a cookie cutter
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  18/02/2014
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
#include <stdio.h>
#include "PolyClip.h"

enum side leftright(Line l, Coordinate c){
    // Use the sign of the determinant of Vectors (AB,AM) where
    // M(X,Y) is the query point.
    // Position left or right = sign( (Bx-Ax)*(Y-Ay) - (By-Ay)*(X-Ax) )
    Coordinate point1 = l.point1;
    Coordinate point2 = l.point2;
    float sign = ((point2.x-point1.x)*(c.y-point1.y)) - ((point2.y-point1.y)*(c.x-point1.x));
    /*if(sign == 0){
        return(centre);
    }else */if (sign <= 0){
        return(rightinside);
    }else{
        return(leftoutside);
    }
}

Coordinate intersectionPoint(Line line, Coordinate start, Coordinate end){
    Line l;
    l.point1 = start ;
    l.point2 = end ;
    
    /*
     From http://paulbourke.net/geometry/lineline2d/
     
     The equations of the lines are
     Pa = P1 + ua ( P2 - P1 )
     Pb = P3 + ub ( P4 - P3 )
     
     Solving for the point where Pa = Pb gives the following two
     equations in two unknowns (ua and ub)
     
     x1 + ua (x2 - x1) = x3 + ub (x4 - x3)
     and
     y1 + ua (y2 - y1) = y3 + ub (y4 - y3)
     
     Solving gives the following expressions for ua and ub
     
     d  = (y4-y3)(x2-x1) - (x4-x3)(y2-y1)
     n1 = (x4-x3)(y1-y3) - (y4-y3)(x1-x3)
     n2 = (x2-x1)(y1-y3) - (y2-y1)(x1-x3)
     ua = n1/d
     ub = n2/d
     
     Substituting either of these into the corresponding equation for
     the line gives the intersection point. For example the intersection point (x,y) is
     
     x = x1 + ua (x2 - x1)
     y = y1 + ua (y2 - y1)
     */
    
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
void addToPoly(Poly *p, Coordinate c){
    if(p->nVerts == 8){
        //        printf("Error: request for extra vertices: %d\n",p->nVerts);
        return ;
    }
    p->coordinates[p->nVerts].x = c.x;
    p->coordinates[p->nVerts].y = c.y;
    p->nVerts++;
    return ;
}

// Uses polygon p1 to 'cookie-cut' polygon p2.
// The ordering is important if p1 has no area is inside p2 then
// no intersection will be found. If p2 has no area and is inside
// p1 then the correct union will be found
// Polygons must be clockwise
//
Poly intersection(Poly p1, Poly p2){
    
    Poly outputPoly = p2;
    int nVertices = p1.nVerts;
    
    for(int i=0; i<nVertices; i++){
        Coordinate c1 = p1.coordinates[i];
        Coordinate c2 = p1.coordinates[(i+1)%nVertices];
        Line clipEdge;
        clipEdge.point1 = c1;
        clipEdge.point2 = c2;
        Poly inputPoly = outputPoly;
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


