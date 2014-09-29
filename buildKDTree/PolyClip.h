/***************************************************************************
 *
 *       Module:    PolyClip.h
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

#ifndef Sadilac_PolyClip_h
#define Sadilac_PolyClip_h

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

typedef struct Poly {
    int nVerts;
    Coordinate coordinates[8]; // Max points for two intersecting quadrilaterals
} Poly;

enum side leftright(Line l, Coordinate c) ;
Coordinate intersectionPoint(Line line, Coordinate start, Coordinate end) ;
void addToPoly(Poly *p, Coordinate c) ;
Poly intersection(Poly p1, Poly p2) ;

#endif
