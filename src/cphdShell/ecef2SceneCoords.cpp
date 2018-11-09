/***************************************************************************
 * 
 *           Module :  ecef2SceneCoords.cpp
 *          Program :  cphdShell
 *       Created by :  Darren Muff on 02/08/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *      Description :  Program to convert a set of ECEF coordniate to 
 *                     ENU (East North Up) centred on a lat,;pn,alt
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
 ****************************************************************************/

#include "ecef2SceneCoords.hpp"

void ecef2SceneCoords(int nPoints, SPVector *points, SPVector geoPoint){
    SPVector point;
    SPVector Z, E, N, A, v ;
    double e,n,a ;
    
    VECT_CREATE(0, 0, 1., Z) ;           // Create Z-axis - Earth axis of rotation
    VECT_CROSS(Z, geoPoint, E) ;         // Create local East vector - tangent to Earth's surface
    VECT_NORM(E, E) ;                    // and make it a unit vector
    VECT_CROSS(geoPoint, E, N);          // Create local North vector - tangent to Earth's surface
    VECT_NORM(N, N) ;                    // and make it a unit vector
    VECT_NORM(geoPoint, A) ;             // Create local height vector
 //    double graz = DEG2RAD(0);
    for (int i=0; i<nPoints; i++){
        point = points[i];
        VECT_SUB(point, geoPoint, v) ;   // v is vector from point in space to geopoint
        e = VECT_DOT(v, E) ;             // e is component of v in local East direction
        n = VECT_DOT(v, N) ;             // n is component of v in local North direction
        a = VECT_DOT(v, A) ;             // a is component of v in altitude direction
 //        a= fabs(e) * tan(graz);
 //        if(i==0)printf("grazing angle is %fdeg a is %f\n",RAD2DEG(atan2(a,fabs(e))),a);
 //        n=0;
        VECT_CREATE(e, n, a, points[i]) ;
    }
   
    return ;
}




