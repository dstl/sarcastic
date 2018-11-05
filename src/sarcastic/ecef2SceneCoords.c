/***************************************************************************
 * 
 *           Module :  ecef2SceneCoords.c
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 02/08/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *    converts an array of coordinates that are in the ECEF system into
 *    local coordinates for a scene. geoPoint is ceneter of the scene and
 *    is specified in ECEF coordinates.
 *    Note scene coordinates are assumed to be X is East and Y is North
 *    with Z being altitude.
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

#include "ecef2SceneCoords.h"
//#define OVERRIDEGRAZINGANGLEWITH ((double)(5.4))

void ecef2SceneCoords(int nPoints, SPVector *points, SPVector geoPoint){
    SPVector point;
    SPVector Z, E, N, A, v ;
    double e,n,a ;
    
    VECT_CREATE(0, 0, 1., Z) ;          // Create Z-axis - Earth axis of rotation
    VECT_CROSS(Z, geoPoint, E) ;        // Create local East vector - tangent to Earth's surface
    VECT_NORM(E, E) ;                   // and make it a unit vector
    VECT_CROSS(geoPoint, E, N);         // Create local North vector - tangent to Earth's surface
    VECT_NORM(N, N) ;                   // and make it a unit vector
    VECT_NORM(geoPoint, A) ;            // Create local height vector
#ifdef OVERRIDEGRAZINGANGLEWITH
    double graz = DEG2RAD(OVERRIDEGRAZINGANGLEWITH);
#endif
    for (int i=0; i<nPoints; i++){
        point = points[i];
        VECT_SUB(point, geoPoint, v) ;  // v is vector from point in space to geopoint
        e = VECT_DOT(v, E) ;            // e is component of v in local East direction
        n = VECT_DOT(v, N) ;            // n is component of v in local North direction
        a = VECT_DOT(v, A) ;            // a is component of v in altitude direction
#ifdef OVERRIDEGRAZINGANGLEWITH
        a= fabs(e) * tan(graz);
#endif
//        if(i==0)printf("grazing angle is %fdeg a is %f\n",RAD2DEG(atan2(a,fabs(e))),a);
//        n=0;
        VECT_CREATE(e, n, a, points[i]) ;
    }
   
    return ;
}




