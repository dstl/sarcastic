/***************************************************************************
 * 
 *           Module :  reflect.cpp
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 25/01/2015
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      Function to calculate the reflection rays from an input of
 *      hit locations
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

#include "reflect.hpp"
#include "materialProperties.h"

void reflect(int    nRays,              // Number of rays to reflect
             Ray    *rays,              // Array of rays to consider
             Hit    *hits,              // Array of hit points for each ray
             ATS    *accelTriangles,    // Array of triangles in 'accelerated' format
             Ray    *reflectedRays      // output array of reflected rays
)
{
    int t, k, ku, kv ;
    unsigned int modulo[5];
    ATS T ;
    SPVector N, v, R, I, hp, perpol;
    float ks;
    
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    
    for(int ind=0; ind<nRays; ++ind){
        t = hits[ind].trinum;
        T = accelTriangles[t];
        k = T.k;
        ku = modulo[k+1];
        kv = modulo[k+2];
        N.cell[k] = 1;
        N.cell[ku] = T.nd_u;
        N.cell[kv] = T.nd_v;
        VECT_NORM(N,N);
        
        I = rays[ind].dir;
        if(VECT_DOT(I,N)>0){
            VECT_MINUS(N,N);
        }
        VECT_SCMULT(N, 2.0 * VECT_DOT(I, N), v);
        VECT_SUB(I, v, R);
        
        VECT_SCMULT(rays[ind].dir, hits[ind].dist, hp);
        
        VECT_ADD(rays[ind].org, hp, reflectedRays[ind].org);
        reflectedRays[ind].dir = R ;
        reflectedRays[ind].len = rays[ind].len + hits[ind].dist ;
        reflectedRays[ind].id  = rays[ind].id ;

        
        // Now calculate forward scattered ray power by only considering the
        // specular component of the reflective surface texture.
        // (Diffuse and shinyness components are only taken into account on rays
        // returning to the sensor from a visible hit point
        // Note that this is because we use PO for last bounce but phong shding for all others
        //
        ks =  globalMatProps[T.textureInd].specular ;
        // power at reflected point is Pt*Gtx / 4 PI R^2
        //
        double islf = (4.0 * 3.1415926536 * hits[ind].dist * hits[ind].dist);
        islf = (islf < 1 ) ? 1 : islf ;
        
        // Calculate the polarisation of the reflected ray based upon the polarisation
        // of the incident ray
        //
        VECT_CROSS(rays[ind].pol, N, perpol);
        
        // The reflected ray is calculated using Geometrical optics. PO is only used
        // on the shadowrays back to the receiver. A diasadvantage of this is that
        // the reflected power does not take into account the polarisation. To approximate
        // the reduction in E-field intensity for a parallel polarised E-field multiply it
        // by the sine of the angle between the normal and the polarisation.
        //
        double sintheta = VECT_MAG(perpol);
        reflectedRays[ind].pow = rays[ind].pow * ks * sintheta * sintheta ;
        VECT_CROSS(perpol, reflectedRays[ind].dir, reflectedRays[ind].pol);
        VECT_NORM(reflectedRays[ind].pol,reflectedRays[ind].pol);

    }
    
    return ;
}
