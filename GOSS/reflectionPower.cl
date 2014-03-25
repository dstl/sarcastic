//***************************************************************************
//
//  reflectPower.cl
//  GOSS
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//  OpenCl kernel code to calculate reflection power for an
//  array of rays back at the receiver
//
//
// CLASSIFICATION        :  UNCLASSIFIED
// Date of CLASSN        :  02/08/2012
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
// USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
// SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
// IRELAND.
//
//***************************************************************************
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#include "SPVector.cl"
#include "structures.cl"

__kernel void reflectPower(__global Triangle * Triangles,       // Array of Triangles. Each triangle references a texture
                           __global Texture  * Textures,        // Array of Textures containing specular, diffuse adn shinyness constants
                           __global Hit *hits,                  // Array of hit locations to x-ref with triangles (and then Textures) for material props
                           const SPVector RxPos,                // Location of Receiver in x,y,z
                           const double GrxOverFourPi,          // Receiver antenna gain / 4Pi.
                           int nRays,
                           __global Ray *rays,                  // unit vector rays arriving at hitpoint
                           __global Ray *Rrays,                 // Reflection unit vector rays
                           __global Ray *Vrays,                 // Viewpoint unit vector. The vector direction to calculate power for (usually to receiver)
                           __global double *ranges,             // Range to receiver for each shadow ray (precalculated in shadowRay generation)
                           __global rangeAndPower *rnp           // Output array of ray power at, and range to reciever
                           )
{
    int ind, k, ku, kv ;
    Triangle T ;
    Texture tex ;
    SPVector N, R, L, V;
    double is, kd, ks, n, LdotN, RdotV, pow, rng;
    /* double id; */
    unsigned int modulo[5];
    
    ind = get_global_id(0) ;
    
    if (ind >=0 && ind < nRays ) {
        modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;

        T   = Triangles[hits[ind].trinum];
        tex = Textures[T.textureInd];
        
        // Create surface normal for taxture associated with this hitpoint
        //
        k  = T.k;
        ku = modulo[k+1];
        kv = modulo[k+2];
        
        N.cell[k]  = 1;
        N.cell[ku] = T.nd_u;
        N.cell[kv] = T.nd_v;
        VECT_NORM(N,N);
        
        // Set up Phong shading parameters for this texture
        //
        /*ia = 0 ;  */          // Ambient RF power. Assumed 0 for RF frequencies
        /*id = rays[ind].pow;*/ // Diffuse component incident energy on hitpoint
        is = rays[ind].pow ;    // Specular component of incident energy on hitpoint
        kd = tex.kd ;           // Diffuse constant of hitpoint's texture
        ks = tex.ks ;           // Specular constant of hitpoint's texture
        n  = tex.n ;            // Shinyness constant of hitpoints texture
        
        VECT_MINUS(rays[ind].dir, L); // L is the unit vector illumination ray (ray point to source of illumination
        R = Rrays[ind].dir;
        V = Vrays[ind].dir;
        
        LdotN = VECT_DOT(L,N);
        RdotV = VECT_DOT(R, V);
        
        LdotN = (LdotN < 0) ? -LdotN : LdotN;   // removes 'one-way mirror' effect
        RdotV = (RdotV < 0) ? 0 : RdotV;        // if ray bouncing away from radar then there is no specular component
        
        // This is the equation for the power but as there is no ambient Rf energy and is = id then we can refine it a bit
        //   pow = (ia * ka ) + (kd * LdotN * id) + ( ks * pow( RdotV , n) * is );
        //
        pow = is * ((kd * LdotN) + ( ks * pow( RdotV, n) )) ;
        
        // Now calculate the amount of power recieved at the receiver, taking into account
        // the inverse square law for return path and receiver antenna gain
        //
        rng = ranges[ind] ;
        rnp[ind].power = pow * GrxOverFourPi / (rng * rng) ;
        rnp[ind].range = rng ;
    }
    return ;
}