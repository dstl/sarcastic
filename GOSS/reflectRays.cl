//***************************************************************************
//
//  reflectRays.cl
//  GOSS
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//  OpenCl kernel code to reflect an array of rays. Inputs are
//  the rays and an array of hitpoints (which contains the triangle id
//  for this hit) and an array of triangles
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

__kernel void reflect(__global Triangle * Triangles,
                      __global Texture  * Textures,
                      int nTextures,
                      __global Ray *rays,
                      __global Hit *hits,
                      __global Ray *reflectedRays,
                      int nRays
                      )
{
    int ind, t, k, ku, kv ;
    unsigned int modulo[5];
    Triangle T ;
    SPVector N, v, R, I, hp;
 
    modulo[0] = 0; modulo[1] = 1; modulo[2] = 2; modulo[3] = 0; modulo[4]=1;
    ind = get_global_id(0) ;
    
    if (ind >=0 && ind < nRays ) {
        
        t = hits[ind].trinum;
        T = Triangles[t];
        k = T.k;
        ku = modulo[k+1];
        kv = modulo[k+2];
        N.cell[k] = 1;
        N.cell[ku] = T.nd_u;
        N.cell[kv] = T.nd_v;
        VECT_NORM(N,N);
        
        I = rays[ind].dir;
        VECT_SCMULT(N, 2.0 * VECT_DOT(I, N), v);
        VECT_SUB(I, v, R);
//        VECT_NORM(R,R) ;
        
        VECT_SCMULT(rays[ind].dir, hits[ind].dist, hp);
        
        VECT_ADD(rays[ind].org, hp, reflectedRays[ind].org);
        reflectedRays[ind].dir = R ;
        reflectedRays[ind].len = rays[ind].len + hits[ind].dist ;
        
        
        // Now calculate forward scattered ray power by only considering the
        // specular component of the reflective surface texture.
        // (Diffuse and shinyness components are only taken into account on rays
        // returning to the sensor from a visible hit point
        //
        reflectedRays[ind].pow = rays[ind].pow * Textures[T.textureInd].ks ;
//        printf("power in %e, power out %e\n",rays[ind].pow,reflectedRays[ind].pow);
    }
    
    return ;
}

