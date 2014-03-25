//***************************************************************************
//
//  buildShadowRays.cl
//  GOSS
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//  OpenCl kernel code to build an array of shadow rays from a reflection point
//  back to the receiver. As we have to do this for every ray we can do a bit
//  of optimisation here by also returning the range from the reflection point
// to the receiver. This measn we dont have to calculate it later on.
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

__kernel void buildShadowRays(const    int nRays,           // The number of reflected rays being considered
                              const    SPVector RxPos,      // The Receiver location in x,y,z
                              __global Ray *reflectedRays,  // Array of reflected rays - used for their origin as its the reflection point to Rx
                              __global Ray *shadowRays,     // Output - shadow rays to be tested for occlusion later
                              __global double *ranges       // Output - the range to the receiver for each shadowray. Calculated here to save calcs later
                           )
{
    int ind = get_global_id(0) ;
    
    SPVector dir, origin, rxp ;
    double unitvect_tmp, rng ;
    
    if (ind >=0 && ind < nRays ) {
        
        origin              = reflectedRays[ind].org ;
        rxp                 = RxPos ;
        
        VECT_SUB(rxp, origin, dir);
        
        rng                 = VECT_MOD(dir);
        unitvect_tmp        = 1.0 / rng;
        
        VECT_SCMULT(dir, unitvect_tmp, dir) ;
        
        ranges[ind]         = rng ;
        shadowRays[ind].org = origin ;
        shadowRays[ind].dir = dir ;
        shadowRays[ind].pow = reflectedRays[ind].pow ;
        shadowRays[ind].len = reflectedRays[ind].len ;
        
    }
    return ;
    
}