/***************************************************************************
 * 
 *           Module :  shadowRays.cpp
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 15/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      routine to generate an array of 'shadow rays' that travel back
 *      to the sensor
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

#include "shadowRays.hpp"
#include <sarclib/sarclib.h>

void buildShadowRays(int                 nRays,              // The number of reflected rays being considered
                     SPVector            RxPos,              // The Receiver location in x,y,z
                     Ray                 *reflectedRays,     // Array of reflected rays - used for their origin as its the reflection point to Rx
                     Ray                 *shadowRays,        // Output - shadow rays to be tested for occlusion later
                     double              *ranges             // Output - the range to the receiver for each shadowray. Calculated here to save calcs later
)
{
    
    SPVector dir, origin, rxp ;
    double unitvect_tmp, rng ;
    
    for(int ind=0; ind<nRays; ++ind){
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
        shadowRays[ind].pol = reflectedRays[ind].pol ;
        shadowRays[ind].id  = reflectedRays[ind].id  ;

    }
    return ;
}
