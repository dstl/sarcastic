//***************************************************************************
//
//  randomRays.cl
//  GOSS
//
//  Created by Darren Muff on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//  OpenCL kernel to generate a random gaussian distribution of rays. The
//  rays 
//
// CLASSIFICATION        :  UNCLASSIFIED
// Date of CLASSN        :  26/02/2014
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


//! Represents the state of a particular generator
typedef struct{ uint x; uint c; } mwc64x_state_t;

enum{ MWC64X_A = 4294883355U };
enum{ MWC64X_M = 18446383549859758079UL };

float box_muller(ulong seed, float m, float s, float *y2, float *use_last);
float ranf(mwc64x_state_t *s);
void MWC64X_Step(mwc64x_state_t *s);
uint MWC64X_NextUint(mwc64x_state_t *s);

__kernel void randomRays(const ulong  seed,         // Random number generator seed
                         const double xStdev,       // x standard deviation of ray distribution
                         const double yStdev,       // y standard deviation of ray distribution
                         const SPVector txPos,      // Transmit location of rays
                         const SPVector aimpoint,   // aimpoint (mean) of beam
                         const int nx,              // number of rays in x
                         const int ny,              // number of rays in y
                         const double pow,          // Ray Power (Pp = (Pt * Gtx / (4*PI)) * dAz * dEl)
                         __global Ray *rayArray)    // output of rays
{
    
    float xang, yang;
    SPVector zHat, elAxis, azAxis, azRayDir, elRayDir, aimDir, rayDir ;
    float use_last = 0;
    float y2;
    
    int xId = get_global_id(0);
    int yId = get_global_id(1);
    
    if (xId >=0 && xId < nx && yId >=0 && yId < ny) {

        // Calculate beam az and el axes for this aim direction
        //
        VECT_SUB(aimpoint, txPos, aimDir);
        VECT_CREATE(0, 0, 1, zHat);
        VECT_CROSS(aimDir, zHat,   elAxis);
        VECT_CROSS(elAxis, aimDir, azAxis);
        xang = box_muller(seed, 0, xStdev, &y2, &use_last);
        yang = box_muller(seed, 0, yStdev, &y2, &use_last);
        vectRotateAxis(aimDir, azAxis, xang, &azRayDir);
        vectRotateAxis(azRayDir, elAxis, yang, &elRayDir);
        VECT_NORM(elRayDir, rayDir) ;
        rayArray[yId*nx+xId].dir = rayDir;
        rayArray[yId*nx+xId].org = txPos ;
        rayArray[yId*nx+xId].pow = pow ;
        rayArray[yId*nx+xId].len = 0;
    }
    return ;
}

float box_muller(ulong rand, float m, float s, float *y2, float *use_last)	/* normal random variate generator */
{                                   /* mean m, standard deviation s */
	float x1, x2, w, y1;
//	static float y2;
//	static int use_last = 0;
    mwc64x_state_t seed={get_global_id(0)+rand,get_global_id(1)+rand};
    
	if (*use_last)		        /* use value from previous call */
	{
		y1 = *y2;
		*use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * ranf(&seed) - 1.0;
			x2 = 2.0 * ranf(&seed) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
        
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		*y2 = x2 * w;
		*use_last = 1;
	}
    
	return( m + y1 * s );
}

//float ranf2(ulong randoms){
//    ulong seed = randoms + get_global_id(1) * get_global_size(0) + get_global_id(0) ;
//    seed = (seed * 0x5DEECE66DL + 0xBL) & ((1L << 48) - 1);
//    uint iresult = seed >> 16 ;
//    return (iresult  / 4294967296.0f);
//}

float ranf(mwc64x_state_t *s){
    float rn = (float)MWC64X_NextUint(s);
    float ans = rn / 4294967296;
    return (ans);
}

//! Return a 32-bit integer in the range [0..2^32)
uint MWC64X_NextUint(mwc64x_state_t *s)
{
	uint res=s->x ^ s->c;
	MWC64X_Step(s);
	return res;
}

void MWC64X_Step(mwc64x_state_t *s)
{
	uint X=s->x, C=s->c;
	uint Xn=MWC64X_A*X+C;
	uint carry=(uint)(Xn<C);				// The (Xn<C) will be zero or one for scalar
	uint Cn=mad_hi(MWC64X_A,X,carry);
	
	s->x=Xn;
	s->c=Cn;
}
