/***************************************************************************
 *  boxMuller_Random.c
 *  triDecomp
 *
 *  Created by Muff Darren on 27/09/2014.
 *  Copyright (c) 2014 [dstl]. All rights reserved.
 *
 * Original source code 
 *  (c) Copyright 1994, Everett F. Carter Jr.
 *  Permission is granted by the author to use
 *  this software for any application provided this
 *  copyright notice is preserved.
 *
 *   Description:
 *      Implements the Polar form of the Box-Muller Transformation
 *      This module is taken from http://www.taygeta.com/random/gaussian.html
 *      There is a useful discussion on generating Gaussian random numbers
 *      at http://www.design.caltech.edu/erik/Misc/Gaussian.html
 *
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   25/09/2014
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
 ***************************************************************************/

#include "boxMullerRandom.h"

extern float ranf();         /* ranf() is uniform in 0..1 */

float box_muller(float m, float s)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	float x1, x2, w, y1;
	static float y2;
	static int use_last = 0;
    
	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * Ranf() - 1.0;
			x2 = 2.0 * Ranf() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
        
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
    
	return( m + y1 * s );
}
