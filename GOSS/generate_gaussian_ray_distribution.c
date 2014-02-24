/***************************************************************************
 *
 *       Module:    generate_gaussian_ray_distribution.c
 *      Program:    GOSS
 *   Created by:    Darren on 01/08/2013.
 *                  Copyright (c) 2013 Dstl. All rights reserved.
 *
 *   Description:
 *      function to generate random distribution of ray directions
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  22/02/2014
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

#include <stdio.h>
#include "GOSS.h"
#include <stdlib.h>
#define BIGNUMBER ((int) 16384) //

float ranf(){
    double rand = arc4random_uniform(16384) ;
    double r = rand / (double)16384;
    return r;
}

float box_muller(float m, float s)	/* normal random variate generator */
{                                   /* mean m, standard deviation s */
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
			x1 = 2.0 * ranf() - 1.0;
			x2 = 2.0 * ranf() - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );
        
		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}
    
	return( m + y1 * s );
}

// Function to generate a gaussian distribution of rays that originate from 'position'
// xStdev, yStdev are the standard deviation of the beam in x and y
// txPos is the transmit position for the beam
// GRP is the ground reference point
// nx,ny are the number of rays required to be stored in rayArray (which must be allocated first)
//
void generate_gaussian_ray_distribution(double xStdev,double yStdev, SPVector txPos, SPVector GRP, int nx,int ny, Ray *rayArray){
   
    float xang, yang;
    SPVector zHat, elAxis, azAxis, azRayDir, elRayDir, aimDir, rayDir ;
    
    // Calculate beam az and el axes for this aim direction
    //
    VECT_SUB(GRP, txPos, aimDir);
    VECT_CREATE(0, 0, 1, zHat);
    VECT_CROSS(aimDir, zHat,   elAxis);
    VECT_CROSS(elAxis, aimDir, azAxis);
    for (int y=0; y<ny; y++) {
        for (int x=0; x<nx; x++) {
            xang = box_muller(0, xStdev);
            yang = box_muller(0, yStdev);
            vectRotateAxis(aimDir, azAxis, xang, &azRayDir);
            vectRotateAxis(azRayDir, elAxis, yang, &elRayDir);
            VECT_NORM(elRayDir, rayDir)
            rayArray[y*nx+x].dir = rayDir;
            rayArray[y*nx+x].org = txPos ;
        }
    }
}


