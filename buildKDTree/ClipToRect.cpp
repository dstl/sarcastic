/***************************************************************************
 *  CohenSutherland.cpp
 *  Sadilac
 *
 *  Created by Muff Darren on 04/06/2012.
 *  Copyright (c) 2012 [dstl]. All rights reserved.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   04/06/2012
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

#include <iostream>
#include "ClipToRect.h"

// Use cohen-sutherland algorithm
// http://en.wikipedia.org/wiki/Cohen-Sutherland_algorithm
extern "C" {
    
    OutCode ComputeOutCode(float x, float y, float xmin, float xmax, float ymin, float ymax);
    
    // Compute the bit code for a point (x, y) using the clip rectangle
    // bounded diagonally by (xmin, ymin), and (xmax, ymax)
    
    OutCode ComputeOutCode(float x, float y, float xmin, float xmax, float ymin, float ymax)
    {
        OutCode code;
        
        code = INSIDE;          // initialised as being inside of clip window
        
        if (x < xmin)           // to the left of clip window
            code |= LEFT;
        else if (x > xmax)      // to the right of clip window
            code |= RIGHT;
        if (y < ymin)           // below the clip window
            code |= BOTTOM;
        else if (y > ymax)      // above the clip window
            code |= TOP;
        
        return code;
    }
    
    // Cohen–Sutherland clipping algorithm clips a line from
    // P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with 
    // diagonal from (xmin, ymin) to (xmax, ymax).
    void ClipToRect(float *x0, float *y0, float *x1, float *y1, float xmin, float xmax, float ymin, float ymax, int *status)
    {
        // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
        OutCode outcode0 = ComputeOutCode(*x0, *y0, xmin, xmax, ymin, ymax);
        OutCode outcode1 = ComputeOutCode(*x1, *y1, xmin, xmax, ymin, ymax);
        bool accept = false;
        
        while (true) {
            if (!(outcode0 | outcode1)) { // Bitwise OR is 0. Trivially accept and get out of loop
                accept = true;
                break;
            } else if (outcode0 & outcode1) { // Bitwise AND is not 0. Trivially reject and get out of loop
                break;
            } else {
                // failed both tests, so calculate the line segment to clip
                // from an outside point to an intersection with clip edge
                double x, y;
                
                // At least one endpoint is outside the clip rectangle; pick it.
                OutCode outcodeOut = outcode0 ? outcode0 : outcode1;
                
                // Now find the intersection point;
                // use formulas y = y0 + slope * (x - x0), x = x0 + (1 / slope) * (y - y0)
                if (outcodeOut & TOP) {           // point is above the clip rectangle
                    x = *x0 + (*x1 - *x0) * (ymax - *y0) / (*y1 - *y0);
                    y = ymax;
                } else if (outcodeOut & BOTTOM) { // point is below the clip rectangle
                    x = *x0 + (*x1 - *x0) * (ymin - *y0) / (*y1 - *y0);
                    y = ymin;
                } else if (outcodeOut & RIGHT) {  // point is to the right of clip rectangle
                    y = *y0 + (*y1 - *y0) * (xmax - *x0) / (*x1 - *x0);
                    x = xmax;
                } else {   // point is to the left of clip rectangle
                    y = *y0 + (*y1 - *y0) * (xmin - *x0) / (*x1 - *x0);
                    x = xmin;
                }
                
                // Now we move outside point to intersection point to clip
                // and get ready for next pass.
                if (outcodeOut == outcode0) {
                    *x0 = x;
                    *y0 = y;
                    outcode0 = ComputeOutCode(*x0, *y0, xmin, xmax, ymin, ymax);
                } else {
                    *x1 = x;
                    *y1 = y;
                    outcode1 = ComputeOutCode(*x1, *y1, xmin, xmax, ymin, ymax);
                }
            }
        }
        if (accept) {
            return ;
        }else{
            *status = NOINTERSECTION ;
            return ;
        }
    }
}
