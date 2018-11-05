/***************************************************************************
 * 
 *           Module :  clipToAABB.cl
 *          Program :  kernels
 *       Created by :  Darren Muff on 02/08/2012
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *  OpenCl code for Cohen–Sutherland clipping algorithm that clips a line from
 *  P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with
 *  diagonal from (xmin, ymin) to (xmax, ymax).
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
 ****************************************************************************/

#ifndef CLIPTOAABB_CL
#define CLIPTOAABB_CL
#include "SPVector.cl"
#include "structures.cl"

typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000
#define NOINTERSECTION -1
#define TRUE 1
#define FALSE 0

int     clipToAABB    (AABB boundingBox, SPVector *lineStart, SPVector *lineEnd);
void    ClipToBox     (SPVector *p0,     SPVector *p1,        SPVector min, SPVector max, int *status);
OutCode ComputeOutCode(SPVector p,       SPVector min,        SPVector max);

int clipToAABB(AABB boundingBox, SPVector *lineStart, SPVector *lineEnd){
    
    int status=0;
    ClipToBox(lineStart, lineEnd, boundingBox.AA, boundingBox.BB, &status);
    if (status == NOINTERSECTION) return 0;     // No intersect in any dim means no intersection with volume
    
    return 1;
    
}

// Cohen–Sutherland clipping algorithm clips a line from
// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with
// diagonal from (xmin, ymin) to (xmax, ymax).
void ClipToBox(SPVector *p0, SPVector *p1, SPVector min, SPVector max, int *status)
{
    // compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
    OutCode outcode0 = ComputeOutCode(*p0, min, max);
    OutCode outcode1 = ComputeOutCode(*p1, min, max);
    int accept = FALSE;
    
    while (TRUE) {
        if (!(outcode0 | outcode1)) { // Bitwise OR is 0. Trivially accept and get out of loop
            accept = TRUE;
            break;
        } else if (outcode0 & outcode1) { // Bitwise AND is not 0. Trivially reject and get out of loop
            break;
        } else {
            // failed both tests, so calculate the line segment to clip
            // from an outside point to an intersection with clip edge
            double x, y, z, m;
            
            // At least one endpoint is outside the clip rectangle; pick it.
            OutCode outcodeOut = outcode0 ? outcode0 : outcode1;
            
            // Now find the intersection point;
            
            if (outcodeOut & FARAABB) {             // Point is beyond the box (y-axis)
                m = (max.y - p0->y) / (p1->y - p0->y);
                x = p0->x + m * (p1->x - p0->x) ;
                y = max.y;
                z = p0->z + m * (p1->z - p0->z) ;
            }else if (outcodeOut & NEARAABB){     // Point is before the box (y-axis)
                m = (min.y - p0->y) / (p1->y - p0->y);
                x = p0->x + m * (p1->x - p0->x) ;
                y = min.y;
                z = p0->z + m * (p1->z - p0->z) ;
            }else if (outcodeOut & LEFTAABB){       // Point is left of box (x-axis)
                m = (min.x - p0->x) / (p1->x - p0->x);
                x = min.x;
                y = p0->y + m * (p1->y - p0->y);
                z = p0->z + m * (p1->z - p0->z);
            }else if (outcodeOut & RIGHTAABB){      // Point is right of box (x-axis)
                m = (max.x - p0->x) / (p1->x - p0->x);
                x = max.x;
                y = p0->y + m * (p1->y - p0->y);
                z = p0->z + m * (p1->z - p0->z);
            }else if (outcodeOut & TOPAABB){        // Point is above box (z-axis)
                m = (max.z - p0->z) / (p1->z - p0->z);
                x = p0->x + m * (p1->x - p0->x);
                y = p0->y + m * (p1->y - p0->y);
                z = max.z;
            }else{                              // Point is below box (z-axis)
                m = (min.z - p0->z) / (p1->z - p0->z);
                x = p0->x + m * (p1->x - p0->x);
                y = p0->y + m * (p1->y - p0->y);
                z = min.z;
            }
            
            // Now we move outside point to intersection point to clip
            // and get ready for next pass.
            if (outcodeOut == outcode0) {
                p0->x = x;
                p0->y = y;
                p0->z = z;
                outcode0 = ComputeOutCode(*p0, min, max);
            } else {
                p1->x = x;
                p1->y = y;
                p1->z = z;
                outcode1 = ComputeOutCode(*p1, min, max);
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

// Compute the bit code for a point (x, y) using the clip rectangle
// bounded diagonally by (xmin, ymin), and (xmax, ymax)

OutCode ComputeOutCode(SPVector p, SPVector min, SPVector max)
{
    OutCode code;
    
    code = INSIDEAABB;      // initialised as being inside of clip window
    
    if (p.x < min.x)        // to the left of box
        code |= LEFTAABB;
    else if (p.x > max.x)   // to the right of box
        code |= RIGHTAABB;
    if (p.y < min.y)        // nearer that the box
        code |= NEARAABB;
    else if (p.y > max.y)   // farther than the box
        code |= FARAABB;
    if (p.z < min.z)        // Below bottom of box
        code |= BOTTOMAABB;
    else if(p.z > max.z)    // Above top of box
        code |= TOPAABB;
    
    return code;
}
#endif