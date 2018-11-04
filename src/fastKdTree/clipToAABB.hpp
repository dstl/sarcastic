/***************************************************************************
 * 
 *           Module :  clipToAABB.hpp
 *          Program :  fastKdTree
 *       Created by :  Darren Muff on 09/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  04-Nov-2018
 *      Description :  Header for functions that clip a three dimensional ray to a
 *                     axis-aligned bounding box (AABB). The algorithm
 *                     is an extension of the Cohen-Sutherland line splitting
 *                     algorithm used in computer graphics. [1]
 *
 *             [1] https://en.wikipedia.org/wiki/Cohen–Sutherland_algorithm
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

#ifndef clipToAABB_hpp
#define clipToAABB_hpp

#include <stdio.h>
#include <sarclib/sarclib.h>
#include "AABB.hpp"

typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000
#define NOINTERSECTION -1


int     clipToAABB    (AABB boundingBox, SPVector *lineStart, SPVector *lineEnd);
void    ClipToBox     (SPVector *p0,     SPVector *p1,        SPVector min, SPVector max, int *status);
OutCode ComputeOutCode(SPVector p,       SPVector min,        SPVector max);

#endif /* clipToAABB_hpp */
