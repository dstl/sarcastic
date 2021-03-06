/***************************************************************************
 * 
 *           Module :  AABB.hpp
 *          Program :  fastKdTree
 *       Created by :  Darren Muff on 05/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  04-Nov-2018
 *   Description:
 *      Programme to build a K-Dimensional tree quicker and more scalably than the pevious implementation
 *      This version uses the approach detailed in [1]. The previous approach uses that of [2]
 *
 *      1. Zhou, Kun, et al. "Real-time kd-tree construction on graphics hardware."
 *         ACM Transactions on Graphics (TOG) 27.5 (2008): 126.
 *
 *      2. Wald, Ingo, and Vlastimil Havran. "On building fast kd-trees for ray tracing,
 *         and on doing that in O (N log N)." Interactive Ray Tracing 2006, IEEE Symposium
 *         on. IEEE, 2006.
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

#ifndef AABB_hpp
#define AABB_hpp

#include <stdio.h>
#include <sarclib/sarclib.h>
#include <vector>
#include <algorithm>

typedef int OutCode;

#define INSIDEAABB  0   // 000000
#define LEFTAABB    1   // 000001
#define RIGHTAABB   2   // 000010
#define NEARAABB    4   // 000100
#define FARAABB     8   // 001000
#define TOPAABB     16  // 010000
#define BOTTOMAABB  32  // 100000

OutCode ComputeOutCode(SPVector * p, SPVector * min, SPVector * max);

class AABB {
    
public:
    SPVector AA ;   // Lower coord of box
    SPVector BB ;   // Upper coord of box
    
    AABB(){} ;
    AABB(SPVector min, SPVector max): AA(min), BB(max) {}
    AABB(std::vector<SPVector> points) ;
    AABB(std::vector<SPVector> points, int maxDim, float boundMin, float boundMax);
    
    void clipToTriangle(SPVector vertA, SPVector vertB, SPVector vertC, double splitPos, int splitDim, AABB &left, AABB &right);
    float surfaceArea();
    void print();
    bool isEmpty();
    
};



#endif /* AABB_hpp */
