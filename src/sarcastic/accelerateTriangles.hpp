/***************************************************************************
 *
 *       Module:    accelerateTriangles.hpp
 *      Program:    SARCASTIC
 *   Created by:    Darren on 15/03/2017.
 *                  Copyright (c) 2017 Dstl. All rights reserved.
 *
 *   Description:
 *      Function to convert triangle coordinates into barycentric coordinates
 *      in order to accelerate the ray-triangle intersection calculations
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  15/03/2017
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

#ifndef accelerateTriangles_hpp
#define accelerateTriangles_hpp

#include <stdio.h>
#include "TriangleMesh.hpp"

typedef struct ATS {
    int  triNum;    // Triangle ID
    double d;         // Constant of plane equation
    double nd_u;      // Normal.u / normal.k
    double nd_v;      // normal.v / normal.k
    int k;          // projection dimension
    double kbu;
    double kbv;
    double kbd;
    double kcu;
    double kcv;
    double kcd;
    int textureInd;
} ATS;

void accelerateTriangles(TriangleMesh *mesh, ATS **accelTriangles) ;

#endif /* accelerateTriangles_hpp */
