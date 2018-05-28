/***************************************************************************
 *
 *       Module:    sceneExtent
 *      Program:    SARCASTIC
 *   Created by:    Darren on 28/05/2018.
 *                  Copyright (c) 2018 Dstl. All rights reserved.
 *
 *   Description:
 *      function to calculate the maximum and minimum extent of a scene when viewed
 *      from a specific location
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  21/04/2017
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


#ifndef sceneExtent_h
#define sceneExtent_h
#include <SIlib2/SIlib2.h>
#include "TriangleMesh.hpp"
void sceneExtent(SPVector Pos, TriangleMesh &mesh, double &maxEl, double &maxAz, double &minEl, double &minAz, AABB &sceneBoundingBox);

#endif /* sceneExtent_h */
