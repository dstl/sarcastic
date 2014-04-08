//***************************************************************************
//
//  readKdTree.h
//  Sadilac
//
//  Created by Darren on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//
//
// CLASSIFICATION        :  UNCLASSIFIED
// Date of CLASSN        :  02/08/2012
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

#ifndef __GOSS__readKdTree__
#define __GOSS__readKdTree__
#include "sarcastic.h"

void readKdTree(const char *filename,           // Filename for KdTree
                    AABB * SceneBoundingBox,    // Bounding box for scene
                    int * nTriangles,           // Number of triangles in scene
                    Triangle ** Triangles,      // Array containing triangles
                    int * nTextures,            // number if textures in scene
                    Texture ** textures,        // array containing textures
                    int * nLeaves,              // number of leaf nodes in KdTree
                    int *** triangleLists,      // Array of int arrays containing triangles for each leaf
                    int * nTreeNodes,           // Number of nodes in KdTree
                    KdData ** KdTree,           // KdTree returned.
                    TriCoords **tricos );       // Triangle coordinates

#endif
