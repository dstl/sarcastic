/***************************************************************************
 * 
 *           Module :  buildRays.hpp
 *          Program :  sarcastic
 *       Created by :  Darren Muff on 15/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *      Description :
 *      routine to generate an array of rays ready to be cast at a scene
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

#ifndef buildRays_hpp
#define buildRays_hpp

#include <stdio.h>
#include "AABB.hpp"
#include "rayTrace.hpp"
#include "threadCore.hpp"

enum  BUILDRAYMETHOD {
        TRIANGLECENTRE  = 1,    // 1 - each ray aimed at triangle centre
        RANDOMRAYS      = 2,    // 2 - random rays on each call across scene
        FIRSTTIMERANDOM = 3,    // 3 - random rays created first time but the same hitpoints used for each subsequent call
        PARALLELRANDOM  = 4     // 4 - like 2 (random rays on each call across the scene) but rays are parallel from Tx
} ;


void buildRays(Ray **rayArray, int *nRays, int nAzRays, int nElRays, TriangleMesh *mesh, SPVector TxPos,
               double PowPerRay, AABB SceneBoundingBox,SPVector **rayAimPoints, int method, int Pol);

#endif /* buildRays_hpp */
