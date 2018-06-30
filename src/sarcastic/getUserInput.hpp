/***************************************************************************
 *
 *       Module:    getUserInput.hpp
 *      Program:    SARCASTIC
 *   Created by:    Darren on 15/03/2017.
 *                  Copyright (c) 2017 Dstl. All rights reserved.
 *
 *   Description:
 *      Simple routone to interrogate the user and see if teh programme can get
 *      some sense out of them regarding what they would like to achieve
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

#ifndef getUserInput_hpp
#define getUserInput_hpp

#include <stdio.h>
#include "TriangleMesh.hpp"
#include <sarclib/sarclib.h>
#include "threadCore.hpp"

int getUserInput(CPHDHeader *hdr, TriangleMesh *baseMesh, TriangleMesh *moverMesh, char **outCPHDFile,
                 int *startPulse, int *nPulses,
                 int *bounceToShow, int *nAzBeam, int *nElBeam,
                 int *interrogate, SPVector *interogPt, double *interogRad, int *interogX, int *interogY,
                 FILE **interrogateFP, int *pulseUndersampleFactor, int *polarisation, int *rayGenMethod, SPStatus *status) ;

#endif /* getUserInput_hpp */
