/***************************************************************************
 *
 *       Module:    trianglecpp.h
 *      Program:    SARCASTIC
 *   Created by:    Darren Muff on 21/07/2015.
 *                  Copyright (c) 2015 Dstl. All rights reserved.
 *
 *   Description:
 *      Class header for trianglecpp class
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  21/07/2015
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

#ifndef __sarcastic__trianglecpp__
#define __sarcastic__trianglecpp__

#include <iostream>
#include <stdio.h>
#include <SIlib/SIlib.h>

class Triangle {
   
    
public:
    Triangle(SPVector aa, SPVector bb, SPVector cc) ;
    Triangle(SPVector aa, SPVector bb, SPVector cc, int material) ;
    Triangle(SPVector aa, SPVector bb, SPVector cc, std::string material) ;
    ~Triangle();
    
    void setMaterial(int matId);
    void setMaterial(std::string material);
    
    SPVector   AA ;
    SPVector   BB ;
    SPVector   CC ;
    SPVector   NN ;
    double     area ;
    int        matId;
    double     globalToLocalMat[9];
    double     localToGlobalMat[9];
};


#endif /* defined(__sarcastic__trianglecpp__) */
