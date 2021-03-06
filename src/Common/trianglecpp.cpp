/***************************************************************************
 * 
 *           Module :  trianglecpp.cpp
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 21/07/2015
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *      Class header for trianglecpp class
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

#include "trianglecpp.h"
#include "materialProperties.h"
extern "C"
{
#include "matrixMultiplication.h"
}

Triangle::Triangle(){}

Triangle::Triangle(SPVector aa, SPVector bb, SPVector cc){
    AA = aa ;
    BB = bb ;
    CC = cc ;
    matId = 0;
    
    // Calculate Normal. Calculated here in case its incorrect in collada file
    //
    SPVector ab; VECT_SUB(bb, aa, ab);
    SPVector bc; VECT_SUB(cc, bb, bc);
    SPVector x;  VECT_CROSS(ab, bc, x);
    double l ; l = VECT_MAG(x);
    area = l * 0.5;
    VECT_SCMULT(x, (1/l), NN);
    d = VECT_DOT(AA, NN) ;
    
    SPVector zhat; VECT_CREATE(0, 0, 1, zhat);
    double alpha; alpha = atan2(NN.y, NN.x);
    double beta; beta = acos(VECT_DOT(zhat, NN));
    double T_dash[9], T_dashdash[9];
    
    T_dash[0]     = cos(alpha);
    T_dash[1]     = sin(alpha);
    T_dash[2]     = 0;
    T_dash[3]     = -sin(alpha);
    T_dash[4]     = cos(alpha);
    T_dash[5]     = 0;
    T_dash[6]     = 0;
    T_dash[7]     = 0;
    T_dash[8]     = 1;
    T_dashdash[0] = cos(beta);
    T_dashdash[1] = 0;
    T_dashdash[2] = -sin(beta);
    T_dashdash[3] = 0;
    T_dashdash[4] = 1;
    T_dashdash[5] = 0;
    T_dashdash[6] = sin(beta);
    T_dashdash[7] = 0;
    T_dashdash[8] = cos(beta);
    
    matmul(T_dashdash, T_dash, globalToLocalMat, 3, 3, 3, 3);
    mat3by3inv(globalToLocalMat, localToGlobalMat);
    return ;
}

Triangle::Triangle(SPVector aa, SPVector bb, SPVector cc, int material) {
    
    *this = Triangle(aa, bb, cc) ;
    setMaterial(material) ;
    return ;
}

Triangle::Triangle(SPVector aa, SPVector bb, SPVector cc, std::string material) {
    
    *this = Triangle(aa, bb, cc) ;
    setMaterial(material) ;
    return ;
}

Triangle::~Triangle(){
    
}

void Triangle::setMaterial(int mat){ matId = mat; }

void Triangle::setMaterial(std::string material){
  
    std::transform(material.begin(), material.end(),material.begin(), ::toupper);
    
    for(int i=0; i<NMATERIALS; i++){
        if (material.find(std::string(globalMatProps[i].matname)) != std::string::npos ) {
            matId = i ;
            i = NMATERIALS ;
        }
    }
    return ;
}

void Triangle::print(){
    std::cout << AA.x << " " << AA.y << " " << AA.z << std::endl;
    std::cout << BB.x << " " << BB.y << " " << BB.z << std::endl;
    std::cout << CC.x << " " << CC.y << " " << CC.z << std::endl;
    std::cout << AA.x << " " << AA.y << " " << AA.z << std::endl;

    return ;
    
}
