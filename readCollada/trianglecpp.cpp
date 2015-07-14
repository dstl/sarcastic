//
//  trianglecpp.cpp
//  sarcastic
//
//  Created by Darren on 05/07/2015.
//  Copyright (c) 2015 Dstl. All rights reserved.
//

#include "trianglecpp.h"
#include "materialProperties.h"
extern "C"
{
#include "matrixMultiplication.h"
}

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
    
    // Convert material to be upper case
    // 'just in case ;-)'
    //
    char c[material.length()];
    for(int i=0; i<material.length(); i++)c[i] = std::toupper(material[i]);
    
    for(int i=0; i<NMATERIALS; i++){
        if (!strcmp(materialProperties[i].matname,c)) {
            matId = i ;
            i = NMATERIALS ;
        }
    }
    return ;
}