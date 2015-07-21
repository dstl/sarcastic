//
//  trianglecpp.h
//  sarcastic
//
//  Created by Darren on 05/07/2015.
//  Copyright (c) 2015 Dstl. All rights reserved.
//

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
