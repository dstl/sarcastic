//
//  splitCandidate.cpp
//  sarcastic
//
//  Created by Darren on 27/10/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#include "splitCandidate.hpp"

splitCandidate::splitCandidate(float pos, int o, int dim, int ntris): pos(pos), owner(o), dim(dim), ntris(ntris) {
    
    leftTris = new unsigned char [ntris] ;
    rghtTris = new unsigned char [ntris] ;
}

splitCandidate::~splitCandidate() {
    delete [] leftTris ;
    delete [] rghtTris ;
}

// Copy constructor - required as constructor dynamically allocates memory
//
splitCandidate::splitCandidate(const splitCandidate &split) {
    pos = split.pos ;
    owner = split.owner ;
    dim = split.dim  ;
    ntris = split.ntris ;
    leftTris = new unsigned char [ntris] ;
    memcpy(leftTris, split.leftTris, ntris*sizeof(unsigned char)) ;
    rghtTris = new unsigned char [ntris] ;
    memcpy(rghtTris, split.rghtTris, ntris*sizeof(unsigned char)) ;
}
