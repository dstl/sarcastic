//
//  splitCandidate.hpp
//  sarcastic
//
//  Created by Darren on 27/10/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#ifndef splitCandidate_hpp
#define splitCandidate_hpp

#include <stdio.h>
#include <vector>
#include <SIlib2/SIlib2.h>

class splitCandidate {
    
public:
    
    float pos;                              // position for this event
    int dim;                                // dimension of this splitcandidate
    int owner;                              // index of owning triangle
    int ntris;
    unsigned char  *leftTris ;              // Array mask of triangles to left of this event
    unsigned char  *rghtTris ;              // Array mask of triangles to right of this event
    
    splitCandidate(){};
    splitCandidate(float pos, int o, int dim, int ntris);
    ~splitCandidate();
    splitCandidate(const splitCandidate &split) ;   // Copy constructor
};

#endif /* splitCandidate_hpp */
