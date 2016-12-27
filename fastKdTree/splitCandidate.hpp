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

class splitCandidate {
    
public:
    
    float pos;                              // position for this event
    int dim;                                // dimension of this splitcandidate
    int owner;                              // index of owning triangle
    std::vector<unsigned char>  leftTris ;  // Array mask of triangles to left of this event
    std::vector<unsigned char>  rghtTris ;  // Array mask of triangles to right of this event
    
    splitCandidate(){};
    splitCandidate(float pos, int o, int dim, int ntris): pos(pos), owner(o), dim(dim){ leftTris.reserve(ntris); rghtTris.reserve(ntris);};
};

#endif /* splitCandidate_hpp */
