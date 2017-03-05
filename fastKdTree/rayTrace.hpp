//
//  rayTrace.hpp
//  sarcastic
//
//  Created by Darren Muff on 05/03/2017.
//  Copyright © 2017 Dstl. All rights reserved.
//

#ifndef rayTrace_hpp
#define rayTrace_hpp

#include <stdio.h>
#include "buildTree.hpp"

void  rayTrace(TriangleMesh *mesh, kdTree::KdData **kdTree, int *numNodesInTree) ;

#endif /* rayTrace_hpp */
