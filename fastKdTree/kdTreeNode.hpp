//
//  kdTreeNode.hpp
//  sarcastic
//
//  Created by Darren on 23/10/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#ifndef kdTreeNode_hpp
#define kdTreeNode_hpp

#include <stdio.h>
#include "TriangleMesh.hpp"
#include "AABB.hpp"
#include "splitCandidate.hpp"

TriangleMesh globalMesh;

typedef struct packedNode {
    int splitDim;
    float splitPos;
    int leftIndex;
    int rghtIndex;
    int nTriangles;
    int *triIndices;
} packedNode ;

class kdTreeNode {
    
public:
    int size;
    int leftAddress;                // Array index of leftchild
    AABB                            aabb ;
    std::vector<AABB>               triAABBs ;
    std::vector<splitCandidate>     splitList ;
    std::vector<int>                triangles ;     // index of triangle in global Mesh
    std::vector<int>                triangleMask;   // bit mask indicating if triangles[x] is in this node

    bool    isLeaf=false;
    float   splitPos;
    int     dim;        // Dimension of split axis, 0,1,2 = x,y,z
    int     level;
    
    kdTreeNode(){};
    kdTreeNode(std::vector<int> tris):triangles(tris) {} ;
    kdTreeNode(std::string plyFileName) ;
    
    AABB boundingVolume();
    void medianSplit(kdTreeNode &left, kdTreeNode &rght);
    void split(int dim, float pos, kdTreeNode &left, kdTreeNode &rght) ;
    
};

#endif /* kdTreeNode_hpp */
