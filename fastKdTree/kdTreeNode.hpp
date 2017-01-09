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

extern TriangleMesh globalMesh;
extern int * kdTreeTriangleIndicesOutput ;

typedef union {
    struct KdTreeLeaf {
        unsigned int flagDimAndOffset;
        // bits 0..1        : splitting dimension
        // bits 2..30       : offset bits to 'leaf list'
        // bits 31 (sign)   : flag whether node is a leaf
        float splitPosition;
    } leaf;
    
    struct KdTreeBranch {
        unsigned int flagDimAndOffset;
        // bits 0..30       : offset to first child
        // bit 31 (sign)    : flag whether node is a leaf
        float splitPosition;
    } branch;
} KdData ;

extern KdData * kdTreeOutput ;

class kdTreeNode {
    
public:
    int size=0;                       // Number of triangles in this node +1 (used to calc memory reqs)
    int leftAddress=0;                // Array index of leftchild node
    int triangleIndex=0 ;             // index of first triangle
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
    kdTreeNode(std::vector<int> tris) ;
    kdTreeNode(std::string plyFileName) ;
    
    AABB BVforAllTris();
    void medianSplit(kdTreeNode &left, kdTreeNode &rght);
    void split(int dim, float pos, kdTreeNode &left, kdTreeNode &rght) ;
    
};

#endif /* kdTreeNode_hpp */
