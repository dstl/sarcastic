//
//  kdTreeNode.cpp
//  sarcastic
//
//  Created by Darren on 23/10/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#include "kdTreeNode.hpp"

kdTreeNode::kdTreeNode(std::string plyFileName)
// Initialise the root KdTreeNode using a .plyFile
//
{
    globalMesh.readPLYFile(plyFileName);
    globalMesh.checkIntegrityAndRepair();
    globalMesh.buildTriangleAABBs();
    
    leftChild = NULL;
    rghtChild = NULL;
    triAABBs = globalMesh.AABBs ;
    for(int i=0; i< globalMesh.triangles.size(); ++i){
        triangles.push_back(i) ;
    }
    level = 0;
}

AABB kdTreeNode::boundingVolume()
{
    SPVector triaa, tribb;
    AABB ans;
    
    VECT_CREATE( 9e9,  9e9,  9e9, ans.AA);
    VECT_CREATE(-9e9, -9e9, -9e9, ans.BB);
    
    for(int i=0; i<triangles.size(); ++i){
        triaa = triAABBs[i].AA ;
        tribb = triAABBs[i].BB ;
        
        ans.AA.x = (triaa.x < ans.AA.x) ? triaa.x : ans.AA.x ;
        ans.AA.y = (triaa.y < ans.AA.y) ? triaa.y : ans.AA.y ;
        ans.AA.z = (triaa.z < ans.AA.z) ? triaa.z : ans.AA.z ;
        
        ans.BB.x = (tribb.x > ans.BB.x) ? tribb.x : ans.BB.x ;
        ans.BB.y = (tribb.y > ans.BB.y) ? tribb.y : ans.BB.y ;
        ans.BB.z = (tribb.z > ans.BB.z) ? tribb.z : ans.BB.z ;
    }
    return ans;
}

void kdTreeNode::split(int k, float pos, kdTreeNode &left, kdTreeNode &rght)
{
    // Set the AABBs for the kids
    //
    left.aabb.AA = aabb.AA ;
    left.aabb.BB = aabb.BB ;
    left.aabb.BB.cell[k] = pos ;
    rght.aabb.AA = aabb.AA ;
    rght.aabb.BB = aabb.BB ;
    rght.aabb.AA.cell[k] = pos ;
    
    // Sort triangles from this node into child nodes using bit mask
    //
    left.triangleMask.reserve(triangleMask.size());
    rght.triangleMask.reserve(triangleMask.size());
    
    for(int t=0; t<triangles.size(); ++t){
        
        left.triangleMask[t] = triangleMask[t] & (triAABBs[t].AA.cell[k] <= pos) ;
        rght.triangleMask[t] = triangleMask[t] & (triAABBs[t].BB.cell[k] >= pos) ;
        
    }
    
    left.splitList = splitList ;
    left.triAABBs  = triAABBs  ;
    rght.splitList = splitList ;
    rght.triAABBs  = triAABBs  ;
    
    // Set level for the children
    //
    left.level = level+1;
    rght.level = level+1;
    
    return ;
}

void kdTreeNode:: medianSplit(kdTreeNode &left, kdTreeNode &rght)
{
    
    // Set level for the children
    //
    left.level = level+1;
    rght.level = level+1;
    
    // find longest axis
    //
    double len[3], pos;
    int maxAxis ;
    len[0] = aabb.BB.x - aabb.AA.x ;
    len[1] = aabb.BB.y - aabb.AA.y ;
    len[2] = aabb.BB.z - aabb.AA.z ;
    if(len[0] > len[1]){
        if(len[0] > len[2]){ maxAxis = 0;}else{ maxAxis = 2;}
    }else if (len[1] > len[2]){
        maxAxis = 1;
    }else{
        maxAxis = 2;
    }
    pos = aabb.AA.cell[maxAxis]+(len[maxAxis] /2) ;
    
    splitPos = pos ;
    dim      = maxAxis;
    
    // Set the AABBs for the kids by splitting this node
    // down the middle
    //
    left.aabb.AA = aabb.AA ;
    left.aabb.BB = aabb.BB ;
    left.aabb.BB.cell[maxAxis] = pos ;
    
    rght.aabb.AA = aabb.AA;
    rght.aabb.BB = aabb.BB;
    rght.aabb.AA.cell[maxAxis] = pos ;
    
    // sort the triangles in the mesh in this node into the mesh for each of the child nodes.
    //
    
    SPVector vertA,vertB,vertC;
    AABB ablft, abrgt;
    int triIdx ;
    
    for (int t=0; t<triangles.size(); ++t){
        triIdx = triangles[t] ;
        if(triAABBs[t].AA.cell[maxAxis] > pos ){
            rght.triangles.push_back(triIdx);
            rght.triAABBs.push_back(triAABBs[t]) ;
        }else if(triAABBs[t].BB.cell[maxAxis] <= pos){
            left.triangles.push_back(triIdx);
            left.triAABBs.push_back(triAABBs[t]) ;
        }else{
            // Triangle straddles the split
            //
            left.triangles.push_back(triIdx);
            rght.triangles.push_back(triIdx);
            
            vertA = globalMesh.vertices[globalMesh.triangles[triIdx].a].asSPVector() ;
            vertB = globalMesh.vertices[globalMesh.triangles[triIdx].b].asSPVector() ;
            vertC = globalMesh.vertices[globalMesh.triangles[triIdx].c].asSPVector() ;
            
            triAABBs[t].clipToTriangle(vertA, vertB, vertC, pos, maxAxis, ablft, abrgt);
            
            left.triAABBs.push_back(ablft) ;
            rght.triAABBs.push_back(abrgt) ;
        }
    }
    
    return;
}



