//
//  kdTreeNode.cpp
//  sarcastic
//
//  Created by Darren on 23/10/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#include "kdTreeNode.hpp"
TriangleMesh globalMesh;
int * kdTreeTriangleIndicesOutput ;
KdData * kdTreeOutput ;
int numOfTriangleIndices ;
int numOfKdTreeNodes ;

kdTreeNode::kdTreeNode(){
    size = 0 ;
    leftAddress = 0;
    triangleIndex = 0 ;
    isLeaf = 0 ;
    dim = 0 ;
    level = 0 ;
    smallroot = -1 ;
    smallntris = 0 ;
    triangleMask = NULL;
}
kdTreeNode::kdTreeNode(std::vector<int> tris):triangles(tris) {
    aabb = BVforAllTris() ;
    level = 0;
    smallroot = -1;
    smallntris = 0;
}

kdTreeNode::kdTreeNode(std::string plyFileName)
// Initialise the root KdTreeNode using a .plyFile
//
{
    globalMesh.readPLYFile(plyFileName);
    globalMesh.checkIntegrityAndRepair();
    globalMesh.buildTriangleAABBs();
    
    triAABBs = globalMesh.AABBs ;
    for(int i=0; i< globalMesh.triangles.size(); ++i){
        triangles.push_back(i) ;
    }
    
    aabb = BVforAllTris() ;
    level = 0;
    smallroot = -1;
    smallntris = 0;
}

// Copy constructor
//
kdTreeNode::kdTreeNode(const kdTreeNode &node){
    size = node.size ;
    leftAddress = node.leftAddress ;
    triangleIndex = node.triangleIndex ;
    aabb = node.aabb ;
    for (int i=0; i<node.triAABBs.size(); ++i){
        triAABBs.push_back(node.triAABBs[i]);
    }
    for (int i=0; i<node.splitList.size(); ++i){
        splitList.push_back(node.splitList[i]) ;
    }
    for (int i=0; i<node.triangles.size(); ++i){
        triangles.push_back(node.triangles[i]);
    }
    isLeaf = node.isLeaf ;
    splitPos = node.splitPos ;
    dim = node.dim ;
    level = node.level ;
    smallroot = node.smallroot ;
    smallntris = node.smallntris ;
    
    if (smallntris > 0) {
        triangleMask = new unsigned char [smallntris] ;
        memcpy(triangleMask, node.triangleMask, smallntris * sizeof(unsigned char)) ;
    }
}

kdTreeNode::~kdTreeNode(){
    if (smallntris > 0) {
        delete [] triangleMask ;
    }
}

AABB kdTreeNode::BVforAllTris()
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
    if (smallroot == -1) {
        printf("Error: split() function only works on small nodes. Did you mean medianSplit()?\n");
        exit(1) ;
    }
    
    splitPos = pos;
    dim = k;
    
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
    left.triangleMask = new unsigned char [smallntris] ;
    rght.triangleMask = new unsigned char [smallntris] ;
    
//    left.triangleMask = (unsigned char *)sp_malloc(smallntris * sizeof(unsigned char) ) ;
//    rght.triangleMask = (unsigned char *)sp_malloc(smallntris * sizeof(unsigned char) ) ;
    
    for(int t=0; t<smallntris; ++t){
        left.triangleMask[t] = triangleMask[t] & (triAABBs[t].AA.cell[k] <= pos) ;
        rght.triangleMask[t] = triangleMask[t] & (triAABBs[t].BB.cell[k] >= pos) ;
    }
    
    left.splitList  = splitList  ;
    left.triAABBs   = triAABBs   ;
    left.smallroot  = smallroot  ;
    left.smallntris = smallntris ;
    rght.splitList  = splitList  ;
    rght.triAABBs   = triAABBs   ;
    rght.smallroot  = smallroot  ;
    rght.smallntris = smallntris ;
    
    // Set level for the children
    //
    left.level = level+1;
    rght.level = level+1;
    
    return ;
}

void kdTreeNode:: medianSplit(kdTreeNode &left, kdTreeNode &rght)
// medianSplit() splits the node at the median point of the largest axis
// and returns two child nodes that have their AABB's set and the triAABB's correctly clipped to the split plane.
// In addition the triangles array in each child correctly holds the index of each triangle
//
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
    pos = aabb.AA.cell[maxAxis]+(len[maxAxis] /2.0) ;
    
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
        if(triAABBs[t].AA.cell[maxAxis] >= pos ){
            rght.triangles.push_back(triIdx);
            rght.triAABBs.push_back(triAABBs[t]) ;
        }
        if(triAABBs[t].BB.cell[maxAxis] <= pos){
            left.triangles.push_back(triIdx);
            left.triAABBs.push_back(triAABBs[t]) ;
        }
        if( (triAABBs[t].AA.cell[maxAxis] < pos) && (triAABBs[t].BB.cell[maxAxis] > pos) )
        {
            // Triangle straddles the split
            //
            vertA = globalMesh.vertices[globalMesh.triangles[triIdx].a].asSPVector() ;
            vertB = globalMesh.vertices[globalMesh.triangles[triIdx].b].asSPVector() ;
            vertC = globalMesh.vertices[globalMesh.triangles[triIdx].c].asSPVector() ;
            
            triAABBs[t].clipToTriangle(vertA, vertB, vertC, pos, maxAxis, ablft, abrgt);
            
            if(ablft.surfaceArea() > 0.0001 ){
                left.triangles.push_back(triIdx);
                left.triAABBs.push_back(ablft) ;
            }
            if(abrgt.surfaceArea() > 0.0001 ){
                rght.triangles.push_back(triIdx);
                rght.triAABBs.push_back(abrgt) ;
            }
    
        }
    }
    
    return;
}

void printKdTreeNodes(std::vector<kdTreeNode> nodelist){
    
    if(nodelist[0].size == 0)return ;
    
    printf("Triangle Indexing Information\n");
    printf("-----------------------------\n");
    for(int i=0; i<nodelist[0].size; ++i){
        printf("[%03d]  %02d\n",i,kdTreeTriangleIndicesOutput[i]);
    }
    printf("KdTree\n");
    printf("------\n");
    for(int i=0; i<nodelist.size(); ++i){
        printf("[%02d]",i);
        for(int j=0;j<nodelist[i].level;++j){
            printf("-");
        }
        if(nodelist[i].isLeaf){
            printf("X");
            printf(" <%02d> #%02d = [",nodelist[i].triangleIndex,kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex]);
            for(int k=0; k<kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex]; ++k){
                printf(" %02d",kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex+1+k]);
            }
            printf(" ]");
        }else{
            printf("|");
            printf("  [%02d][%02d]",nodelist[i].leftAddress,nodelist[i].leftAddress+1);
        }
        printf("\n");
    }
    return;
}

void printKdTreeData(){
    
    KdData * node;
    
    printf("Triangle Indexing Information\n");
    printf("-----------------------------\n");
    for(int i=0; i<numOfTriangleIndices; ++i){
        printf("[%03d]  %02d\n",i,kdTreeTriangleIndicesOutput[i]);
    }
    printf("KdTree\n");
    printf("------\n");
    for(int i=0; i<numOfKdTreeNodes; ++i){
        
        node = &(kdTreeOutput[i]) ;
        printf("[%02d]",i);
        if(KDT_ISLEAF( node ) ){
            printf("X");
            printf(" <%02d> #%02d = [",KDT_OFFSET(node),kdTreeTriangleIndicesOutput[KDT_OFFSET(node)]);
            for(int k=0; k<kdTreeTriangleIndicesOutput[KDT_OFFSET(node)]; ++k){
                printf(" %02d",kdTreeTriangleIndicesOutput[KDT_OFFSET(node)+1+k]);
            }
            printf(" ]");
        }else{
            printf("|");
            printf("  [%02d][%02d]",KDT_OFFSET(node),KDT_OFFSET(node)+1);
        }
        printf("\n");
    }

}

