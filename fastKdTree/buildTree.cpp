/***************************************************************************
 *
 *       Module:    buildTree.cpp
 *      Program:    fastKdTree
 *   Created by:    Darren on 05/03/2017.
 *                  Copyright (c) 2017 Dstl. All rights reserved.
 *
 *   Description:
 *      Programme to build a K-Dimensional tree quicker and more scalably than the pevious implementation
 *      This version uses the approach detailed in [1]. The previous approach uses that of [2]
 *
 *      1. Zhou, Kun, et al. "Real-time kd-tree construction on graphics hardware."
 *         ACM Transactions on Graphics (TOG) 27.5 (2008): 126.
 *
 *      2. Wald, Ingo, and Vlastimil Havran. "On building fast kd-trees for ray tracing,
 *         and on doing that in O (N log N)." Interactive Ray Tracing 2006, IEEE Symposium
 *         on. IEEE, 2006.
 *
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  15/01/2014
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 * THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
 * USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
 * SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
 * IRELAND.
 *
 ***************************************************************************/

#include "buildTree.hpp"
TriangleMesh kdTree::globalMesh ;
std::vector<kdTree::kdTreeNode *> kdTree::nodelist ;

void kdTree::buildTree(TriangleMesh *mesh, KdData **kdTree, int *numNodesInTree, TREEOUTPUT printOutput){
    
    kdTreeNode *p, *q;      // used for indexing through the linked list
    
    // Perform the Zhou KdTree build
    //
    treeList *activelist = new treeList ;
    treeList *smalllist  = new treeList ;
    treeList *nextlist   = new treeList ;
    
    // Set up nodelist which is an array of pointers to nodes.
    //
    int nodelistSize = 1;
    for(int i=0; i<=MAXTREELEVEL; ++i){
        nodelistSize *= 2 ;
    }
    nodelist.reserve(nodelistSize) ;
    
    kdTreeNode rootnode(mesh) ;
    
    printf("Input Mesh has %ld triangles\n",globalMesh.triangles.size()) ;
    
    activelist->push_back(&rootnode);
    
    int dpth = 0;
    printf("Processing Large Nodes...\n");
    while(!activelist->empty()){
        
        nextlist->clear() ;
        
        processLargeNodes(&activelist, &smalllist, &nextlist) ;
        
        // Copy activelist into nodelist.
        //
        q = activelist->front();
        while (q){
            p=q ;
            nodelist.push_back(p) ;
            q=p->next ;
        }
        
        printf("[%3d] active list size: %4d\tSmalllist size %4d\tnextlist size: %4d\n",dpth++,activelist->size(),smalllist->size(),nextlist->size());
        swapLists(&nextlist, &activelist) ;
    }
    printf("Done!\n");

    printf("PreProcessing Small Nodes ...\n");
    preProcessSmallNodes(&smalllist);
    
    swapLists(&smalllist, &activelist);
    delete smalllist ;
    
    printf("Done!\n");
    
    printf("Processing Small Nodes...\n");
    dpth=0;
    while(!activelist->empty()){
        
        nextlist->clear() ;
        processSmallNodes(&activelist, &nextlist);
        
        // Copy activelist into nodelist.
        //
        q = activelist->front();
        while (q){
            p=q ;
            nodelist.push_back(p) ;
            q=p->next ;
        }
        
        printf("[%3d] active list size: %4d\tnextlist size: %4d\n",dpth++,activelist->size(),nextlist->size());
        swapLists(&nextlist, &activelist);
    }
    
    if(printOutput >= 4){
        struct stat sb;
        
        if (stat("/tmp/AABBs", &sb) == 0 && S_ISDIR(sb.st_mode)) {
            printf("Directory '/tmp/AABBs' exists. Using it for AABB output\n");
        }else{
            printf("Directory '/tmp/AABBs' does not exist. Creating it...\n");
            if (mkdir("/tmp/AABBs", 0777) != 0 ){
                printf("Failed to create directory '/tmp/AABBs'\n");
                exit(0);
            }
        }
        for (int i=0; i<nodelist.size(); ++i) {
            AABB ab = nodelist[i]->data.aabb ;
            char fn[255] ;
            sprintf(fn, "/tmp/AABBs/aabb_d%02d_%02d.ply",nodelist[i]->data.level,i);
            writeAABBtoPlyFile(ab, std::string(fn));
        }
    }
    
    preOrderTraversalNode(&nodelist, kdTree, numNodesInTree) ;

    if ((printOutput & 0x2) == OUTPUTNODE) {
        printKdTreeNodes(nodelist) ;
    }
    
    if ((printOutput & 0x1) == OUTPUTDATA) {
        printKdTreeData(kdTree, numNodesInTree) ;
    }
    
    if ((printOutput & 0x8) == OUTPUTSUMM) {
        printSummary(nodelist) ;
    }
    
    nodelist.clear() ;
    activelist->erase() ;
    
    return;
}

void kdTree::processLargeNodes(treeList **activelist, treeList **smalllist, treeList **nextlist)
{
    
    // BV4TrisInNode is the AABB that tightly bounds all the triangles in the node
    //
    std::vector<AABB> BV4TrisInNode;
    AABB aabb ;
    // Compute the AABB for each node in activelist
    //
    // As activelist is a ptr to a linked list we can be more efficient at stepping through it
    // than asking for each node at() a given location as the at() function steps through the
    // list from the head each time.
    //
    kdTreeNode *p, *q;
    p = (*activelist)->front() ; // First node in list (head is at index 0)
    q = (*activelist)->front() ; // First node in list (head is at index 0)
    while (q) {
        p = q ;
        aabb = p->BVforAllTris() ;
        BV4TrisInNode.push_back(aabb) ;
        q = p->next ;
    }
    
    // split large nodes
    //
    // Start by clipping away empty space in node around triangles that are in the node
    //
    double s;               // size of node in a given dimansion
    double ldist, hdist;    // low and high distances from the BV containing the triangles and the node BV
    p = (*activelist)->front() ; // First node in list (head is at index 0)
    q = (*activelist)->front() ; // First node in list (head is at index 0)
    int i=0;
    while(q){
        p = q ;
        
        for(int j=0; j<3; ++j){
            s = p->data.aabb.BB.cell[j] - p->data.aabb.AA.cell[j] ;
            aabb = BV4TrisInNode[i] ;
            ldist = aabb.AA.cell[j] - p->data.aabb.AA.cell[j] ;
            if( (ldist/s) > Ce){
                p->data.aabb.AA.cell[j] = aabb.AA.cell[j] ;
            }
            hdist = p->data.aabb.BB.cell[j] - aabb.BB.cell[j] ;
            if( (hdist/s) > Ce){
                p->data.aabb.BB.cell[j] = aabb.BB.cell[j] ;
            }
        }
        
        kdTreeNode *leftNode = new kdTreeNode ;
        kdTreeNode *rghtNode = new kdTreeNode ;
        p->medianSplit(&leftNode, &rghtNode);
        int numLeftChild  = (int)leftNode->data.triangles.size();
        int numRightChild = (int)rghtNode->data.triangles.size();
        p->leftChild = leftNode ;
        p->rghtChild = rghtNode ;
        
        if (numLeftChild < SMALLSIZE) {
            leftNode->smallroot = leftNode ;
            (*smalllist)->push_back(leftNode) ;
        }else{
            (*nextlist)->push_back(leftNode) ;
        }
        if (numRightChild < SMALLSIZE) {
            rghtNode->smallroot = rghtNode ;
            (*smalllist)->push_back(rghtNode) ;
        }else{
            (*nextlist)->push_back(rghtNode) ;
        }
        
        i++;
        q = p->next ;
    }
    
    return ;
}

void kdTree::preProcessSmallNodes(treeList **smalllist)
// This function generates a list of split candidates for each node
// in smalllist.
//
{
    kdTreeNode *p, *q;
    
    p = (*smalllist)->front() ;
    q = (*smalllist)->front() ;
    int i=0;
    while (q){
        p=q ;   // p is the kdTreeNode in smalllist under consideration
        
        // Set the triangle mask to be all '1's to start with showing that each triangle in
        // triangles is actually in this smallnode
        //
        p->data.smallntris = (int)p->data.triangles.size() ;
        p->smallroot = p ;
        
        if(p->data.smallntris <=0){
            printf("error : smallroot with zero triangles in preProcessSmallNodes\n");
            exit(1);
        }
        if(p->smallroot == NULL){
            printf("Error: node %d in activelist is not a small node in function preProcessSmallNodes()\n",i) ;
            exit(1);
        }
        p->data.triangleMask = new unsigned char [p->data.smallntris] ;
        for (int j=0; j<p->data.triangles.size(); ++j) {
            p->data.triangleMask[j] = 1 ;
        }
        // Create split candidates based upon the AABBs of the triangles in this node
        // two for each triangle (low and high), and for each dimension.
        // Store them all in splitList for this node
        //
        for (int k=0; k<3; ++k){
            
            for(int j=0; j<p->data.smallntris; ++j){
                splitCandidate eA(p->data.triAABBs[j].AA.cell[k], j, k, p->data.smallntris) ;
                splitCandidate eB(p->data.triAABBs[j].BB.cell[k], j, k, p->data.smallntris) ;
                
                p->data.splitList.push_back(eA);
                p->data.splitList.push_back(eB);
            }
        }
        for(int j=0; j<p->data.splitList.size(); ++j){
            
            int k = p->data.splitList[j].dim ;
            
            for(int t=0; t<p->data.smallntris; ++t){
                p->data.splitList[j].leftTris[t] = (p->data.triAABBs[t].AA.cell[k] <= p->data.splitList[j].pos) ;
                p->data.splitList[j].rghtTris[t] = (p->data.triAABBs[t].BB.cell[k] >= p->data.splitList[j].pos) ;
            }
        }
        
        q=p->next ;
        i++;
    }
    return ;
}

void kdTree::processSmallNodes(treeList **activelist, treeList **nextlist)
{
    int Cl, Cr ;
    SPVector AA, BB ;
    AABB Vleft,Vrght ;
    float Al,Ar,ProbL,ProbR,SAHp ;
    
    kdTreeNode *p,*q ;
    
    p = (*activelist)->front();
    q = (*activelist)->front();
    int i=0;
    while (q){
        p=q ;
        
        if(p->smallroot == NULL){
            printf("Error: node %d in activelist is not a small node in function processSmallNodes()\n",i) ;
            exit(1);
        }
        
        float A0 = p->data.aabb.surfaceArea() ;
        int minId = -1;
        
        int ntrisInThisNode = reduce(p->data.triangleMask, p->data.smallntris );
        float SAH0 = INTERSECTIONCOST * ntrisInThisNode ;
        float minSAH = SAH0;
        
        for (int j=0; j<p->data.splitList.size(); ++j){
            
            // For this split candidate find the number of triangles to the
            // left and right of the split plane that are also in this node
            //
            int tri = p->data.splitList[j].owner ;
            if(p->data.triangleMask[tri] != 0){
                
                Cl = reduce(p->data.splitList[j].leftTris, p->data.smallntris) ;
                Cr = reduce(p->data.splitList[j].rghtTris, p->data.smallntris) ;
                
                AA = p->data.aabb.AA;
                BB = p->data.aabb.BB;
                AA.cell[p->data.splitList[j].dim] = p->data.splitList[j].pos ;
                BB.cell[p->data.splitList[j].dim] = p->data.splitList[j].pos ;
                Vleft = AABB(p->data.aabb.AA, BB);
                Vrght = AABB(AA, p->data.aabb.BB);
                Al = Vleft.surfaceArea();
                Ar = Vrght.surfaceArea();
                ProbL = Al / A0 ;
                ProbR = Ar / A0 ;
                SAHp = TRAVERSALCOST + ((Cl * ProbL + Cr * ProbR) * INTERSECTIONCOST);
                if(SAHp < minSAH){
                    minSAH = SAHp ;
                    minId = j;
                }
            }
        }
        
        if(minId == -1){
            p->data.isLeaf = true ;
        }else{
            kdTreeNode *leftNode = new kdTreeNode ;
            kdTreeNode *rghtNode = new kdTreeNode ;
            
            p->split(p->data.splitList[minId].dim, p->data.splitList[minId].pos, leftNode, rghtNode) ;
            
            p->leftChild = leftNode ;
            p->rghtChild = rghtNode ;
            
            (*nextlist)->push_back(leftNode);
            (*nextlist)->push_back(rghtNode);
            
        }
        
        
        q=p->next ;
        ++i;
    }
    
    return ;
}

void kdTree::preOrderTraversalNode(std::vector<kdTreeNode *> *nodelist, KdData **kdTreeOutput, int *numOfKdTreeNodes)
{
    // There are two ways to store the tree:
    // Bread-first : the data structure contains all the nodes for each level
    //             : Useful as each level is a multiple of two. Quicker to generate
    // depth-first : The data structure contains the depth of the tree for teh first node before the depth of subsequent nodes
    //             : Quicker to traverse as the data most likely to be accessed next is in nearby memory
    //
    
    // We cant use the bread-first way to store data easily as the process of processing large nodes and then small nodes
    // means that the nodes in nodelist connot be guaranteeed to be in an ascending level/depth order - an earlier node
    // in nodelist could be assigned into small list of it has a few number of triangles which is then added to nodelist
    // later on.
    //
    
    int level ;
    kdTreeNode *p;
    
    // Find the lowest level in the tree
    // and assign the index of each node as the address
    //
    if (! nodelist->empty()) {
        level = 0 ;
        for(int i=0; i<nodelist->size(); ++i){
            p = (*nodelist)[i] ;
            level = ( p->data.level > level) ? p->data.level : level ;
            
        }
    }else{
        printf("Error : Empty nodelist in preOrderTraversalNode. Exiting...\n");
        exit(1);
    }
    
    // Calculate the size of each node in the tree working from the lowest level up
    // As well as provding the size of the kdtree, the info can also be used to
    // store the node info.
    //
    for (int i=level; i>= 0; --i){
        buildSizes(nodelist, i) ;
    }
    
    // Now work down the levels of the tree storing the address of each child in each node
    // (this is a depth-first addressing system
    //
    for (int i=0; i<=level; ++i){
        buildAddresses(nodelist, i);
    }
    
    *numOfKdTreeNodes = nodelist->front()->data.size ;
    *kdTreeOutput = new KdData [*numOfKdTreeNodes] ;
    
    int address ;
    for(int i=0; i<nodelist->size(); ++i){
        p = (*nodelist)[i] ;
        address = p->data.address ;
        // check
        if( address < 0 || address >= *numOfKdTreeNodes){
            printf("Error : calculated address is out of bounds for kdTreeOutput\n");
            exit(1);
        }
        
        //  leafDim (unsigned char)
        //
        //           _______third LSB = flag for leaf 1 = leaf
        //           | ____ 2 LSB = dimension 1=x,2=y or 3=z
        //           | | |
        //           v v v
        // x x x x x 1 1 1
        //
        if (p->data.isLeaf) {
            (*kdTreeOutput)[address].leaf.splitPosition = p->data.splitPos ;
            (*kdTreeOutput)[address].leaf.leafDim = (unsigned char)0x4 ;
            (*kdTreeOutput)[address].leaf.leafDim = ((*kdTreeOutput)[address].leaf.leafDim) | p->data.dim ;
            (*kdTreeOutput)[address].leaf.triangleIndex = i ;
            kdTreeNode *s = p->smallroot ;
            std::vector<int> triangles = s->data.triangles ;
            // check
            if (triangles.size() != p->data.smallntris) {
                printf("uh oh! different number of smallroot triangles in mask\n");
                exit(1);
            }
            int tcount = 0;
            for (int t=0; t<p->data.smallntris; ++t) {
                if (p->data.triangleMask[t] == 1){
                    (*kdTreeOutput)[address+tcount+1].leaf.leafDim = (unsigned char)0x4 ;
                    (*kdTreeOutput)[address+tcount+1].leaf.leafDim = ((*kdTreeOutput)[address].leaf.leafDim) | p->data.dim ;
                    (*kdTreeOutput)[address+tcount+1].leaf.triangleIndex = triangles[t] ;
                    (*kdTreeOutput)[address+tcount+1].leaf.ntriangles = 0 ;
                    tcount++;
                }
            }
            (*kdTreeOutput)[address].leaf.ntriangles =  tcount;
            
        }else{
            (*kdTreeOutput)[address].brch.splitPosition = p->data.splitPos ;
            (*kdTreeOutput)[address].brch.leafDim = (unsigned char)0x0 ;
            (*kdTreeOutput)[address].brch.leafDim = ((*kdTreeOutput)[address].leaf.leafDim)| p->data.dim ;
            (*kdTreeOutput)[address].brch.leftaddress = p->leftChild->data.address ;
            (*kdTreeOutput)[address].brch.rghtaddress = p->rghtChild->data.address ;
        }
        
    }
    
    return;
}

void kdTree::buildSizes(std::vector<kdTreeNode *> *nodelist, int level){
    
    kdTreeNode *p,*l,*r;
    
    for(int i=0; i<nodelist->size(); ++i){
        p = (*nodelist)[i] ;
        if(p->data.level == level){
            if (!p->data.isLeaf) {
                l = p->leftChild;
                r = p->rghtChild;
                p->data.size = l->data.size + r->data.size + 1 ;
            }else{
                p->data.size = reduce(p->data.triangleMask,p->data.smallntris) + 1 ;
            }
        }
    }
    
    return ;
}

void kdTree::buildAddresses(std::vector<kdTreeNode *> *nodelist, int level){
    
    kdTreeNode *p,*l,*r;
    
    for(int i=0; i<nodelist->size(); ++i){
        p = (*nodelist)[i] ;
        if(p->data.level == level){
            if (!p->data.isLeaf) {
                l = p->leftChild;
                r = p->rghtChild;
                l->data.address = p->data.address + 1;
                r->data.address = p->data.address + 1 + l->data.size ;
            }
        }
    }
    return ;
}

int kdTree::reduce(unsigned char *list, long int size){
    int sum =0;
    for(long int i=0; i<size; ++i){
        sum += list[i] ;
    }
    return sum;
}

int kdTree::reduce(std::vector<int> list){
    int sum =0;
    for(int i=0; i<list.size(); ++i){
        sum += list[i] ;
    }
    return sum;
}

// Kernels for reductions
//

void kdTree::scanInclusive(int *in, int *out, int n){
    out[0] = in[0];
    for (int k=1; k<n; ++k){
        out[k] = in[k] + out[k-1] ;
    }
}

void kdTree::scanExclusive(int *in, int *out, int n){
    out[0] = 0;
    for (int k=1; k<n; ++k){
        out[k] = in[k-1] + out[k-1] ;
    }
}

//float kdTree::reduceSum(float *arr, int nProcs)
//{
//    int a1=0,a2 ;
//    for(int p=0; p<nProcs; ++p){
//        for(int i=0; i<log2(nProcs)-1; i++){
//            a1 = pow(2,i+1)*p;
//            a2 = a1 + pow(2,i);
//            if(a2<nProcs){
//                arr[a1] = arr[a1] + arr[a2] ;
//            }
//        }
//    }
//    return arr[a1] ;
//}
//
//int kdTree::reduceSum(int *arr, int nProcs)
//{
//    int a1 = 0,a2 ;
//    for(int p=0; p<nProcs; ++p){
//        for(int i=0; i<log2(nProcs)-1; i++){
//            a1 = pow(2,i+1)*p;
//            a2 = a1 + pow(2,i);
//            if(a2<nProcs){
//                arr[a1] = arr[a1] + arr[a2] ;
//            }
//        }
//    }
//    return arr[a1] ;
//}
//
//float kdTree::reduceMin(float *arr, int nProcs)
//{
//    int a1 = 0,a2 ;
//    for(int p=0; p<nProcs; ++p){
//        for(int i=0; i<log2(nProcs)-1; i++){
//            a1 = pow(2,i+1)*p;
//            a2 = a1 + pow(2,i);
//            if(a2<nProcs){
//                arr[a1] = (arr[a1] < arr[a2]) ? arr[a1] : arr[a2] ;
//            }
//        }
//    }
//    return arr[a1] ;
//}
//
//float kdTree::reduceMax(float *arr, int nProcs)
//{
//    int a1= 0,a2 ;
//    for(int p=0; p<nProcs; ++p){
//        for(int i=0; i<log2(nProcs)-1; i++){
//            a1 = pow(2,i+1)*p;
//            a2 = a1 + pow(2,i);
//            if(a2<nProcs){
//                arr[a1] = (arr[a1] > arr[a2]) ? arr[a1] : arr[a2] ;
//            }
//        }
//    }
//    return arr[a1] ;
//}
//void kdTree::segReduceMin(float *data, int *owner, int nElements, float **results)
//{
//    int p0, p1, w0,w1 ;
//    for (int d=0; d<log2(nElements)-1; ++d){
//        p0 = pow(2,d) ;
//        p1 = pow(2,d+1) ;
//        for (int i=0; i < ((nElements-1) / p1); ++i) {
//            w0 = owner[p1*i] ;
//            w1 = owner[p1*i+p0] ;
//            if (w0 != w1){
//                (*results)[w1] = ((*results)[w1] < data[p1*i + p0]) ? (*results)[w1] : data[p1*i + p0] ;
//            }else{
//                data[p1*i] = (data[p1*i] < data[p1*i + p0]) ? data[p1*i] : data[p1*i + p0] ;
//            }
//        }
//    }
//    return ;
//}
//
//void kdTree::segReduceMax(float *data, int *owner, int nElements, float **results)
//{
//    int p0, p1, w0,w1 ;
//    for (int d=0; d<log2(nElements)-1; ++d){
//        p0 = pow(2,d) ;
//        p1 = pow(2,d+1) ;
//        for (int i=0; i < ((nElements-1) / p1); ++i) {
//            w0 = owner[p1*i] ;
//            w1 = owner[p1*i+p0] ;
//            if (w0 != w1){
//                (*results)[w1] = ((*results)[w1] > data[p1*i + p0]) ? (*results)[w1] : data[p1*i + p0] ;
//            }else{
//                data[p1*i] = (data[p1*i] > data[p1*i + p0]) ? data[p1*i] : data[p1*i + p0] ;
//            }
//        }
//    }
//    return ;
//}

void kdTree::writeAABBtoPlyFile(AABB bv, std::string filename){
    
    int verticesInAABB = 8;
    int facesInAABB = 6;
    
    FILE *fp = fopen(filename.c_str(), "w");
    
    if (fp == NULL) {
        printf("Error : could not open file %s for writing\n",filename.c_str());
        exit(1);
    }
    
    fprintf(fp,"ply\n");
    fprintf(fp,"format ascii 1.0\n");
    fprintf(fp,"comment AABB PLY File by fastKDTree\n");
    fprintf(fp,"element vertex %d\n",verticesInAABB);
    fprintf(fp,"property float x\n");
    fprintf(fp,"property float y\n");
    fprintf(fp,"property float z\n");
    fprintf(fp,"element face %d\n",facesInAABB);
    fprintf(fp,"property list uchar int vertex_index\n");
    fprintf(fp,"end_header\n");
    
    // vertex info for AABB
    //
    fprintf(fp,"%4.4f %4.4f %4.4f\n",bv.AA.x,bv.AA.y,bv.AA.z);  // Vertex 0
    fprintf(fp,"%4.4f %4.4f %4.4f\n",bv.BB.x,bv.AA.y,bv.AA.z);  // Vertex 1
    fprintf(fp,"%4.4f %4.4f %4.4f\n",bv.BB.x,bv.BB.y,bv.AA.z);  // Vertex 2
    fprintf(fp,"%4.4f %4.4f %4.4f\n",bv.AA.x,bv.BB.y,bv.AA.z);  // Vertex 3
    fprintf(fp,"%4.4f %4.4f %4.4f\n",bv.AA.x,bv.AA.y,bv.BB.z);  // Vertex 4
    fprintf(fp,"%4.4f %4.4f %4.4f\n",bv.BB.x,bv.AA.y,bv.BB.z);  // Vertex 5
    fprintf(fp,"%4.4f %4.4f %4.4f\n",bv.BB.x,bv.BB.y,bv.BB.z);  // Vertex 6
    fprintf(fp,"%4.4f %4.4f %4.4f\n",bv.AA.x,bv.BB.y,bv.BB.z);  // Vertex 7
    
    // face info for AABB
    //
    fprintf(fp, "4 0 1 2 3\n");
    fprintf(fp, "4 4 5 6 7\n");
    fprintf(fp, "4 0 1 5 4\n");
    fprintf(fp, "4 1 2 6 5\n");
    fprintf(fp, "4 2 3 7 6\n");
    fprintf(fp, "4 3 0 4 7\n");
    
    fclose(fp) ;
    return ;
    
}

void kdTree::dumpsmall(treeList **list){
    kdTreeNode *p, *n;
    
    if ((*list)->size() == 0) {
        return;
    }
    
    p = (*list)->front();
    int i=0;
    while(p){
        n = p;
        printf("[%02d]<%p>  triMask=<%p>[%02d] ",i,p,p->data.triangleMask,p->data.smallntris) ;
        for(int j=0; j<p->data.smallntris; ++j){
            printf("%d",p->data.triangleMask[j]) ;
        }
        printf("\n");
        p=n->next;
        ++i;
    }
}

void kdTree::swapLists(treeList **a, treeList **b) {
    treeList *t;
    t = *a ;
    *a = *b ;
    *b = t ;
    return  ;
}

void kdTree::printKdTreeNodes(std::vector<kdTreeNode *> nodelist){
    
    kdTreeNode *p;
    
    for (int i=0; i<nodelist.size(); ++i){
        p = nodelist[i] ;
        printf("[<%p>]",p);
        for(int j=0;j<p->data.level;++j){
            printf("-");
        }
        if(p->data.isLeaf){
            printf("X");
            printf(" #%02d = [",p->data.smallntris);
            for(int k=0; k<p->data.smallntris; ++k){
                if (p->data.triangleMask[k] == 1) {
                    printf(" %02d", p->smallroot->data.triangles[k]);
                }
            }
            printf(" ]");
        }else{
            printf("|");
            printf("  [%p][%p]",p->leftChild ,p->rghtChild);
        }
        printf("\n");
    }
    
    return;
}

void kdTree::printKdTreeData(KdData **kdTree, int *numOfKdTreeNodes){
    
    KdData * node;
    
    printf("\n\n");
    printf("                        KdTree Data Output (post Pre-Traversal Correction)\n");
    printf("| --------------------------------------------------------------------------------------------------|\n");
    printf("| --------------------------------------------------------------------------------------------------|\n");
    printf("|                                       Key                                                         |\n");
    printf("| --------------------------------------------------------------------------------------------------|\n");
    printf("| [IDX]         : Index/address of each node in the tree                                            |\n");
    printf("| [A]           : If a branch then dimension of split plane. If Leaf then **LEAF**                  |\n");
    printf("| [BBBB]        : If branch the the position ofthe split plane                                      |\n");
    printf("| [CCCCCCCCCCC] : If branch then the index/address of left and right children.                      |\n");
    printf("|               : If Leaf then the number of triangles in the leaf followed by the triangle indices |\n");
    printf(" -------------------------------------------------------------------------------------------------- |\n");
    printf("[IDX] [A][BBBB]  = [CCCCCCCCCCC]\n");
    printf("--------------------------------\n");
    
    for(int i=0; i<*numOfKdTreeNodes; ++i){
        
        node = &((*kdTree)[i]) ;
        if(KDT_ISLEAF( node ) ){
            if(KDT_NUMTRIS(node) > 0){
                printf("[%03d]  ",KDT_INDEX(node));
                printf("**LEAF**  = # %d ",KDT_NUMTRIS(node));
                printf("[");
                for(int j=0; j<KDT_NUMTRIS(node);++j){
                    printf(" %d",(*kdTree)[i+j+1].leaf.triangleIndex);
                }
                printf(" ]");
                printf("\n");
            }
            
        }else{
            printf("[%03d]  ",i);
            switch (KDT_DIMENSION(node)) {
                case 0:
                    printf("X ");
                    break;
                case 1:
                    printf("Y ");
                    break;
                case 2:
                    printf("Z ");
                    break;
                default:
                    printf("Unusual dimension ! %d\n ", (unsigned int)KDT_DIMENSION(node)) ;
                    break;
            }
            printf(" %-6.2f = ",KDT_SPLITPOS(node)) ;
            printf("LEFT= [%02d], RGHT=[%02d]\n",KDT_LEFTCHILD(node),KDT_RGHTCHILD(node));
            
        }
        
    }
    return ;
}

void kdTree::printSummary(std::vector<kdTreeNode *> nodelist)
// Note that this complexity analysis is taken from [1]
//
//  1. Wald, Ingo, and Vlastimil Havran. "On building fast kd-trees for ray tracing,
//     and on doing that in O (N log N)." Interactive Ray Tracing 2006, IEEE
//     Symposium on. IEEE, 2006.
//
{
    int Nl, Nne,ntri,nt;
    double  Nat, Et, El, Ei, Ct, SAs, x, Kt, Ki;
    
    Kt = TRAVERSALCOST ;
    Ki = INTERSECTIONCOST ;

    SAs = nodelist[0]->data.aabb.surfaceArea() ;
    Nl = Nne = ntri = Et = El = Ei = 0;
    for (int i=0; i<nodelist.size(); ++i) {
        Et += nodelist[i]->data.aabb.surfaceArea() / SAs ;
        if (nodelist[i]->data.isLeaf) {
            x = nodelist[i]->data.aabb.surfaceArea() / SAs ;
            El += x ;
            Nl++;
            nt = reduce(nodelist[i]->data.triangleMask, nodelist[i]->data.smallntris);
            if (nt != 0 ) {
                Ei += (x * nt);
                Nne++;
                ntri += nt;
            }
        }
    }
    
    Nat = (double)ntri / (double)Nne;
    Ct  = (Et * Kt) + (Ei * Ki) ;
    printf("\n            KDTree Summary\n");
    printf("==========================================\n");
    printf("  Triangles                       : %4ld\n",globalMesh.triangles.size());
    printf("  Number of Nodes                 : %4ld\n",nodelist.size());
    printf("  Number of leaves (NL)           : %4d \n",Nl);
    printf("  Number of non-empty leaves      : %4d\n",Nne);
    printf("  Ave triangles / non-empty leaf  : %4.1f\n",Nat);
    printf("  Expected Traversals             : %4.1f\n",Et);
    printf("  Expected leaves visited         : %4.1f\n",El);
    printf("  Expected triangle intersections : %4.1f\n",Ei);
    printf("------------------------------------------\n");
    printf("  Total tree traversal cost       : %4.1f\n\n",Ct);
    
    return ;
    
}
