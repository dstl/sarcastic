/***************************************************************************
 * 
 *           Module :  kdTreeNode.cpp
 *          Program :  fastKdTree
 *       Created by :  Darren Muff on 05/03/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  04-Nov-2018
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
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
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
 ***************************************************************************/

#include "kdTreeNode.hpp"


kdTree::kdTreeNode::kdTreeNode(){
    // default constructor
    data.size = 0 ;
    data.address = 0 ;
    data.triangleIndex = 0 ;
    data.isLeaf = false ;
    data.dim = -1 ;
    data.level = 0 ;
    data.smallntris = 0 ;
    data.triangleMask = NULL;
    next = NULL ;
    leftChild = NULL ;
    rghtChild = NULL ;
    smallroot = NULL ;
}

kdTree::kdTreeNode::kdTreeNode(const kdTreeNode *node){
    data.size = node->data.size;
    data.address = node->data.address ;
    data.triangleIndex = node->data.triangleIndex;
    data.aabb = node->data.aabb ;
    data.triAABBs = node->data.triAABBs ;
    data.splitList = node->data.splitList ;
    data.triangles = node->data.triangles ;
    data.isLeaf = node->data.isLeaf ;
    data.dim = node->data.dim;
    data.splitPos = node->data.splitPos ;
    data.level = node->data.level ;
    data.smallntris = node->data.smallntris ;
    if (node->data.smallntris > 0) {
        data.triangleMask = new unsigned char [node->data.smallntris] ;
        memcpy(data.triangleMask, node->data.triangleMask, node->data.smallntris * sizeof(unsigned char)) ;
    }else{
        data.triangleMask = NULL;
    }
    
    next = node->next ;
    leftChild = node->leftChild ;
    rghtChild = node->rghtChild ;
    smallroot = node->smallroot ;
    
}

kdTree::kdTreeNode::kdTreeNode(std::vector<int> tris) {
    data.triangles = tris ;
    data.aabb = BVforAllTris() ;
}

// Initialise the root KdTreeNode using a .plyFile
//
kdTree::kdTreeNode::kdTreeNode(std::string plyFileName, TriangleMesh &globalMesh)
{
    globalMesh.readPLYFile(plyFileName);
    globalMesh.checkIntegrityAndRepair();
    globalMesh.buildTriangleAABBs();
    
    data.triAABBs = globalMesh.AABBs ;
    for(int i=0; i< globalMesh.triangles.size(); ++i){
        data.triangles.push_back(i) ;
    }
    data.aabb = BVforAllTris() ;
}

kdTree::kdTreeNode::kdTreeNode(TriangleMesh &mesh)
{
    if (mesh.AABBs.size()==0) {
        mesh.buildTriangleAABBs();
    }
    data.triAABBs = mesh.AABBs ;
    for(int i=0; i< mesh.triangles.size(); ++i){
        data.triangles.push_back(i) ;
    }
    data.aabb = BVforAllTris() ;
    
}


//kdTree::kdTreeNode::~kdTreeNode(){
//    
//    if (data.smallntris > 0) {
//        delete [] data.triangleMask ;
//    }
//    
//}

// Some routines for handling the kdTreeNode as a list
//


AABB kdTree::kdTreeNode::BVforAllTris()
{
    SPVector triaa, tribb;
    AABB ans;
    
    VECT_CREATE( 9e9,  9e9,  9e9, ans.AA);
    VECT_CREATE(-9e9, -9e9, -9e9, ans.BB);
    
    for(int i=0; i<data.triangles.size(); ++i){
        triaa = data.triAABBs[i].AA ;
        tribb = data.triAABBs[i].BB ;
        
        ans.AA.x = (triaa.x < ans.AA.x) ? triaa.x : ans.AA.x ;
        ans.AA.y = (triaa.y < ans.AA.y) ? triaa.y : ans.AA.y ;
        ans.AA.z = (triaa.z < ans.AA.z) ? triaa.z : ans.AA.z ;
        
        ans.BB.x = (tribb.x > ans.BB.x) ? tribb.x : ans.BB.x ;
        ans.BB.y = (tribb.y > ans.BB.y) ? tribb.y : ans.BB.y ;
        ans.BB.z = (tribb.z > ans.BB.z) ? tribb.z : ans.BB.z ;
    }
    return ans;
}

void kdTree::kdTreeNode::split(int k, double pos, std::shared_ptr<kdTreeNode> left, std::shared_ptr<kdTreeNode> rght)
{
    if (smallroot == NULL) {
        printf("Error: split() function only works on small nodes. Did you mean medianSplit()?\n");
        exit(1) ;
    }
    
    data.splitPos = pos;
    data.dim = k;
    
    // Set the AABBs for the kids
    //
    (left)->data.aabb.AA = data.aabb.AA ;
    (left)->data.aabb.BB = data.aabb.BB ;
    (left)->data.aabb.BB.cell[k] = pos ;
    for(int i=0;i<3; ++i){
        if ((left)->data.aabb.BB.cell[i] < (left)->data.aabb.AA.cell[i]) {
            printf("ERROR: Inside out AABB for dimension %d\n",i);
            printf("left node AABB (dimension %d) : %f - %f\n",i,(left)->data.aabb.AA.cell[i],(left)->data.aabb.BB.cell[i]);
            printf("\n");
        }
    }
    (rght)->data.aabb.AA = data.aabb.AA ;
    (rght)->data.aabb.BB = data.aabb.BB ;
    (rght)->data.aabb.AA.cell[k] = pos ;
    for(int i=0;i<3; ++i){
        if ((rght)->data.aabb.BB.cell[i] < (rght)->data.aabb.AA.cell[i]) {
            printf("ERROR: Inside out AABB for dimension %d\n",i);
            printf("rght node AABB (dimension %d) : %f - %f\n",i,(rght)->data.aabb.AA.cell[i],(rght)->data.aabb.BB.cell[i]);
            printf("\n");
        }
    }
    
    // Sort triangles from this node into child nodes using bit mask
    //
    (left)->data.triangleMask = new unsigned char [data.smallntris] ;
    (rght)->data.triangleMask = new unsigned char [data.smallntris] ;
    
    for(int t=0; t<data.smallntris; ++t){
        (left)->data.triangleMask[t] = data.triangleMask[t] & ((data.triAABBs[t].AA.cell[k] <= pos)) ;
        (rght)->data.triangleMask[t] = data.triangleMask[t] & ((data.triAABBs[t].BB.cell[k] >= pos)) ;
    }
    
    (left)->data.splitList  = data.splitList  ;
    (left)->data.triAABBs   = data.triAABBs   ;
    (left)->smallroot       = smallroot       ;
    (left)->data.smallntris = data.smallntris ;
    (left)->data.dim        = data.dim        ;
    (left)->data.splitPos   = data.splitPos   ;
    (rght)->data.splitList  = data.splitList  ;
    (rght)->data.triAABBs   = data.triAABBs   ;
    (rght)->smallroot       = smallroot       ;
    (rght)->data.smallntris = data.smallntris ;
    (rght)->data.dim        = data.dim        ;
    (rght)->data.splitPos   = data.splitPos   ;

    // Set level for the children
    //
    (left)->data.level = data.level+1;
    (rght)->data.level = data.level+1;
    
    return ;
}

void kdTree::kdTreeNode::medianSplit(std::shared_ptr<kdTreeNode> *left, std::shared_ptr<kdTreeNode> *rght, TriangleMesh &globalMesh)
// medianSplit() splits the node at the median point of the largest axis
// and returns two child nodes that have their AABB's set and the triAABB's correctly clipped to the split plane.
// In addition the triangles array in each child correctly holds the index of each triangle
//
{
    
    // Set level for the children
    //
    (*left)->data.level = data.level+1;
    (*rght)->data.level = data.level+1;
    
    // find longest axis
    //
    double len[3], pos;
    int maxAxis ;
    len[0] = data.aabb.BB.x - data.aabb.AA.x ;
    len[1] = data.aabb.BB.y - data.aabb.AA.y ;
    len[2] = data.aabb.BB.z - data.aabb.AA.z ;
    if(len[0] > len[1]){
        if(len[0] > len[2]){ maxAxis = 0;}else{ maxAxis = 2;}
    }else if (len[1] > len[2]){
        maxAxis = 1;
    }else{
        maxAxis = 2;
    }
    pos = data.aabb.AA.cell[maxAxis]+(len[maxAxis] /2.0) ;
    
    data.splitPos = pos ;
    data.dim      = maxAxis;
    
    // Set the AABBs for the kids by splitting this node
    // down the middle
    //
    (*left)->data.aabb.AA = data.aabb.AA ;
    (*left)->data.aabb.BB = data.aabb.BB ;
    (*left)->data.aabb.BB.cell[maxAxis] = pos ;
    for(int i=0;i<3; ++i){
        if ((*left)->data.aabb.BB.cell[i] < (*left)->data.aabb.AA.cell[i]) {
            printf("ERROR: Inside out AABB for dimension %d\n",i);
            printf("left node AABB (dimension %d) : %f - %f\n",i,(*left)->data.aabb.AA.cell[i],(*left)->data.aabb.BB.cell[i]);
            printf("\n");
        }
    }
    (*rght)->data.aabb.AA = data.aabb.AA;
    (*rght)->data.aabb.BB = data.aabb.BB;
    (*rght)->data.aabb.AA.cell[maxAxis] = pos ;
    for(int i=0;i<3; ++i){
        if ((*rght)->data.aabb.BB.cell[i] < (*rght)->data.aabb.AA.cell[i]) {
            printf("ERROR: Inside out AABB for dimension %d\n",i);
            printf("rght node AABB (dimension %d) : %f - %f\n",i,(*rght)->data.aabb.AA.cell[i],(*rght)->data.aabb.BB.cell[i]);
            printf("\n");
        }
    }
    
    // sort the triangles in the mesh in this node into the mesh for each of the child nodes.
    //
    
    SPVector vertA,vertB,vertC;
    AABB ablft, abrgt;
    int triIdx ;
    
    for (int t=0; t<data.triangles.size(); ++t){
        triIdx = data.triangles[t] ;
        if(data.triAABBs[t].AA.cell[maxAxis] >= pos ){
            (*rght)->data.triangles.push_back(triIdx);
            (*rght)->data.triAABBs.push_back(data.triAABBs[t]) ;
        }
        if(data.triAABBs[t].BB.cell[maxAxis] <= pos){
            (*left)->data.triangles.push_back(triIdx);
            (*left)->data.triAABBs.push_back(data.triAABBs[t]) ;
        }
        if( (data.triAABBs[t].AA.cell[maxAxis] < pos) && (data.triAABBs[t].BB.cell[maxAxis] > pos) )
        {
            // Triangle straddles the split
            //
            vertA = globalMesh.vertices[globalMesh.triangles[triIdx].a].asSPVector() ;
            vertB = globalMesh.vertices[globalMesh.triangles[triIdx].b].asSPVector() ;
            vertC = globalMesh.vertices[globalMesh.triangles[triIdx].c].asSPVector() ;
            
            data.triAABBs[t].clipToTriangle(vertA, vertB, vertC, pos, maxAxis, ablft, abrgt);
            
            if(ablft.surfaceArea() > 1.0e-10 ){
                (*left)->data.triangles.push_back(triIdx);
                (*left)->data.triAABBs.push_back(ablft) ;
            }
            if(abrgt.surfaceArea() > 1.0e-10 ){
                (*rght)->data.triangles.push_back(triIdx);
                (*rght)->data.triAABBs.push_back(abrgt) ;
            }
    
        }
    }
    
    return;
}

// Member functions related to handling the linked list follow
//
kdTree::treeList::treeList(){
    head = NULL ;
    length = 0 ;
}


int kdTree::treeList::size(){ return length; };

bool kdTree::treeList::indexCheck(int position){
    if ( (position < 0) || (position >= length) ){
        std::cout << "\nIndex error : position " << position << " out of range ( list size:"<<size()<<" )\n" ;
        return false ;
    }
    return	true;
}

bool kdTree::treeList::insertNode(std::shared_ptr<kdTreeNode> node, int position) {
    if( !indexCheck(position) ) {
        exit(1);
    }
    if (head->next == NULL) {
        head->next = node;
        length++;
        return true;
    }
    int count = 0;
    std::shared_ptr<kdTreeNode> p = head;
    std::shared_ptr<kdTreeNode> q = head;
    while (q) {
        if (count==position) {
            p->next = node;
            node->next = q;
            length++;
            return true;
        }
        p = q;
        q = p->next;
        count++;
    }
    std::cout << "\nError: node was not added to list.\n";
    return false;
} ;

bool kdTree::treeList::removeNode(int position) {
    if( !indexCheck(position) ) {
        exit(1);
    }
    if (head->next == NULL) {
        std::cout << "\nError: there is nothing to remove.\n";
        return false;
    }
    int count = 0;
    std::shared_ptr<kdTreeNode> p = head;
    std::shared_ptr<kdTreeNode> q = head;
    while (q) {
        if (count == position) {
            p->next = q->next;
            q.reset();
            length--;
            return true;
        }
        p = q;
        q = p->next;
        count++;
    }
    std::cout << "n\Error: nothing was removed from the list.\n";
    return false;
}

std::shared_ptr<kdTree::kdTreeNode> kdTree::treeList::at(int position){
    if( !indexCheck(position) ) {
        exit(1);
    }
    if (position == 0) {
        return head ;
    }
    
    int count = 1;
    std::shared_ptr<kdTreeNode> p;
    std::shared_ptr<kdTreeNode> q = head;
    while (q) {
        if (count == position) {
            if (q->next == NULL) {
                std::cout << "\nError : Null pointer at location "<<count<<"\n" ;
                exit (1);
            }
            return q->next ;
        }
        p = q ;
        q = p->next ;
        count++;
    }
    
    std::cout << " Index not found in list\n" ;
    return NULL;
}

std::shared_ptr<kdTree::kdTreeNode> kdTree::treeList::front(){
    return head;
}

std::shared_ptr<kdTree::kdTreeNode> kdTree::treeList::back(){
    return tail;
}

void kdTree::treeList::push_back(std::shared_ptr<kdTreeNode> node){
    
    // Create a new node pointer with allocated memory
    // so that when this function scope ends the memory
    // allocated remains.
    //
    
    if(node->next != NULL){
        printf("warning : adding list to list (not just a node)\n");
    }
    
    if (head == NULL) {
        head = node ;
        tail = node ;
        length++;
        return ;
    }
    
    std::shared_ptr<kdTreeNode> p = tail ;
    std::shared_ptr<kdTreeNode> q = tail ;
    while (q){
        p = q ;
        q = p->next ;
    }
    p->next = node ;
    tail = node ;
    length++ ;
    return ;
}

void kdTree::treeList::clear() {
    std::shared_ptr<kdTreeNode> p = head;
    while (p) {
        head = tail->next;
        p = head ;
    }
    length = 0;
    return ;
}

void kdTree::treeList::erase() {
    
    std::shared_ptr<kdTreeNode> p = head;
    while (p) {
        head = head->next;
        p.reset();
        p = head;
    }
    length = 0;
    return ;
}

void kdTree::treeList::printList() {
    std::shared_ptr<kdTreeNode> p ;
    std::shared_ptr<kdTreeNode> q = head;
    int count = 0 ;
    std::cout << "\n---------------------------\n";
    std::cout << "List output ( length "<<length<<" )"<<std::endl;
    while (q)
    {
        p = q;
        printf("[%2d] <%p> left:<%p> right:<%p> smallroot:<%p>",count, p.get(), p->leftChild.get(),p->rghtChild.get(),p->smallroot.get());
        if(p->next == NULL){
            printf("next is NULL\n");
        }else{
            printf("next is %p\n", p->next.get()) ;
        }
        q = p -> next;
        count++;
    }
    return ;
}

bool kdTree::treeList::empty(){
    if (length == 0) return true;
    return false;
}


kdTree::treeList::~treeList(){
    clear();
    return ;
}



