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
std::vector<kdTreeNode *> nodelist ;

kdTreeNode::kdTreeNode(){
    // default constructor
    // this is to allow us to create an object without any initialization
}

kdTreeNode::kdTreeNode(const kdTreeNode *node){
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

// Copy constructor
//
//kdTreeNode::kdTreeNode(const kdTreeNode &node){
//    
//    
//    data = node.data ;
//    next = node.next ;
//    leftChild = node.leftChild ;
//    rghtChild = node.rghtChild ;
//    smallroot = node.smallroot ;
//    
//    if (node.data.smallntris > 0) {
//        data.triangleMask = new unsigned char [node.data.smallntris] ;
//        memcpy(data.triangleMask, node.data.triangleMask, node.data.smallntris * sizeof(unsigned char)) ;
//    }
//}


kdTreeNode::kdTreeNode(std::vector<int> tris) {
    data.triangles = tris ;
    data.aabb = BVforAllTris() ;
}

// Initialise the root KdTreeNode using a .plyFile
//
kdTreeNode::kdTreeNode(std::string plyFileName)
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


kdTreeNode::~kdTreeNode(){
    
    if (data.smallntris > 0) {
        delete [] data.triangleMask ;
    }
    
}

// Some routines for handling the kdTreeNode as a list
//


AABB kdTreeNode::BVforAllTris()
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

void kdTreeNode::split(int k, float pos, kdTreeNode *left, kdTreeNode *rght)
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
    (rght)->data.aabb.AA = data.aabb.AA ;
    (rght)->data.aabb.BB = data.aabb.BB ;
    (rght)->data.aabb.AA.cell[k] = pos ;
    
    // Sort triangles from this node into child nodes using bit mask
    //
    (left)->data.triangleMask = new unsigned char [data.smallntris] ;
    (rght)->data.triangleMask = new unsigned char [data.smallntris] ;
    
    for(int t=0; t<data.smallntris; ++t){
        (left)->data.triangleMask[t] = data.triangleMask[t] & (data.triAABBs[t].AA.cell[k] <= pos) ;
        (rght)->data.triangleMask[t] = data.triangleMask[t] & (data.triAABBs[t].BB.cell[k] >= pos) ;
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

void kdTreeNode::medianSplit(kdTreeNode **left, kdTreeNode **rght)
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
    
    (*rght)->data.aabb.AA = data.aabb.AA;
    (*rght)->data.aabb.BB = data.aabb.BB;
    (*rght)->data.aabb.AA.cell[maxAxis] = pos ;
    
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
            
            if(ablft.surfaceArea() > 0.0001 ){
                (*left)->data.triangles.push_back(triIdx);
                (*left)->data.triAABBs.push_back(ablft) ;
            }
            if(abrgt.surfaceArea() > 0.0001 ){
                (*rght)->data.triangles.push_back(triIdx);
                (*rght)->data.triAABBs.push_back(abrgt) ;
            }
    
        }
    }
    
    return;
}

void printKdTreeNodes(std::vector<kdTreeNode *> nodelist){
    
    kdTreeNode *p,*q,*n;

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

void printKdTreeData(){
    
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

    for(int i=0; i<numOfKdTreeNodes; ++i){
        
        node = &(kdTreeOutput[i]) ;
        if(KDT_ISLEAF( node ) ){
            if(KDT_NUMTRIS(node) > 0){
                printf("[%03d]  ",i);
                printf("**LEAF**  = # %d ",KDT_NUMTRIS(node));
                printf("[");
                for(int j=0; j<KDT_NUMTRIS(node);++j){
                    printf(" %d",kdTreeOutput[i+j+1].leaf.triangleIndex);
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

// Member functions related to handling the linked list follow
//
treeList::treeList(){
    head = NULL ;
    length = 0 ;
}


int treeList::size(){ return length; };

bool treeList::indexCheck(int position){
    if ( (position < 0) || (position >= length) ){
        std::cout << "\nIndex error : position " << position << " out of range ( list size:"<<size()<<" )\n" ;
        return false ;
    }
    return	true;
}

bool treeList::insertNode(kdTreeNode *node, int position) {
    if( !indexCheck(position) ) {
        exit(1);
    }
    if (head->next == NULL) {
        head->next = node;
        length++;
        return true;
    }
    int count = 0;
    kdTreeNode * p = head;
    kdTreeNode * q = head;
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

bool treeList::removeNode(int position) {
    if( !indexCheck(position) ) {
        exit(1);
    }
    if (head->next == NULL) {
        std::cout << "\nError: there is nothing to remove.\n";
        return false;
    }
    int count = 0;
    kdTreeNode * p = head;
    kdTreeNode * q = head;
    while (q) {
        if (count == position) {
            p->next = q->next;
            delete q;
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

kdTreeNode * treeList::at(int position){
    if( !indexCheck(position) ) {
        exit(1);
    }
    if (position == 0) {
        return head ;
    }
    
    int count = 1;
    kdTreeNode * p;
    kdTreeNode * q = head;
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

kdTreeNode * treeList::front(){
    return head;
}

kdTreeNode * treeList::back(){
    return tail;
}

void treeList::push_back(kdTreeNode *node){
    
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
    
    kdTreeNode * p = tail ;
    kdTreeNode * q = tail ;
    while (q){
        p = q ;
        q = p->next ;
    }
    p->next = node ;
    tail = node ;
    length++ ;
    return ;
}

void treeList::clear() {
    kdTreeNode *p = head;
    while (p) {
        head = tail->next;
        p = head ;
    }
    length = 0;
    return ;
}

void treeList::erase() {
    kdTreeNode * p = head;
    while (p) {
        head = head->next;
        delete p;
        p = head;
    }
    length = 0;
    return ;
}

void swapLists(treeList **a, treeList **b) {
    treeList *t;
    t = *a ;
    *a = *b ;
    *b = t ;
    return  ;
}

void treeList::printList() {
    kdTreeNode * p ;
    kdTreeNode * q = head;
    int count = 0 ;
    std::cout << "\n---------------------------\n";
    std::cout << "List output ( length "<<length<<" )"<<std::endl;
    while (q)
    {
        p = q;
        printf("[%2d] <%p> left:<%p> right:<%p> smallroot:<%p>",count, p, p->leftChild,p->rghtChild,p->smallroot);
        if(p->next == NULL){
            printf("next is NULL\n");
        }else{
            printf("next is %p\n", p->next) ;
        }
        q = p -> next;
        count++;
    }
    return ;
}

bool treeList::empty(){
    if (length == 0) return true;
    return false;
}


treeList::~treeList(){
    clear();
    return ;
}



