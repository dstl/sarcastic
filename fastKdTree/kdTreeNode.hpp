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
extern int numOfTriangleIndices ;


#define KDT_ISLEAF(n)    ((n->leaf.leafDim & (unsigned char)(0x4))>>2)
#define KDT_DIMENSION(n) (n->leaf.leafDim & 0x3)
#define KDT_NUMTRIS(n)   (n->leaf.ntriangles)
//#define KDT_OFFSET(n)    ((n->leaf.leafDim & (0x7FFFFFFC))>>2)
#define KDT_LEFTCHILD(n)    (n->brch.leftaddress)
#define KDT_RGHTCHILD(n)    (n->brch.rghtaddress)
#define KDT_SPLITPOS(n)     (n->leaf.splitPosition)
#define KDT_TRIINDEX(n)     (n->leaf.triangleIndex)
#define TRAVERSALCOST ((float)(15.0))
#define INTERSECTIONCOST ((float)(20.0))
#define SMALLSIZE (72)
#define Ce (0.25)
#define MAXTREELEVEL 10

//typedef union {
//    struct KdTreeLeaf {
//        unsigned int flagDimAndOffset;
//        // bits 0..1        : splitting dimension
//        // bits 2..30       : offset bits to 'leaf list'
//        // bits 31 (sign)   : flag whether node is a leaf
//        float splitPosition;
//    } leaf;
//    
//    struct KdTreeBranch {
//        unsigned int flagDimAndOffset;
//        // bits 0..30       : offset to first child
//        // bit 31 (sign)    : flag whether node is a leaf
//        float splitPosition;
//    } branch;
//} KdData ;

//  leafDim (unsigned char)
//
//           _______third LSB = flag for leaf 1 = leaf
//           | ____ 2 LSB = dimension 1=x,2=y or 3=z
//           | | |
//           v v v
// x x x x x 1 1 1
//
typedef union {
    struct KdTreeLeaf {
        float           splitPosition ; // splitposition for node
        unsigned char   leafDim ;       // lsb 2 = leaf flag. Bits 0 + 1 are for dimension (above)
        unsigned int    ntriangles ;    // num of triangles in this leaf
        unsigned int    triangleIndex ; // used to store the triangle index info
    } leaf;
    struct KdTreeBrch {
        float           splitPosition ; // splitposition for node
        unsigned char   leafDim ;       // lsb 2 = leaf flag. Bits 0 + 1 are for dimension (above)
        unsigned int    leftaddress ;   // index of left child in structure
        unsigned int    rghtaddress ;   // index of right child in structure
    } brch;
    
} KdData ;


extern KdData * kdTreeOutput ;
extern int numOfKdTreeNodes ;
//
// kdTreeNode Notes
// ----------------
// The kdTreeNode contains all the information for a given node in a kdTree.
// Each node can be a leaf (which contains triangles) or a branch. isLeaf is true for a leaf node otherwise its a branch
// It references triangles stored in a global mesh defined in 'globalMesh' using the index of the triangle in the
// mesh rather than the triangle data itself. (for efficiency)
// Each node can be part of a linked list with the next node referenced by nextPtr
// If the node is not a leaf (a branch) then it has two children pointed to by leftPtr and rghtPtr.
// The left child is on the numerically lower side of the splitplane defined by dim (the dimension x,y,z) and splitPos.
// The right child is the higher side of the split position splitPos
// When a node contains only a few triangles it is defined as 'small' As there are many 'small' nodes it is more efficient to store
//      the triangle information as a mask of the triangles (triangleMask, size=smallntris) in the first of its ancestors that was
//      defined as 'small'. This ancestor is pointed to by smallroot.
//
class kdTreeNode {
private:
public:
    struct data {
        int                         size            = 0 ;       // Number of triangles in this node (used to calc memory reqs)
        int                         address         = 0 ;       // Array index of this node in output
        int                         triangleIndex   = 0 ;       // index of first triangle
        AABB                        aabb ;                      // Axis aligned bounding box for this node
        std::vector<AABB>           triAABBs ;                  // vector array of axis-aligned bounding boxes-one for each triangle in this node
        std::vector<splitCandidate> splitList ;                 // vector array of split candidates. One for each dimension and for each triangle vertex in each dim
        std::vector<int>            triangles ;                 // index of triangle in global Mesh
        bool                        isLeaf          = false ;   // Is this node a leaf (containing triangles) or a branch ?
        int                         dim             = -1 ;      // Dimension of split axis, 0,1,2 = x,y,z
        float                       splitPos ;                  // position of splitplane for this node along axis defined by dim
        int                         level           = 0 ;       // depth down the tree for this node (for denug purposes)
        int                         smallntris      = 0 ;       // Number of triangles in this small node. Determines length of triangleMask
        unsigned char               *triangleMask   = NULL ;    // bit mask indicating if triangles[x] is in this node
    } data ;
  
    kdTreeNode *next	        = NULL ;    // Next kdTreeNode in the list
    kdTreeNode *leftChild	    = NULL ;    // ptr to left child node
    kdTreeNode *rghtChild	    = NULL ;    // ptr ro rght child node
    kdTreeNode *smallroot       = NULL ;    // This is the index into nodelist of the root of the smalllist tree for this node

    // Define class methods here
    //
    kdTreeNode();                           // Default constructor
    kdTreeNode(const kdTreeNode *node);            // constructor for initialising wth another node
    kdTreeNode(std::vector<int> tris) ;		// constructor using a vector array of ints representing triangle indices
    kdTreeNode(std::string plyFileName) ;	// constructor using a triangle mesh in a .ply file
//    kdTreeNode(const kdTreeNode &node) ;    // copy constructor as we malloc in this class
    ~kdTreeNode() ;							// default destructor
    
    AABB BVforAllTris();
    void medianSplit(kdTreeNode **left, kdTreeNode **rght);
    void split(int dim, float pos, kdTreeNode *left, kdTreeNode *rght) ;
    
};

extern std::vector<kdTreeNode *> nodelist ;


// Linked list to hande the kdTree nodes
// We are using this instead of the std vector<> or list<> classes as we
// need a bit more 'flexibility' (sic) to track nodes in memory
//
class treeList {
private :
    kdTreeNode *head ;
    kdTreeNode *tail ;
    int length = 0 ;
    bool indexCheck(int position) ;
    
public :
    treeList() ;
    int size() ;
    bool insertNode(kdTreeNode *node, int position) ;
    bool removeNode(int position) ;
    kdTreeNode * at(int ind) ;
    kdTreeNode * front() ;
    kdTreeNode * back() ;
    void push_back(kdTreeNode *node);
    void attach(kdTreeNode *node);
    void clear() ;
    void erase() ;
    void printList() ;
    bool empty() ;
    ~treeList() ;
} ;

void swapLists(treeList **a, treeList **b) ;
void printKdTreeNodes(std::vector<kdTreeNode *> nodelist);
void printKdTreeData() ;

#endif /* kdTreeNode_hpp */
