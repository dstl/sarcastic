/***************************************************************************
 *
 *       Module:    buildTree.hpp
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
#ifndef buildTree_hpp
#define buildTree_hpp

#include <stdio.h>
#include <sys/stat.h>
#include "kdTreeNode.hpp"
#include "splitCandidate.hpp"
#include "TriangleMesh.hpp"

namespace kdTree {
    
    
#define TRAVERSALCOST ((float)(15.0))
#define INTERSECTIONCOST ((float)(20.0))
#define SMALLSIZE (64)
#define Ce (0.25)           // Percentage empty space in a node
#define MAXTREELEVEL 10
#define KDT_ISLEAF(n)       ((n->leaf.leafDim & (unsigned char)(0x4))>>2)
#define KDT_DIMENSION(n)    (n->leaf.leafDim & 0x3)
#define KDT_NUMTRIS(n)      (n->leaf.ntriangles)
#define KDT_LEFTCHILD(n)    (n->brch.leftaddress)
#define KDT_RGHTCHILD(n)    (n->brch.rghtaddress)
#define KDT_SPLITPOS(n)     (n->leaf.splitPosition)
#define KDT_TRIINDEX(n)     (n->leaf.triangleIndex)
#define KDT_INDEX(n)        (n->leaf.triangleIndex)
#define NILROPE ((int) -666666 )
    
    enum  TREEOUTPUT {
        OUTPUTNO   = 0,
        OUTPUTDATA = 1,
        OUTPUTNODE = 2,
        OUTPUTAABB = 4,
        OUTPUTSUMM = 8
    } ;
    
    void processLargeNodes(treeList **activelist, treeList **smalllist, treeList **nextlist, TriangleMesh &globalMesh);
    void preProcessSmallNodes(treeList **smalllist) ;
    int reduce(unsigned char *list, long int size) ;
    int reduce(std::vector<int> list);
    void processSmallNodes(treeList **activelist, treeList **nextlist);
    void buildSizes(std::vector<kdTreeNode *> *nodelist, int level) ;
    void buildAddresses(std::vector<kdTreeNode *> *nodelist, int level);
    void scanInclusive(int *in, int *out, int n) ;
    void scanExclusive(int *in, int *out, int n) ;
    void preOrderTraversalNode(std::vector<kdTreeNode *> *nodelist, KdData **tree, int *numNodesInTree) ;
    void dumpsmall(treeList **list) ;
    void writeAABBtoPlyFile(AABB bv, std::string filename);
    void writeAABBtoPlyFile(std::vector<AABB> abs, std::string filename);
    void writeAABBtoPlyFile(std::vector<AABB> abs, std::string filename, int r, int g, int b);
    void swapLists(treeList **a, treeList **b) ;
    void printKdTreeNodes(std::vector<kdTreeNode *> nodelist);
    void printKdTreeData(KdData **kdTree, int *numNodesInTree) ;
    void printSummary(std::vector<kdTreeNode *> nodelist, TriangleMesh &mesh) ;

    void buildTree(TriangleMesh &mesh, KdData **kdTree, int *numNodesInTree, TREEOUTPUT output) ;
    
//    extern TriangleMesh globalMesh;
//    extern std::vector<kdTree::kdTreeNode *> nodelist ;
    
}

#endif /* buildTree_hpp */
