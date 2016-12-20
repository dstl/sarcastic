//
//  main.cpp
//  fastKdTree
//
//  Created by Darren on 23/10/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#include <iostream>
#include <SIlib2/SIlib2.h>
#include "colourCodes.h"
#include "TriangleMesh.hpp"
#include "kdTreeNode.hpp"
#include "splitCandidate.hpp"

#define ROOTPATH "/tmp"
#define TRAVERSALCOST ((float)(15.0))
#define INTERSECTIONCOST ((float)(20.0))
#define SMALLSIZE (64)
#define Ce (0.25)

using namespace std;

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def);
void processLargeNodes(vector<kdTreeNode> *activelist, vector<kdTreeNode> *smalllist, vector<kdTreeNode> *nextlist);
void preProcessSmallNodes(vector<kdTreeNode> &smalllist) ;
int reduce(std::vector<unsigned char> list);
int reduce(std::vector<int> list);
void processSmallNodes(vector<kdTreeNode> &activelist, vector<kdTreeNode> nextlist);
void UpPass(vector<kdTreeNode> &activelist, int level) ;
void scanInclusive(int *in, int *out, int n) ;

int main(int argc, const char * argv[]) {
    
    // Initialise SILib
    //
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, (char *)"buildKDTree", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status);
    
    // Read in triangle data
    //
    char *instr = tryReadFile((char *)"Input triangle filename", (char *)"infilename",
                              (char *)"The name of a triangle .ply file containing triangle data",
                              (char *) ROOTPATH"/delaunay.ply");
    
    char *oustr  = input_string((char *)"KdTree Filename", (char *)"kdTreeFname",
                                (char *)"Full path name where KdTree will be stored",
                                (char *) ROOTPATH"/KdTree.kdt");
    FILE *fpout ;
    fpout = fopen(oustr, "w");
    if (fpout == NULL) {
        printf("Error : failed to open input triangle file %s\n",oustr);
        exit(1);
    }
    
    // Perform the Zhou KdTree build
    //
    vector<kdTreeNode> nodelist ;
    vector<kdTreeNode> activelist ;
    vector<kdTreeNode> smalllist ;
    vector<kdTreeNode> nextlist ;
    kdTreeNode rootnode((string(instr))) ;
    activelist.push_back(rootnode);
    
    while(!activelist.empty()){
        for(auto it=activelist.begin(); it<activelist.end(); ++it){
            nodelist.push_back(*it);
        }
        nextlist.clear() ;
        processLargeNodes(&activelist, &smalllist, &nextlist) ;
        nextlist.swap(activelist);
    }
    
    preProcessSmallNodes(smalllist);
    activelist = smalllist ;
    
    while(!activelist.empty()){
        for(auto it=activelist.begin(); it<activelist.end(); ++it){
            nodelist.push_back(*it);
        }
        nextlist.clear() ;
        processSmallNodes(activelist, nextlist);
        nextlist.swap(activelist);
    }
    
    
    preOrderTraversal(nodelist);
    
    // Write out data to KdTree file
    //
    
    return 0;
}

void processLargeNodes(vector<kdTreeNode> *activelist, vector<kdTreeNode> *smalllist, vector<kdTreeNode> *nextlist)
{
    vector<AABB> tightAABBs;
    AABB aabb ;
    // Compute the AABB for each node in activelist
    //
    for(auto it=activelist->begin(); it!= activelist->end(); ++it){
        aabb = it->boundingVolume() ;
        tightAABBs.push_back(aabb) ;
    }
    
    // split large nodes
    //
    // Start by clipping away empty space in node around triangles that are in the node
    //
    double s; // size of node
    double ldist, hdist;
    kdTreeNode node;
    for (int i=0; i<activelist->size(); ++i) {
        node = activelist->at(i);
        for(int j=0; j<3; ++j){
            s = node.aabb.BB.cell[j] - node.aabb.AA.cell[j] ;
            aabb = tightAABBs[i] ;
            ldist = aabb.AA.cell[j] - node.aabb.AA.cell[j] ;
            if( (ldist/s) > Ce){
                node.aabb.AA.cell[j] = aabb.AA.cell[j] ;
            }
            hdist = node.aabb.BB.cell[j] - aabb.BB.cell[j] ;
            if( (hdist/s) > Ce){
                node.aabb.BB.cell[j] = aabb.BB.cell[j] ;
            }
        }
        
        kdTreeNode leftNode;
        kdTreeNode rghtNode;
        node.medianSplit(leftNode, rghtNode) ;
        
        int numLeftChild  = (int)leftNode.triangles.size();
        int numRightChild = (int)rghtNode.triangles.size();
        
        if (numLeftChild < SMALLSIZE) {
            smalllist->push_back(leftNode) ;
        }else{
            nextlist->push_back(leftNode) ;
        }
        if (numRightChild < SMALLSIZE) {
            smalllist->push_back(rghtNode) ;
        }else{
            nextlist->push_back(rghtNode) ;
        }
    }
    return ;
}

void preProcessSmallNodes(vector<kdTreeNode> &smalllist)
// This function generates a list of split candidates for each node
// in smalllist.
{
    
    int ntris;
    
    for(int i=0; i< smalllist.size(); ++i){
        ntris = (int)smalllist[i].triangles.size() ;
        
        smalllist[i].triangleMask.reserve(ntris) ;
        
        for(int j=0; j<ntris; ++j){
            smalllist[i].triangleMask[j] = 1;
        }
        
        for (int k=0; k<3; ++k){
            
            for(int j=0; j<ntris; ++j){
                splitCandidate eA(smalllist[i].triAABBs[j].AA.cell[k], j, k, ntris) ;
                splitCandidate eB(smalllist[i].triAABBs[j].BB.cell[k], j, k, ntris) ;
                
                smalllist[i].splitList.push_back(eA);
                smalllist[i].splitList.push_back(eB);
            }
        }
        
        for(int j=0; j<smalllist[i].splitList.size(); ++j){
            
            int k = smalllist[i].splitList[j].dim ;
            
            for(int t=0; t<ntris; ++t){
                smalllist[i].splitList[j].leftTris[t] =
                (smalllist[i].triAABBs[t].AA.cell[k] <= smalllist[i].splitList[j].pos) ;
                smalllist[i].splitList[j].rghtTris[t] =
                (smalllist[i].triAABBs[t].BB.cell[k] >= smalllist[i].splitList[j].pos) ;
            }
        }
        
    }
    return ;
}

void processSmallNodes(vector<kdTreeNode> &activelist, vector<kdTreeNode> nextlist)
{
    for(int i=0; i<activelist.size(); ++i){
        float A0 = activelist[i].aabb.surfaceArea() ;
        int minId = -1;
        
        vector<int> s = activelist[i].triangleMask ;
        int ntrisInThisNode = reduce(activelist[i].triangleMask);
        float SAH0 = INTERSECTIONCOST * ntrisInThisNode ;
        float minSAH = SAH0;
        
        for (int j=0; j<activelist[i].splitList.size(); ++j){
            
            // For this split candidate find the number of triangles to the
            // left and right of the split plane that are also in this node
            //
            int tri = activelist[i].splitList[j].owner ;
            if(activelist[i].triangleMask[tri] != 0){
                
                int Cl = reduce(activelist[i].splitList[j].leftTris) ;
                int Cr = reduce(activelist[i].splitList[j].rghtTris) ;
                SPVector BB = activelist[i].aabb.BB;
                BB.cell[activelist[i].splitList[j].dim] = activelist[i].splitList[j].pos ;
                AABB Vleft = AABB(activelist[i].aabb.AA, BB);
                AABB Vrght = AABB(BB, activelist[i].aabb.BB);
                float Al = Vleft.surfaceArea();
                float Ar = Vrght.surfaceArea();
                float ProbL = Al / A0 ;
                float ProbR = Ar / A0 ;
                float SAHp = TRAVERSALCOST + ((Cl * ProbL + Cr * ProbR) * INTERSECTIONCOST);
                if(SAHp < minSAH){
                    minSAH = SAHp ;
                    minId = j;
                }
            }
        }
        
        if(minId == -1){
            activelist[i].isLeaf = true ;
        }else{
            kdTreeNode leftNode;
            kdTreeNode rghtNode;
            activelist[i].split(activelist[i].splitList[minId].dim, activelist[i].splitList[minId].pos, leftNode, rghtNode) ;
            nextlist.push_back(leftNode);
            nextlist.push_back(rghtNode);
        }
    }
    return ;
}

void adoption(vector<kdTreeNode> &nodelist)
{
    int *a = new int [nodelist.size()] ;
    int *b = new int [nodelist.size()] ;
    int *c = new int [nodelist.size()] ;
    
    for (int i=0; i< nodelist.size(); ++i){
        a[i] = ( ! nodelist[i].isLeaf ) ;
    }
    
    scanInclusive(a, b, (int)nodelist.size() );
    
    for (int i=0; i< nodelist.size(); ++i){
        c[i] = a[i] * b[i] ;
        nodelist[i].leftAddress = 2*c[i] - 1;
    }
    
    return ;
}


void preOrderTraversalNode(vector<kdTreeNode> &nodelist)
{
    // Connect nodelist pointers using adoption algorithm
    //
    adoption( nodelist );
    
    kdTreeNode node = *nodelist.end() ;
    int level = node.level ;
    
    for (int i=level-1; i>= 0; --i){
        UpPass(nodelist, i) ;
    }
    // Allocate size of tree based upon size just calculated
    //
    
    
}

void UpPass(vector<kdTreeNode> &nodelist, int level){
    
    for(int i=0; i<nodelist.size(); ++i){
        if (!nodelist[i].isLeaf) {
            nodelist[i].size = nodelist[nodelist[i].leftAddress].size + nodelist[nodelist[i].leftAddress+1].size + 1 ;
        }else{
            nodelist[i].size = (int)nodelist[i].triangles.size() + 1 ;
        }
    }
}

int reduce(std::vector<unsigned char> list){
    int sum =0;
    for(int i=0; i<list.size(); ++i){
        sum += list[i] ;
    }
    return sum;
}

int reduce(std::vector<int> list){
    int sum =0;
    for(int i=0; i<list.size(); ++i){
        sum += list[i] ;
    }
    return sum;
}

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def)
///  prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  help    : the text to display to the user when they just enter '?'
///  def     : the default, ie the value to take if the user just presses return
///
{
    FILE *fp;
    SPStatus fileStat ;
    char *prmpt, *fname;
    
    size_t len = strlen(def);
    prmpt = (char *)sp_malloc(sizeof(char) * len);
    
    do {
        im_init_status(fileStat, 0) ;
        strcpy(prmpt, def);
        fname = input_string(prompt, key, help, def);
        
        if ( (fp = fopen(fname, "r")) == NULL){
            printf(RED "Cannot access file %s\n" RESETCOLOR,fname);
            fileStat.status = BAD_FILE ;
        }
        
    } while(fileStat.status != NO_ERROR);
    
    fclose(fp);
    free(prmpt) ;
    return(fname) ;
}

// Kernels for reductions
//

void scanInclusive(int *in, int *out, int n){
    out[0] = in[0];
    for (int k=1; k<n; ++k){
        out[k] = in[k] + out[k-1] ;
    }
}

float reduceSum(float *arr, int nProcs)
{
    int a1,a2 ;
    for(int p=0; p<nProcs; ++p){
        for(int i=0; i<log2(nProcs)-1; i++){
            a1 = pow(2,i+1)*p;
            a2 = a1 + pow(2,i);
            if(a2<nProcs){
                arr[a1] = arr[a1] + arr[a2] ;
            }
        }
    }
    return arr[a1] ;
}

float reduceMin(float *arr, int nProcs)
{
    int a1,a2 ;
    for(int p=0; p<nProcs; ++p){
        for(int i=0; i<log2(nProcs)-1; i++){
            a1 = pow(2,i+1)*p;
            a2 = a1 + pow(2,i);
            if(a2<nProcs){
                arr[a1] = (arr[a1] < arr[a2]) ? arr[a1] : arr[a2] ;
            }
        }
    }
    return arr[a1] ;
}
float reduceMax(float *arr, int nProcs)
{
    int a1,a2 ;
    for(int p=0; p<nProcs; ++p){
        for(int i=0; i<log2(nProcs)-1; i++){
            a1 = pow(2,i+1)*p;
            a2 = a1 + pow(2,i);
            if(a2<nProcs){
                arr[a1] = (arr[a1] > arr[a2]) ? arr[a1] : arr[a2] ;
            }
        }
    }
    return arr[a1] ;
}
void segReduceMin(float *data, int *owner, int nElements, float **results)
{
    int p0, p1, w0,w1 ;
    for (int d=0; d<log2(nElements)-1; ++d){
        p0 = pow(2,d) ;
        p1 = pow(2,d+1) ;
        for (int i=0; i < ((nElements-1) / p1); ++i) {
            w0 = owner[p1*i] ;
            w1 = owner[p1*i+p0] ;
            if (w0 != w1){
                (*results)[w1] = ((*results)[w1] < data[p1*i + p0]) ? (*results)[w1] : data[p1*i + p0] ;
            }else{
                data[p1*i] = (data[p1*i] < data[p1*i + p0]) ? data[p1*i] : data[p1*i + p0] ;
            }
        }
    }
    return ;
}

void segReduceMax(float *data, int *owner, int nElements, float **results)
{
    int p0, p1, w0,w1 ;
    for (int d=0; d<log2(nElements)-1; ++d){
        p0 = pow(2,d) ;
        p1 = pow(2,d+1) ;
        for (int i=0; i < ((nElements-1) / p1); ++i) {
            w0 = owner[p1*i] ;
            w1 = owner[p1*i+p0] ;
            if (w0 != w1){
                (*results)[w1] = ((*results)[w1] > data[p1*i + p0]) ? (*results)[w1] : data[p1*i + p0] ;
            }else{
                data[p1*i] = (data[p1*i] > data[p1*i + p0]) ? data[p1*i] : data[p1*i + p0] ;
            }
        }
    }
    return ;
}


