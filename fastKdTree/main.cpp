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
#define SMALLSIZE (72)
#define Ce (0.25)

using namespace std;

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def);
void processLargeNodes(treeList **activelist, treeList **smalllist, treeList **nextlist);
void preProcessSmallNodes(treeList **smalllist) ;
int reduce(unsigned char *list, long int size) ;
int reduce(std::vector<int> list);
void processSmallNodes(treeList **activelist, treeList **nextlist);
void buildSizes(std::vector<kdTreeNode *> *nodelist, int level) ;
void buildAddresses(std::vector<kdTreeNode *> *nodelist, int level);
void scanInclusive(int *in, int *out, int n) ;
void scanExclusive(int *in, int *out, int n) ;
void indexTriangles( treeList **nodelist);
void preOrderTraversalNode(std::vector<kdTreeNode *> *nodelist) ;
void writeKdTreeToFile(string filename, std::vector<kdTreeNode *> *nodelist) ;
void writeAABBtoPlyFile(AABB bv, std::string filename);

void dumpsmall(treeList **list){
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

int main(int argc, const char * argv[]) {
    
    kdTreeNode *p, *q;      // used for indexing through the linked list
    
    // Initialise SILib
    //
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, (char *)"fastKDTree", argc, (char **)argv);
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
    
    kdTreeNode rootnode((string(instr))) ;
    
    printf("Input Mesh has %ld triangles\n",globalMesh.triangles.size()) ;
    
    activelist->push_back(&rootnode);
    
    // Initialise and start timer
    //
    Timer runTimer ;
    startTimer(&runTimer, &status) ;
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
    
    printKdTreeNodes(nodelist) ;
    dumpsmall(&smalllist);

    
    printf("PreProcessing Small Nodes ...\n");
    preProcessSmallNodes(&smalllist);
    printf("Smalllist:\n");
    dumpsmall(&smalllist);
    
    swapLists(&smalllist, &activelist);

    
    delete smalllist ;

    printf("Done!\n");
    
    printKdTreeNodes(nodelist) ;

    
    printf("Processing Small Nodes...\n");
    dpth=0;
    while(!activelist->empty()){
        
       

        nextlist->clear() ;
        
        printf("[%02d] activelist is:\n",dpth);
        dumpsmall(&activelist);
        
        processSmallNodes(&activelist, &nextlist);
        
        printf("[%02d] nextlist is:\n",dpth);
        dumpsmall(&nextlist);
        
        // Copy activelist into nodelist.
        //
        q = activelist->front();
        while (q){
            p=q ;
            nodelist.push_back(p) ;
            q=p->next ;
        }
        
        printKdTreeNodes(nodelist);
        
        printf("[%3d] active list size: %4d\tnextlist size: %4d\n",dpth++,activelist->size(),nextlist->size());
        swapLists(&nextlist, &activelist);
    }
    
    printf("Done!\n");
    
    
    /*for (int i=0; i<nodelist.size(); ++i) {
        AABB ab = nodelist[i].aabb ;
        char fn[255] ;
        sprintf(fn, "/tmp/AABBs/aabb_d%02d_%02d.ply",nodelist[i].level,i);
        writeAABBtoPlyFile(ab, string(fn));
    }*/
    
    preOrderTraversalNode(&nodelist) ;
    printKdTreeNodes(nodelist) ;
    printKdTreeData() ;

    // end timer
    //
    endTimer(&runTimer, &status);


    // Write out data to KdTree file
    //
//    writeKdTreeToFile(string(oustr), nodelist) ;
    
#define USECOLOR 0
    if (USECOLOR==1) {
        printf("KdTree constructed in " BOLD BLINK GREEN " %f " RESETCOLOR "seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
    }else{
        printf("KdTree constructed in  %f seconds \n",timeElapsedInSeconds(&runTimer, &status)) ;
    }

    nodelist.clear() ;
    smalllist->clear() ;
    nextlist->clear() ;
    activelist->erase() ;

    
    return 0;
}

//void writeKdTreeToFile(string filename, treeList *nodelist)
//{
//    FILE *fp;
//    double AAx,AAy,AAz,BBx,BBy,BBz;
//    kdTreeNode *p,*q,*n;
//    
//    fp = fopen(filename.c_str(), "wb");
//    if( fp==NULL){
//        printf("ERROR: Failed to open file %s for writing\n",filename.c_str()) ;
//        exit(1);
//    }
//    printf("Writing KdTree to file %s...\n", filename.c_str());
//    
//    n = nodelist->front() ;
//    AABB aabb = n->BVforAllTris() ;
//    AAx = aabb.AA.x ;
//    AAy = aabb.AA.y ;
//    AAz = aabb.AA.z ;
//    BBx = aabb.BB.x ;
//    BBy = aabb.BB.y ;
//    BBz = aabb.BB.z ;
//    
//    fwrite(&AAx,sizeof(double),1,fp);
//    fwrite(&AAy,sizeof(double),1,fp);
//    fwrite(&AAz,sizeof(double),1,fp);
//    fwrite(&BBx,sizeof(double),1,fp);
//    fwrite(&BBy,sizeof(double),1,fp);
//    fwrite(&BBz,sizeof(double),1,fp);
//    
//    int ntri = (int)globalMesh.triangles.size() ;
//    fwrite(&ntri,sizeof(int),1,fp);
//    
//    SPVector aa,bb,cc,nn ;
//    double nx,ny,nz, nd_u, nd_v, d, denom, kbu, kbv, kbd, kcu, kcv, kcd ;
//    int k,u,v,tex;
//    const unsigned int quickmodulo[] = {0,1,2,0,1};
//    
//    printf("Packing %d triangles\n",ntri);
//    
//    for(int i=0; i<ntri; ++i){
//        
//        aa = globalMesh.vertices[globalMesh.triangles[i].a].asSPVector() ;
//        bb = globalMesh.vertices[globalMesh.triangles[i].b].asSPVector() ;
//        cc = globalMesh.vertices[globalMesh.triangles[i].c].asSPVector() ;
//        
//        // Create accelerated triangle structure
//        //
//        SPVector ab; VECT_SUB(bb, aa, ab);
//        SPVector bc; VECT_SUB(cc, bb, bc);
//        SPVector ac; VECT_SUB(cc, aa, ac);
//        SPVector x;  VECT_CROSS(ab, bc, x);
//        double l ; l = VECT_MAG(x);
//        VECT_SCMULT(x, (1/l), nn);
//        
//        nx = fabs(nn.x) ;
//        ny = fabs(nn.y) ;
//        nz = fabs(nn.z) ;
//        
//        if( nx > ny ){
//            if (nx > nz) k = 0; /* X */ else k=2; /* Z */
//        }else{
//            if ( ny > nz) k=1; /* Y */ else k=2; /* Z */
//        }
//        u = quickmodulo[k+1];
//        v = quickmodulo[k+2];
//        nd_u = nn.cell[u] / nn.cell[k] ;
//        nd_v = nn.cell[v] / nn.cell[k] ;
//        d =  (aa.x * ( nn.x / nn.cell[k] )) + (aa.y* ( nn.y / nn.cell[k])) + (aa.z* (nn.z / nn.cell[k])) ;
//        denom = 1./((ac.cell[u]*ab.cell[v]) - (ac.cell[v]*ab.cell[u]));
//        kbu     = -ac.cell[v] * denom;
//        kbv     =  ac.cell[u] * denom;
//        kbd     =  ((ac.cell[v] * aa.cell[u]) - (ac.cell[u] * aa.cell[v])) * denom;
//        kcu     =  ab.cell[v] * denom;
//        kcv     = -ab.cell[u] * denom;
//        kcd     =  ((ab.cell[u] * aa.cell[v]) - (ab.cell[v] * aa.cell[u])) * denom;
//        tex     =  globalMesh.triangles[i].mat ;
//        
//        fwrite(&d,sizeof(double),1,fp);
//        fwrite(&nd_u,sizeof(double),1,fp);
//        fwrite(&nd_v,sizeof(double),1,fp);
//        fwrite(&k,sizeof(int),1,fp);
//        fwrite(&kbu,sizeof(double),1,fp);
//        fwrite(&kbv,sizeof(double),1,fp);
//        fwrite(&kbd,sizeof(double),1,fp);
//        fwrite(&kcu,sizeof(double),1,fp);
//        fwrite(&kcv,sizeof(double),1,fp);
//        fwrite(&kcd,sizeof(double),1,fp);
//        fwrite(&tex, sizeof(int), 1, fp);
//        
//    }
//    
//    // Calculate the number of leaves and then for each leaf write the index of each triangle in the leaf
//    //
//    
//    int nleaves = 0;
//    int indexOfTri ;
//    int nnodetris;
//    
//    p = nodelist->front();
//    q = nodelist->front();
//    while (q){
//        p=q ;
//
//        if (p->data.isLeaf) {
//            nleaves++;
//        }
//        
//        q=p->next ;
//    }
//    
//    fwrite(&nleaves,sizeof(int),1,fp);
//    printf("Packing %d Leaves\n",nleaves);
//    
//    p = nodelist->front();
//    q = nodelist->front();
//    int i = 0;
//    while (q){
//        p=q ;
//        
//        if (p->data.isLeaf) {
//            nnodetris = kdTreeTriangleIndicesOutput[p->data.triangleIndex] ;
//            fwrite(&nnodetris,sizeof(int),1,fp);
//            for(int j=0; j<nnodetris; ++j){
//                indexOfTri = kdTreeTriangleIndicesOutput[p->data.triangleIndex+j+1] ;
//                fwrite(&indexOfTri,sizeof(int),1,fp);
//            }
//        }
//        
//        q=p->next ;
//        i++;
//    }
//    
//    long nnodes = nodelist->size() ;
//    
//    printf("Packing %ld Nodes\n",nnodes);
//    fwrite(&nnodes,sizeof(int),1,fp);
//
//    for(int i=0; i<nodelist->size(); ++i){
//        fwrite(&(kdTreeOutput[i].branch.flagDimAndOffset), sizeof(unsigned int),1,fp);
//        fwrite(&(kdTreeOutput[i].branch.splitPosition),sizeof(float),1,fp);
//    }
//    
//    printf("Writing triangle coords\n");
//    double x,y,z ;
//    double localToGlobal[9] ;
//    double globalToLocal[9] ;
//    for(int i=0; i<ntri; ++i){
//        x = globalMesh.vertices[globalMesh.triangles[i].a].x ;
//        y = globalMesh.vertices[globalMesh.triangles[i].a].y ;
//        z = globalMesh.vertices[globalMesh.triangles[i].a].z ;
//        fwrite(&x,sizeof(double),1,fp);
//        fwrite(&y,sizeof(double),1,fp);
//        fwrite(&z,sizeof(double),1,fp);
//        
//        x = globalMesh.vertices[globalMesh.triangles[i].b].x ;
//        y = globalMesh.vertices[globalMesh.triangles[i].b].y ;
//        z = globalMesh.vertices[globalMesh.triangles[i].b].z ;
//        fwrite(&x,sizeof(double),1,fp);
//        fwrite(&y,sizeof(double),1,fp);
//        fwrite(&z,sizeof(double),1,fp);
//        
//        x = globalMesh.vertices[globalMesh.triangles[i].c].x ;
//        y = globalMesh.vertices[globalMesh.triangles[i].c].y ;
//        z = globalMesh.vertices[globalMesh.triangles[i].c].z ;
//        fwrite(&x,sizeof(double),1,fp);
//        fwrite(&y,sizeof(double),1,fp);
//        fwrite(&z,sizeof(double),1,fp);
//        
//       	// get normal
//       	//
//       	x = globalMesh.triangles[i].N.x ;
//       	y = globalMesh.triangles[i].N.y ;
//       	z = globalMesh.triangles[i].N.z ;
//        fwrite(&x,sizeof(double),1,fp);
//        fwrite(&y,sizeof(double),1,fp);
//        fwrite(&z,sizeof(double),1,fp);
//        
//        // get area
//        //
//        x = globalMesh.triangles[i].Area ;
//        fwrite(&x,sizeof(double),1,fp);        
//        
//        // get global to local matrix and local to global matrix
//        //
//        globalMesh.triangles[i].matrices(localToGlobal, globalToLocal) ;
//        for (int j=0; j<9; j++){
//            fwrite(&(globalToLocal[j]),sizeof(double),1,fp);
//        }
//        for (int j=0; j<9; j++){
//            fwrite(&(localToGlobal[j]),sizeof(double),1,fp);
//        }
//        
//        // get material
//        //
//        fwrite(&(globalMesh.triangles[i].mat), sizeof(int), 1, fp);
//    }
//    
//    fclose(fp);
//    
//    return ;
//}

void processLargeNodes(treeList **activelist, treeList **smalllist, treeList **nextlist)
{
    
    // BV4TrisInNode is the AABB that tightly bounds all the triangles in the node
    //
    vector<AABB> BV4TrisInNode;
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

void preProcessSmallNodes(treeList **smalllist)
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

/*
 void preProcessSmallNodes(vector<kdTreeNode> &smalllist, long int smalllistOffset)
// This function generates a list of split candidates for each node
// in smalllist.
// smalllistOffset is the index of the first small node in the output nodelist
//
{
    for(int i=0; i< smalllist.size(); ++i){
        
        // set the index for the smallroot in the final nodelist so that we can refer back to its
        // triangle contents later. This allows us to just use a triangle mask for all the children
        // of the smallroot
        //
        smalllist[i].smallroot  = (int)smalllistOffset + i ;
        smalllist[i].smallntris = (int)smalllist[i].triangles.size() ;
        
        // Set the triangle mask to be all '1's to start with showing that each triangle in
        // triangles is actually in this smallnode
        //
        if(smalllist[i].smallntris <=0){
            printf("error : smallroot with zero triangles in preProcessSmallNodes\n");
            exit(1);
        }
        if(smalllist[i].smallroot == -1){
            printf("Error: node %d in activelist is not a small node in function preProcessSmallNodes()\n",i) ;
            exit(1);
        }
        smalllist[i].triangleMask = new unsigned char [smalllist[i].smallntris] ;
        for (int j=0; j<smalllist[i].triangles.size(); ++j) {
            smalllist[i].triangleMask[j] = 1 ;
        }
        
        // Create split candidates based upon the AABBs of the triangles in this node
        // two for each triangle (low and high), and for each dimension.
        // Store them all in splitList for this node
        //
        for (int k=0; k<3; ++k){
            
            for(int j=0; j<smalllist[i].smallntris; ++j){
                splitCandidate eA(smalllist[i].triAABBs[j].AA.cell[k], j, k, smalllist[i].smallntris) ;
                splitCandidate eB(smalllist[i].triAABBs[j].BB.cell[k], j, k, smalllist[i].smallntris) ;
                
                smalllist[i].splitList.push_back(eA);
                smalllist[i].splitList.push_back(eB);
            }
        }
        
        for(int j=0; j<smalllist[i].splitList.size(); ++j){
            
            int k = smalllist[i].splitList[j].dim ;
            
            for(int t=0; t<smalllist[i].smallntris; ++t){
                smalllist[i].splitList[j].leftTris[t] = (smalllist[i].triAABBs[t].AA.cell[k] <= smalllist[i].splitList[j].pos) ;
                smalllist[i].splitList[j].rghtTris[t] = (smalllist[i].triAABBs[t].BB.cell[k] >= smalllist[i].splitList[j].pos) ;
            }
        }
    }
    return ;
}
 */

void processSmallNodes(treeList **activelist, treeList **nextlist)
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
            
//            printf("leftnode:");
//            for(int i=0; i<leftNode.data.smallntris; ++i){
//                printf("%d", leftNode.data.triangleMask[i]);
//            }
//            printf("\n");
//            printf("rghtNode:");
//            for(int i=0; i<rghtNode.data.smallntris; ++i){
//                printf("%d", rghtNode.data.triangleMask[i]);
//            }
//            printf("\n");
            
//            kdTreeNode left = leftNode;
//            kdTreeNode rght = rghtNode;
//            
//            printf("left:");
//            for(int j=0; j<left.data.smallntris; ++j){
//                printf("%d",left.data.triangleMask[j]) ;
//            }
//            printf("\n");
//            printf("rght:");
//            for(int j=0; j<rght.data.smallntris; ++j){
//                printf("%d",rght.data.triangleMask[j]) ;
//            }
//            printf("\n");
            
            p->leftChild = leftNode ;
            p->rghtChild = rghtNode ;
            
            (*nextlist)->push_back(leftNode);
            (*nextlist)->push_back(rghtNode);
            
//            dumpsmall(nextlist);
        }

        
        q=p->next ;
        ++i;
    }
    
    return ;
}

//void adoption(vector<kdTreeNode> &nodelist)
//{
//    int *a = new int [nodelist.size()] ;
//    int *b = new int [nodelist.size()] ;
//    int *c = new int [nodelist.size()] ;
//    
//    for (int i=0; i< nodelist.size(); ++i){
//        a[i] = ( ! nodelist[i].isLeaf ) ;
//    }
//    
//    scanInclusive(a, b, (int)nodelist.size() );
//    
//    for (int i=0; i< nodelist.size(); ++i){
//        c[i] = a[i] * b[i] ;
//        nodelist[i].leftAddress = 2*c[i] - 1;
//    }
//    
//    // Useful DEBUG to print out the Adoption algorithm output
//    //
//    /*printf("ind:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",i);
//    }
//    printf("\n");
//    
//    printf("lef:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",nodelist[i].isLeaf);
//    }
//    printf("\n");
//    printf("a[]:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",a[i]);
//    }
//    printf("\n");
//    
//    printf("b[]:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",b[i]);
//    }
//    printf("\n");
//    printf("c[]:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",c[i]);
//    }
//    printf("\n");
//    printf("add:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",nodelist[i].leftAddress);
//    }
//    printf("\n");
//    */
//    
//    return ;
//}


void preOrderTraversalNode(std::vector<kdTreeNode *> *nodelist)
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
    
    printf("\n\n");
    for(int i=0; i<nodelist->size(); ++i){
        printf("[%02d]<%p> leaf:%d levl:%d, address:%d, size:%d, left:<%p> rght:<%p>\n",
               i,(*nodelist)[i],(*nodelist)[i]->data.isLeaf,(*nodelist)[i]->data.level,(*nodelist)[i]->data.address,(*nodelist)[i]->data.size,(*nodelist)[i]->leftChild,(*nodelist)[i]->rghtChild);
    }
    
    numOfKdTreeNodes = nodelist->front()->data.size ;
    kdTreeOutput = new KdData [numOfKdTreeNodes] ;
    
    int address ;
    for(int i=0; i<nodelist->size(); ++i){
        p = (*nodelist)[i] ;
        address = p->data.address ;
        // check
        if( address < 0 || address >= numOfKdTreeNodes){
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
            kdTreeOutput[address].leaf.splitPosition = p->data.splitPos ;
            kdTreeOutput[address].leaf.leafDim = (unsigned char)0x4 ;
            kdTreeOutput[address].leaf.leafDim = (kdTreeOutput[address].leaf.leafDim) | p->data.dim ;
            kdTreeOutput[address].leaf.triangleIndex = 999999999 ;
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
                    kdTreeOutput[address+tcount+1].leaf.leafDim = (unsigned char)0x4 ;
                    kdTreeOutput[address+tcount+1].leaf.leafDim = (kdTreeOutput[address].leaf.leafDim) | p->data.dim ;
                    kdTreeOutput[address+tcount+1].leaf.triangleIndex = triangles[t] ;
                    kdTreeOutput[address+tcount+1].leaf.ntriangles = 0 ;
                    tcount++;
                }
            }
            kdTreeOutput[address].leaf.ntriangles =  tcount;

        }else{
            kdTreeOutput[address].brch.splitPosition = p->data.splitPos ;
            kdTreeOutput[address].brch.leafDim = (unsigned char)0x0 ;
            kdTreeOutput[address].brch.leafDim = (kdTreeOutput[address].leaf.leafDim)| p->data.dim ;
            kdTreeOutput[address].brch.leftaddress = p->leftChild->data.address ;
            kdTreeOutput[address].brch.rghtaddress = p->rghtChild->data.address ;
        }
        
    }
    
    printKdTreeNodes(*nodelist);
    return;
}


//    // Calculate the index to each triangle using the triangleIndex
//    // algorithm and put it in the nodelist
//    //
//    indexTriangles( nodelist ) ;
//    
//    // Allocate size of tree
//    //
//    p = (*nodelist)->front() ;
//    if(p->data.size == 0){
//        printf("Error: tree size is zero\n");
//        exit(1);
//    }
//    
//    numOfTriangleIndices = p->data.size ;
//    kdTreeTriangleIndicesOutput = new int [numOfTriangleIndices] ;
//    
//    // Load triangles into output array
//    //
//    p = (*nodelist)->front() ;
//    q = (*nodelist)->front() ;
//    int i=0;
//    while (q) {
//        p = q;
//    
//        if (p->data.isLeaf) {
//            s = p->smallroot ;
//            // convert triangleMask into indices for each triangle
//            //
//            vector<int> triangleIndices;
//            for(int n=0; n<p->data.smallntris; ++n){
//                if(p->data.triangleMask[n] == 1){
//                    int tri = s->data.triangles[n] ;
//                    triangleIndices.push_back(tri);
//                }
//            }
//            // Now install triangle indices into global output array
//            //
//            kdTreeTriangleIndicesOutput[p->data.triangleIndex] = p->data.size-1 ;
//            if(p->data.triangleIndex > numOfTriangleIndices){
//                printf("ERROR : out of index range exception in kdTreeTriangleIndicesOutput\n");
//                printf("Tried to assign to position %d in %d length array\n",p->data.triangleIndex, numOfTriangleIndices);
//            }
//            printf("node %d has %d triangles\n",i,p->data.size-1) ;
//            // check
//            //
//            if(triangleIndices.size() != p->data.size-1){
//                printf("Something went wrong with triangle counting\n");
//                exit(1);
//            }
//            for(int n=0; n<triangleIndices.size(); ++n){
//                if( p->data.triangleIndex + n + 1 > numOfTriangleIndices){
//                    printf("ERROR Exception : attempting to write beyond end of array\n");
//                    printf("Tried to write to position %d (node[%d].triIndex is %d, n is %d) in array of length %d\n",
//                           p->data.triangleIndex + n + 1, i,p->data.triangleIndex,n,   numOfTriangleIndices ) ;
//                }
//                kdTreeTriangleIndicesOutput[p->data.triangleIndex + n + 1] = triangleIndices[n] ;
//            }
//        }
//        
//        q = p -> next;
//        ++i;
//    }
//
    
//    printKdTreeNodes(*nodelist) ;

//    printf("Triangle Indexing Information\n");
//    printf("-----------------------------\n");
//    for(int i=0; i<nodelist[0].size; ++i){
//        printf("[%03d]  %02d\n",i,kdTreeTriangleIndicesOutput[i]);
//    }
//    printf("KdTree\n");
//    printf("------\n");
//    for(int i=0; i<nodelist.size(); ++i){
//        printf("[%02d]",i);
//        for(int j=0;j<nodelist[i].level;++j){
//            printf("-");
//        }
//        if(nodelist[i].isLeaf){
//            printf("X");
//            printf(" <%02d> #%02d = [",nodelist[i].triangleIndex,kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex]);
//            for(int k=0; k<kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex]; ++k){
//                printf(" %02d",kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex+1+k]);
//            }
//            printf(" ]");
//        }else{
//            printf("|");
//            printf("  [%02d][%02d]",nodelist[i].leftAddress,nodelist[i].leftAddress+1);
//        }
//        printf("\n");
//    }
//    
//    // build KdTree
//    //
//    numOfKdTreeNodes = (int)(*nodelist)->size() ;
//    kdTreeOutput = (KdData *)sp_malloc(sizeof(KdData) * numOfKdTreeNodes) ;
//    unsigned int flagDimAndOffset ;
//    float splitPosition;
//    unsigned int splitDim ;
//    unsigned int offset = 0 ;
//    unsigned int leafoffset = 0;
//    p = (*nodelist)->front() ;
//    q = (*nodelist)->front() ;
//    i=0;
//    while (q) {
//        p = q;
//        
//        splitPosition = p->data.splitPos ;
//        splitDim = p->data.dim ;
//        if (p->data.isLeaf) {
//            flagDimAndOffset = (unsigned int)0x80000000 ;
//            // For compatability with the old buildKdTree Code we just sepcify the
//            // offset as the index of the leaf. We can improve of this by directly writing the
//            // offset as the offset into the triangle array
//            // offset = nodelist[i].triangleIndex ;
//            //
//            flagDimAndOffset = flagDimAndOffset | (leafoffset << 2) ;
//            flagDimAndOffset = flagDimAndOffset | splitDim ;
//            kdTreeOutput[i].leaf.flagDimAndOffset = flagDimAndOffset ;
//            kdTreeOutput[i].leaf.splitPosition    = splitPosition ;
//            leafoffset++;
//        }else{
//            flagDimAndOffset = (unsigned int)0x0 ;
//            offset = p->data.leftAddress ;
//            flagDimAndOffset = flagDimAndOffset | (offset << 2) ;
//            flagDimAndOffset = flagDimAndOffset | splitDim ;
//            kdTreeOutput[i].branch.flagDimAndOffset = flagDimAndOffset ;
//            kdTreeOutput[i].branch.splitPosition    = splitPosition ;
//        }
//        
//        q = p -> next;
//        ++i;
//    }
//    
//    numOfKdTreeNodes = (int)nodelist.size() ;
//    kdTreeOutput = (KdData *)sp_malloc(sizeof(KdData) * numOfKdTreeNodes) ;
//    unsigned int flagDimAndOffset ;
//    float splitPosition;
//    unsigned int splitDim ;
//    unsigned int offset=0 ;
//    unsigned int leafoffset = 0;
//    for(int i=0; i< nodelist.size(); ++i){
//        splitPosition = nodelist[i].splitPos ;
//        splitDim = nodelist[i].dim ;
//        
//        if (nodelist[i].isLeaf) {
//            flagDimAndOffset = (unsigned int)0x80000000 ;
//            // For compatability with the old buildKdTree Code we just sepcify the
//            // offset as the index of the leaf. We can improve of this by directly writing the
//            // offset as the offset into the triangle array
//            // offset = nodelist[i].triangleIndex ;
//            //
//            flagDimAndOffset = flagDimAndOffset | (leafoffset << 2) ;
//            flagDimAndOffset = flagDimAndOffset | splitDim ;
//            kdTreeOutput[i].leaf.flagDimAndOffset = flagDimAndOffset ;
//            kdTreeOutput[i].leaf.splitPosition    = splitPosition ;
//            leafoffset++;
//        }else{
//            flagDimAndOffset = (unsigned int)0x0 ;
//            offset = nodelist[i].leftAddress ;
//            flagDimAndOffset = flagDimAndOffset | (offset << 2) ;
//            flagDimAndOffset = flagDimAndOffset | splitDim ;
//            kdTreeOutput[i].branch.flagDimAndOffset = flagDimAndOffset ;
//            kdTreeOutput[i].branch.splitPosition    = splitPosition ;
//        }
//    }
//    
//    return ;
//}

//void indexTriangles( treeList **nodelist)
//{
//    int *a = new int [(*nodelist).size()] ;
//    int *b = new int [(*nodelist).size()] ;
//    int *c = new int [(*nodelist).size()] ;
//    
//    kdTreeNode *p, *q;
//    p = (*nodelist)->front();
//    q = (*nodelist)->front();
//    int i=0;
//    while(q){
//        p = q ;
//        
//        a[i] = p->data.isLeaf ;
//        b[i] = a[i] * p->data.size ;
//        
//        q = p->next ;
//        ++i;
//    }
//    
//    scanExclusive(b, c, (int)(*nodelist).size() );
//    
//    p = (*nodelist)->front();
//    q = (*nodelist)->front();
//    i=0;
//    while(q){
//        p = q ;
//        p->data.triangleIndex = a[i] * c[i]
//    
//        q=p->next;
//        ++i;
//    }
//    
//    /*
//    printf("siz:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",nodelist[i].size);
//    }
//    printf("\n");
//    printf("b[]:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",b[i]);
//    }
//    printf("\n");
//    printf("c[]:");
//    for (int i=0; i< nodelist.size(); ++i){
//        printf("  %02d",c[i]);
//    }
//    printf("\n");
//    printf("tri:");
//    for (int i=0; i< nodelist.size(); ++i){
//        if(nodelist[i].isLeaf){
//            printf("  %02d",nodelist[i].triangleIndex);
//        }else{
//            printf("  --");
//        }
//    }
//    printf("\n");
//     */
//    
//    return ;
//}


void buildSizes(std::vector<kdTreeNode *> *nodelist, int level){
    
    kdTreeNode *p,*l,*r;
    
    for(int i=0; i<nodelist->size(); ++i){
        p = (*nodelist)[i] ;
        if(p->data.level == level){
            if (!p->data.isLeaf) {
                l = p->leftChild;
                r = p->rghtChild;
                p->data.size = l->data.size + r->data.size + 1 ;
                printf("<%p> node %d size: %d left=%p, right=%p\n",p,i,p->data.size,l,r);
            }else{
                p->data.size = reduce(p->data.triangleMask,p->data.smallntris) + 1 ;
                printf("<%p> node %d size is %d\n",p,i,p->data.size);
            }
        }
    }
    
    return ;
}

void buildAddresses(std::vector<kdTreeNode *> *nodelist, int level){
    
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


int reduce(unsigned char *list, long int size){
    int sum =0;
    for(long int i=0; i<size; ++i){
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

void scanExclusive(int *in, int *out, int n){
    out[0] = 0;
    for (int k=1; k<n; ++k){
        out[k] = in[k-1] + out[k-1] ;
    }
}

float reduceSum(float *arr, int nProcs)
{
    int a1=0,a2 ;
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

int reduceSum(int *arr, int nProcs)
{
    int a1 = 0,a2 ;
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
    int a1 = 0,a2 ;
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
    int a1= 0,a2 ;
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

void writeAABBtoPlyFile(AABB bv, std::string filename){
    
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

