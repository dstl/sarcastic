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
#define SMALLSIZE (16)
#define Ce (0.25)

using namespace std;

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def);
void processLargeNodes(vector<kdTreeNode> *activelist, vector<kdTreeNode> *smalllist, vector<kdTreeNode> *nextlist);
void preProcessSmallNodes(vector<kdTreeNode> &smalllist, long int smalllistOffset) ;
int reduce(unsigned char *list, long int size) ;
int  reduce(std::vector<int> list);
void processSmallNodes(vector<kdTreeNode> &activelist, vector<kdTreeNode> &nextlist);
void buildSizes(vector<kdTreeNode> &nodelist, int level) ;
void scanInclusive(int *in, int *out, int n) ;
void scanExclusive(int *in, int *out, int n) ;
void indexTriangles( vector<kdTreeNode> &nodelist);
void preOrderTraversalNode(vector<kdTreeNode> &nodelist) ;
void writeKdTreeToFile(string filename, vector<kdTreeNode> &nodelist) ;
void writeAABBtoPlyFile(AABB bv, std::string filename);

int main(int argc, const char * argv[]) {
    
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
    vector<kdTreeNode> nodelist ;
    vector<kdTreeNode> activelist ;
    vector<kdTreeNode> smalllist ;
    vector<kdTreeNode> nextlist ;
    kdTreeNode rootnode((string(instr))) ;
//    writeAABBtoPlyFile(rootnode.aabb, string("/tmp/AABBs/root.ply"));

    activelist.push_back(rootnode);
    
    int dpth=0;
    printf("Processing Large Nodes...\n");
    while(!activelist.empty()){
        for(auto it=activelist.begin(); it<activelist.end(); ++it){
            nodelist.push_back(*it);
        }
        nextlist.clear() ;
        processLargeNodes(&activelist, &smalllist, &nextlist) ;
        
//        for (int i=0; i<nextlist.size(); ++i) {
//            AABB ab = nextlist[i].aabb ;
//            char fn[255] ;
//            sprintf(fn, "/tmp/AABBs/aabb_d%02d_%02d.ply",dpth+1,i);
//            writeAABBtoPlyFile(ab, string(fn));
//        }
        printf("[%d] active list size: %ld. Smalllist size %ld, nextlist size: %ld\n",dpth++,activelist.size(),smalllist.size(),nextlist.size());

        nextlist.swap(activelist);
    }
    printf("Done!\n");
    
    printf("PreProcessing Small Nodes ...\n");
    preProcessSmallNodes(smalllist, nodelist.size());
    activelist = smalllist ;
    printf("Done!\n");
    
    for(int i=0; i<activelist.size(); ++i){
        if(activelist[i].smallroot == -1){
            printf("Error: node %d in activelist is not a small node in function main()\n",i) ;
            exit(1);
        }
    }
    
    printf("Processing Small Nodes...\n");
    while(!activelist.empty()){
        
        nextlist.clear() ;
        
        processSmallNodes(activelist, nextlist);
        
        for(auto it=activelist.begin(); it<activelist.end(); ++it){
            nodelist.push_back(*it);
        }
        
        nextlist.swap(activelist);
    }
    
    printf("Done!\n");
    
    
    for (int i=0; i<nodelist.size(); ++i) {
        AABB ab = nodelist[i].aabb ;
        char fn[255] ;
        sprintf(fn, "/tmp/AABBs/aabb_d%02d_%02d.ply",nodelist[i].level,i);
        writeAABBtoPlyFile(ab, string(fn));
    }
    
    preOrderTraversalNode(nodelist) ;

    // Write out data to KdTree file
    //
    writeKdTreeToFile(string(oustr), nodelist) ;
    
    return 0;
}

void writeKdTreeToFile(string filename, vector<kdTreeNode> &nodelist)
{
    FILE *fp;
    double AAx,AAy,AAz,BBx,BBy,BBz;
    
    fp = fopen(filename.c_str(), "wb");
    if( fp==NULL){
        printf("ERROR: Failed to open file %s for writing\n",filename.c_str()) ;
        exit(1);
    }
    printf("Writing KdTree to file %s...\n", filename.c_str());
    
    AABB aabb = nodelist[0].BVforAllTris() ;
    AAx = aabb.AA.x ;
    AAy = aabb.AA.y ;
    AAz = aabb.AA.z ;
    BBx = aabb.BB.x ;
    BBy = aabb.BB.y ;
    BBz = aabb.BB.z ;
    
    fwrite(&AAx,sizeof(double),1,fp);
    fwrite(&AAy,sizeof(double),1,fp);
    fwrite(&AAz,sizeof(double),1,fp);
    fwrite(&BBx,sizeof(double),1,fp);
    fwrite(&BBy,sizeof(double),1,fp);
    fwrite(&BBz,sizeof(double),1,fp);
    
    int ntri = (int)globalMesh.triangles.size() ;
    fwrite(&ntri,sizeof(int),1,fp);
    
    SPVector aa,bb,cc,nn ;
    double nx,ny,nz, nd_u, nd_v, d, denom, kbu, kbv, kbd, kcu, kcv, kcd ;
    int k,u,v,tex;
    const unsigned int quickmodulo[] = {0,1,2,0,1};
    
    printf("Packing %d triangles\n",ntri);
    
    for(int i=0; i<ntri; ++i){
        
        aa = globalMesh.vertices[globalMesh.triangles[i].a].asSPVector() ;
        bb = globalMesh.vertices[globalMesh.triangles[i].a].asSPVector() ;
        cc = globalMesh.vertices[globalMesh.triangles[i].a].asSPVector() ;
        
        // Create accelerated triangle structure
        //
        SPVector ab; VECT_SUB(bb, aa, ab);
        SPVector bc; VECT_SUB(cc, bb, bc);
        SPVector ac; VECT_SUB(cc, aa, ac);
        SPVector x;  VECT_CROSS(ab, bc, x);
        double l ; l = VECT_MAG(x);
        VECT_SCMULT(x, (1/l), nn);
        
        nx = fabs(nn.x) ;
        ny = fabs(nn.y) ;
        nz = fabs(nn.z) ;
        
        if( nx > ny ){
            if (nx > nz) k = 0; /* X */ else k=2; /* Z */
        }else{
            if ( ny > nz) k=1; /* Y */ else k=2; /* Z */
        }
        u = quickmodulo[k+1];
        v = quickmodulo[k+2];
        nd_u = nn.cell[u] / nn.cell[k] ;
        nd_v = nn.cell[v] / nn.cell[k] ;
        d =  (aa.x * ( nn.x / nn.cell[k] )) + (aa.y* ( nn.y / nn.cell[k])) + (aa.z* (nn.z / nn.cell[k])) ;
        denom = 1./((ac.cell[u]*ab.cell[v]) - (ac.cell[v]*ab.cell[u]));
        kbu     = -ac.cell[v] * denom;
        kbv     =  ac.cell[u] * denom;
        kbd     =  ((ac.cell[v] * aa.cell[u]) - (ac.cell[u] * aa.cell[v])) * denom;
        kcu     =  ab.cell[v] * denom;
        kcv     = -ab.cell[u] * denom;
        kcd     =  ((ab.cell[u] * aa.cell[v]) - (ab.cell[v] * aa.cell[u])) * denom;
        tex     =  globalMesh.triangles[i].mat ;
        
        fwrite(&d,sizeof(double),1,fp);
        fwrite(&nd_u,sizeof(double),1,fp);
        fwrite(&nd_v,sizeof(double),1,fp);
        fwrite(&k,sizeof(double),1,fp);
        fwrite(&kbu,sizeof(double),1,fp);
        fwrite(&kbv,sizeof(double),1,fp);
        fwrite(&kbd,sizeof(double),1,fp);
        fwrite(&kcu,sizeof(double),1,fp);
        fwrite(&kcv,sizeof(double),1,fp);
        fwrite(&kcd,sizeof(double),1,fp);
        fwrite(&tex, sizeof(int), 1, fp);
        
    }
    
    // Calculate the number of leaves and then for each leaf write the index of each triangle in the leaf
    //
    
    int nleaves = 0;
    int indexOfTri ;
    for(int i=0; i<nodelist.size(); ++i){
        if (nodelist[i].isLeaf) {
            nleaves++;
        }
    }
    fwrite(&nleaves,sizeof(int),1,fp);
    printf("Packing %d Leaves\n",nleaves);
    
    for(int i=0; i<nodelist.size(); ++i){
        if (nodelist[i].isLeaf) {
            ntri = kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex] ;
            fwrite(&ntri,sizeof(int),1,fp);
            for(int j=0; j<ntri; ++j){
                indexOfTri = kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex+j+1] ;
                fwrite(&indexOfTri,sizeof(int),1,fp);
            }
        }
    }
    
    long nnodes = nodelist.size() ;
    printf("Packing %ld Nodes\n",nnodes);
    
    
    
}

void processLargeNodes(vector<kdTreeNode> *activelist, vector<kdTreeNode> *smalllist, vector<kdTreeNode> *nextlist)
{
    // BV4TrisInNode is the AABB that tightly bounds all the triangles in the node
    //
    vector<AABB> BV4TrisInNode;
    AABB aabb ;
    // Compute the AABB for each node in activelist
    //
    for(auto it=activelist->begin(); it!= activelist->end(); ++it){
        aabb = it->BVforAllTris() ;
        BV4TrisInNode.push_back(aabb) ;
    }
    
    // split large nodes
    //
    // Start by clipping away empty space in node around triangles that are in the node
    //
    double s; // size of node in a given dimansion
    double ldist, hdist; // low and high distances from the BV containing the triangles and the node BV
    kdTreeNode node;
    for (int i=0; i<activelist->size(); ++i) {
        node = activelist->at(i);
        for(int j=0; j<3; ++j){
            s = node.aabb.BB.cell[j] - node.aabb.AA.cell[j] ;
            aabb = BV4TrisInNode[i] ;
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
        printf("tri count- Parent: %ld, leftChild: %ld, rightChild %ld\n",node.triangles.size(),leftNode.triangles.size(),rghtNode.triangles.size());
        
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

void processSmallNodes(vector<kdTreeNode> &activelist, vector<kdTreeNode> &nextlist)
{
    int Cl, Cr ;
    SPVector AA, BB ;
    AABB Vleft,Vrght ;
    float Al,Ar,ProbL,ProbR,SAHp ;
    
    for(int i=0; i<activelist.size(); ++i){
        
        if(activelist[i].smallroot == -1){
            printf("Error: node %d in activelist is not a small node in function processSmallNodes()\n",i) ;
            exit(1);
        }
        
        float A0 = activelist[i].aabb.surfaceArea() ;
        int minId = -1;
        
        int ntrisInThisNode = reduce(activelist[i].triangleMask, activelist[i].smallntris );
        float SAH0 = INTERSECTIONCOST * ntrisInThisNode ;
        float minSAH = SAH0;
        
        for (int j=0; j<activelist[i].splitList.size(); ++j){
            
            // For this split candidate find the number of triangles to the
            // left and right of the split plane that are also in this node
            //
            int tri = activelist[i].splitList[j].owner ;
            if(activelist[i].triangleMask[tri] != 0){
                
                Cl = reduce(activelist[i].splitList[j].leftTris, activelist[i].smallntris) ;
                Cr = reduce(activelist[i].splitList[j].rghtTris, activelist[i].smallntris) ;

                AA = activelist[i].aabb.AA;
                BB = activelist[i].aabb.BB;
                AA.cell[activelist[i].splitList[j].dim] = activelist[i].splitList[j].pos ;
                BB.cell[activelist[i].splitList[j].dim] = activelist[i].splitList[j].pos ;
                Vleft = AABB(activelist[i].aabb.AA, BB);
                Vrght = AABB(AA, activelist[i].aabb.BB);
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
    
    // Useful DEBUG to print out the Adoption algorithm output
    //
    printf("ind:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",i);
    }
    printf("\n");
    
    printf("lef:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",nodelist[i].isLeaf);
    }
    printf("\n");
    printf("a[]:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",a[i]);
    }
    printf("\n");
    
    printf("b[]:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",b[i]);
    }
    printf("\n");
    printf("c[]:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",c[i]);
    }
    printf("\n");
    printf("add:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",nodelist[i].leftAddress);
    }
    printf("\n");
    
    return ;
}


void preOrderTraversalNode(vector<kdTreeNode> &nodelist)
{
    int level ;
    
    // Connect nodelist pointers using adoption algorithm
    //
    adoption( nodelist );
    
    printf("lev:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",nodelist[i].level);
    }
    printf("\n");
    
    // Starting from bottom of tree and working up calculate the size of each node
    //
    if (!nodelist.empty()) {
        level = nodelist.back().level ;
    }else{
        printf("Error : Empty nodelist in preOrderTraversalNode. Exiting...\n");
        exit(1);
    }
    for (int i=level-1; i>= 0; --i){
        buildSizes(nodelist, i) ;
    }
    
    // Calculate the index to each triangle using the triangleIndex
    // algorithm and put it in the nodelist
    //
    indexTriangles( nodelist ) ;
    
    // Allocate size of tree
    //
    if(nodelist[0].size == 0){
        printf("Error: tree size is zero\n");
        exit(1);
    }
    kdTreeTriangleIndicesOutput = new int [nodelist[0].size] ;
    
    // Load triangles into output array
    //
    for(int i=0; i<nodelist.size(); ++i){
        if (nodelist[i].isLeaf) {
            
            kdTreeNode smallrootnode(nodelist[nodelist[i].smallroot]) ;

            // convert triangleMask into indices for each triangle
            //
            vector<int> triangleIndices;
            for(int n=0; n<nodelist[i].smallntris; ++n){
                if(nodelist[i].triangleMask[n] == 1){
                    int tri = smallrootnode.triangles[n] ;
                    triangleIndices.push_back(tri);
                }
            }
            
            // Now install triangle indices into global output array
            //
            kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex] = nodelist[i].size-1 ;
            printf("node %d has %d triangles\n",i,nodelist[i].size-1) ;
            for(int n=0; n<triangleIndices.size(); ++n){
                printf("writing triInd %d to location %d\n",triangleIndices[n],nodelist[i].triangleIndex + n + 1);
                kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex + n + 1] = triangleIndices[n] ;
            }
        }
    }
    
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
//            printf(" <%02d> #%02d = [",nodelist[i].triangleIndex,kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex]);
//            for(int k=0; k<kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex]; ++k){
//                printf(" %02d",kdTreeTriangleIndicesOutput[nodelist[i].triangleIndex+k]);
//            }
//            printf("]");
        }else{
            printf("|");
            printf("  [%02d][%02d]",nodelist[i].leftAddress,nodelist[i].leftAddress+1);
        }
        printf("\n");
    }
    
    // build KdTree
    //
    kdTreeOutput = (KdData *)sp_malloc(sizeof(KdData) * nodelist.size()) ;
    unsigned int flagDimAndOffset ;
    float splitPosition;
    unsigned int splitDim ;
    unsigned int offset ;
    for(int i=0; i< nodelist.size(); ++i){
        splitPosition = nodelist[i].splitPos ;
        splitDim = nodelist[i].dim ;
        
        if (nodelist[i].isLeaf) {
            flagDimAndOffset = (unsigned int)0x80000000 ;
            offset = nodelist[i].triangleIndex ;
            flagDimAndOffset = flagDimAndOffset | (offset << 2) ;
            flagDimAndOffset = flagDimAndOffset | splitDim ;
            kdTreeOutput[i].leaf.flagDimAndOffset = flagDimAndOffset ;
            kdTreeOutput[i].leaf.splitPosition    = splitPosition ;
        }else{
            flagDimAndOffset = (unsigned int)0x0 ;
            offset = nodelist[i].leftAddress ;
            flagDimAndOffset = flagDimAndOffset | (offset << 2) ;
            flagDimAndOffset = flagDimAndOffset | splitDim ;
            kdTreeOutput[i].branch.flagDimAndOffset = flagDimAndOffset ;
            kdTreeOutput[i].branch.splitPosition    = splitPosition ;
        }
    }
    
    return ;
}

void indexTriangles( vector<kdTreeNode> &nodelist)
{
    int *a = new int [nodelist.size()] ;
    int *b = new int [nodelist.size()] ;
    int *c = new int [nodelist.size()] ;
    
    for (int i=0; i< nodelist.size(); ++i){
        a[i] = nodelist[i].isLeaf;
        b[i] = a[i] * nodelist[i].size ;
    }
    scanExclusive(b, c, (int)nodelist.size() );
    for (int i=0; i< nodelist.size(); ++i){
        nodelist[i].triangleIndex = a[i] * c[i] ;
    }
    
    printf("siz:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",nodelist[i].size);
    }
    printf("\n");
    printf("b[]:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",b[i]);
    }
    printf("\n");
    printf("c[]:");
    for (int i=0; i< nodelist.size(); ++i){
        printf("  %02d",c[i]);
    }
    printf("\n");
    printf("tri:");
    for (int i=0; i< nodelist.size(); ++i){
        if(nodelist[i].isLeaf){
            printf("  %02d",nodelist[i].triangleIndex);
        }else{
            printf("  --");
        }
    }
    printf("\n");
    
    return ;
}


void buildSizes(vector<kdTreeNode> &nodelist, int level){
    
    for(int i=0; i<nodelist.size(); ++i){
        if(nodelist[i].level == level){
            if (!nodelist[i].isLeaf) {
                nodelist[i].size = nodelist[nodelist[i].leftAddress].size + nodelist[nodelist[i].leftAddress+1].size ;
            }else{
                nodelist[i].size = reduce(nodelist[i].triangleMask, nodelist[i].smallntris) + 1 ;
            }
        }
    }
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

