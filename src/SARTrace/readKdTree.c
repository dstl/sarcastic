//***************************************************************************
//
//  readKdTree.cpp
//  Sadilac
//
//  Created by Darren on 02/08/2012.
//  Copyright (c) 2012 [dstl]. All rights reserved.
//
//
// Description:
//
//
// CLASSIFICATION        :  UNCLASSIFIED
// Date of CLASSN        :  02/08/2012
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// THE SOFTWARE IN ITS ENTIRETY OR ANY SUBSTANTIAL PORTION SHOULD NOT BE
// USED AS A WHOLE OR COMPONENT PART OF DERIVATIVE SOFTWARE THAT IS BEING
// SOLD TO THE GOVERNMENT OF THE UNITED KINGDOM OF GREAT BRITAIN AND NORTHERN
// IRELAND.
//
//***************************************************************************

#include "readKdTree.h"
//#include "materialProperties.h"

void readKdTree(const char *filename,
                AABB * SceneBoundingBox,
                int * nTriangles,
                ATS ** accelTriangles,
                Triangle **triangles,
                int * nLeaves,
                int *** triangleLists,
                int  *nTreeNodes,
                KdData ** KdTree) {
    
    FILE *fp;
    int i,trisInLeaf,tex;
    double d,nd_u,nd_v,k,kbu,kbv,kbd,kcu,kcv,kcd;
    
    fp=fopen(filename, "r");
    if (fp==NULL){
        printf("ERROR: Failed to open file : %s\n",filename);
        exit(-1);
    }
    
    if(fread(&(SceneBoundingBox->AA.x), sizeof(double), 1, fp)!=1){
        printf("ERROR: Failed to read AAx from file %s\n",filename);
        exit(-1);
    }
    if(fread(&(SceneBoundingBox->AA.y), sizeof(double), 1, fp)!=1){
        printf("ERROR: Failed to read AAy from file %s\n",filename);
        exit(-1);
    }
    if(fread(&(SceneBoundingBox->AA.z), sizeof(double), 1, fp)!=1){
        printf("ERROR: Failed to read AAz from file %s\n",filename);
        exit(-1);
    }
    if(fread(&(SceneBoundingBox->BB.x), sizeof(double), 1, fp)!=1){
        printf("ERROR: Failed to read BBx from file %s\n",filename);
        exit(-1);
    }
    if(fread(&(SceneBoundingBox->BB.y), sizeof(double), 1, fp)!=1){
        printf("ERROR: Failed to read BBy from file %s\n",filename);
        exit(-1);
    }
    if(fread(&(SceneBoundingBox->BB.z), sizeof(double), 1, fp)!=1){
        printf("ERROR: Failed to read BBz from file %s\n",filename);
        exit(-1);
    }
    
    if(fread(nTriangles, sizeof(int), 1, fp)!=1){
        printf("ERROR: Failed to read number of Triangles from file %s\n",filename);
        exit(-1);
    }
    
    *accelTriangles = (ATS *)sp_malloc(sizeof(ATS) * *nTriangles) ;
    
    for (i=0; i< *nTriangles; i++){
        
        if(fread(&d, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read d for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&nd_u, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read nd_u for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&nd_v, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read nd_v for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&k, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read k for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&kbu, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read kbu for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&kbv, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read kbv for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&kbd, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read kbd for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&kcu, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read kcu for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&kcv, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read kcv for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&kcd, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read kcd for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&tex, sizeof(int), 1, fp)!=1){
            printf("ERROR: Failed to read tex for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*accelTriangles)[i].triNum = i;
        (*accelTriangles)[i].d      = d;
        (*accelTriangles)[i].nd_u   = nd_u;
        (*accelTriangles)[i].nd_v   = nd_v;
        (*accelTriangles)[i].k      = k;
        (*accelTriangles)[i].kbu    = kbu;
        (*accelTriangles)[i].kbv    = kbv;
        (*accelTriangles)[i].kbd    = kbd;
        (*accelTriangles)[i].kcu    = kcu;
        (*accelTriangles)[i].kcv    = kcv;
        (*accelTriangles)[i].kcd    = kcd;
        (*accelTriangles)[i].matInd = tex;
    }
    
    if(fread(nLeaves, sizeof(int), 1, fp)!=1){
        printf("ERROR: Failed to read nleaves from file %s\n",filename);
        exit(-1);
    }
    
    *triangleLists = (int **)malloc(sizeof(int *) * *nLeaves);
    if(*triangleLists == NULL){
        printf("ERROR: Failed to malloc triangleLists. Requested %lu bytes\n",sizeof(int *) * *nLeaves);
        exit(-1);
    }
    for(i = 0; i< *nLeaves; i++){
        if(fread(&trisInLeaf, sizeof(int), 1, fp)!=1){
            printf("ERROR: Failed to read number fo tris in leaf %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangleLists)[i] = (int *)malloc(sizeof(int) * (trisInLeaf+1));
        if((*triangleLists)[i] == NULL){
            printf("ERROR: Failed to malloc triangleLists[%d]. Requested %lu bytes\n",i,sizeof(int)* (trisInLeaf+1));
            exit(-1);
        }
        (*triangleLists)[i][0]=trisInLeaf;
        if(fread(&((*triangleLists)[i][1]), sizeof(int), trisInLeaf, fp)!=trisInLeaf){
            printf("ERROR: Failed to read triangle list %d from file %s\n",i,filename);
            exit(-1);
        }
    }
    
    // Now read in the kdTreeNodes
    if(fread(nTreeNodes, sizeof(int), 1, fp)!=1){
        printf("ERROR: Failed to read number of nodes from file %s\n",filename);
        exit(-1);
    }
    
    *KdTree = (KdData *)malloc(sizeof(KdData)* *nTreeNodes);
    if(*KdTree == NULL){
        printf("ERROR: Failed to malloc KdTree. Requested %lu bytes\n",sizeof(KdData)* *nTreeNodes);
        exit(-1);
    }
    
    size_t writtenKdDateSize = sizeof(KdData) - sizeof(AABB) - (sizeof(int)*6);
    for (i=0; i< *nTreeNodes; i++){
        
        if(fread(&((*KdTree)[i]), sizeof(writtenKdDateSize), 1, fp)!= 1){
            printf("ERROR: Failed to read KdData item %d from file %s\n",i,filename);
            exit(-1);
        }
    }
    
    *triangles = (Triangle *)sp_malloc(sizeof(Triangle) * *nTriangles);
    double dub;
    int ind;
    for (i=0; i< *nTriangles; i++){
        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read A for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].AA.x  = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read A for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].AA.y  = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read A for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].AA.z  = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read B for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].BB.x  = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read B for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].BB.y  = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read B for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].BB.z  = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read C for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].CC.x = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read C for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].CC.y = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read C for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].CC.z = dub ;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read N for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].NN.x = dub;

        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read N for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].NN.y = dub;
        
        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read N for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].NN.z = dub ;
        
        if(fread(&dub, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read Area for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].area = dub ;
        
        for(int j=0; j<9; j++){
            fread(&dub, sizeof(double), 1, fp);
            (*triangles)[i].globalToLocalMat[j] = dub;
        }
        for(int j=0; j<9; j++){
            fread(&dub, sizeof(double), 1, fp);
            (*triangles)[i].localToGlobalMat[j] = dub;
        }
        if(fread(&ind, sizeof(int), 1, fp)!=1){
            printf("ERROR: Failed to read material ID for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*triangles)[i].matId = ind;
    }

    fclose(fp);
    
    return;
    
}

void printKdTree(int nNodes, KdData *tree, int *triListData, int *triPtrs){
    
    KdData *node;
    
    for(int n=0; n<nNodes; n++){
        
        node = &(tree[n]) ;
        int isleaf = KDT_ISLEAF(node);
        
        printf("\n                           KdTree Node [%06d] \n",n);
        printf("                        -------------------------------\n");
        if (isleaf) {
            printf("Node is a leaf\n");
            printf("Node has the following ropes: [%06d], [%06d], [%06d], [%06d], [%06d], [%06d] \n",
                   node->leaf.Ropes[0], node->leaf.Ropes[1], node->leaf.Ropes[2], node->leaf.Ropes[3], node->leaf.Ropes[4], node->leaf.Ropes[5] );
            printf("Node is a leaf and contains the following triangles:  ");
    
            int leafIndex = KDT_OFFSET(node) ;
            int indexOfFirstForLeaf = triPtrs[leafIndex] ;
            int trisInLeaf = triListData[indexOfFirstForLeaf] ;
            
            for (int t=0; t<trisInLeaf; t++){
                int triIndex = triListData[indexOfFirstForLeaf+t];
                printf(" %d", triIndex);
            }
            printf("\n");
            
            printf("      AABB:\n");
            printf("   %f,%f,%f \n", node->leaf.aabb.AA.x, node->leaf.aabb.AA.y, node->leaf.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.BB.x, node->leaf.aabb.AA.y, node->leaf.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.BB.x, node->leaf.aabb.BB.y, node->leaf.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.AA.x, node->leaf.aabb.BB.y, node->leaf.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.AA.x, node->leaf.aabb.AA.y, node->leaf.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.AA.x, node->leaf.aabb.AA.y, node->leaf.aabb.BB.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.BB.x, node->leaf.aabb.AA.y, node->leaf.aabb.BB.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.BB.x, node->leaf.aabb.BB.y, node->leaf.aabb.BB.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.AA.x, node->leaf.aabb.BB.y, node->leaf.aabb.BB.z) ;
            printf("   %f,%f,%f \n", node->leaf.aabb.AA.x, node->leaf.aabb.AA.y, node->leaf.aabb.BB.z) ;
            printf("--------------------------------------\n");
        }else{
            printf("Node is a branch\n");
            printf("Split dimension is          : %d", KDT_DIMENSION(node));
            switch(KDT_DIMENSION(node)){
                    case 0: printf(" [x]\n"); break;
                    case 1: printf(" [y]\n"); break;
                    case 2: printf(" [z]\n"); break;
            }
            printf("split position is           : %f\n", node->branch.splitPosition);
            printf("First child node            : [%06d] \n", KDT_OFFSET(node));
            printf("Second child node           : [%06d] \n", KDT_OFFSET(node)+1);
            printf("      AABB:\n");
            printf("   %f,%f,%f \n", node->branch.aabb.AA.x, node->branch.aabb.AA.y, node->branch.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.BB.x, node->branch.aabb.AA.y, node->branch.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.BB.x, node->branch.aabb.BB.y, node->branch.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.AA.x, node->branch.aabb.BB.y, node->branch.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.AA.x, node->branch.aabb.AA.y, node->branch.aabb.AA.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.AA.x, node->branch.aabb.AA.y, node->branch.aabb.BB.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.BB.x, node->branch.aabb.AA.y, node->branch.aabb.BB.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.BB.x, node->branch.aabb.BB.y, node->branch.aabb.BB.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.AA.x, node->branch.aabb.BB.y, node->branch.aabb.BB.z) ;
            printf("   %f,%f,%f \n", node->branch.aabb.AA.x, node->branch.aabb.AA.y, node->branch.aabb.BB.z) ;
            printf("--------------------------------------\n");

            
        }
        
        
    }
}