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
#include "materialProperties.h"

void readKdTree(const char *filename,
                AABB * SceneBoundingBox,
                int * nTriangles,
                Triangle ** Triangles,
                int * nLeaves,
                int *** triangleLists,
                int  *nTreeNodes,
                KdData ** KdTree,
                TriCoords **tricos ) {
    
    FILE *fp;
    int i,trisInLeaf,tex;
    double d,nd_u,nd_v,k,kbu,kbv,kbd,kcu,kcv,kcd;
    double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz ;
    
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
    
    *Triangles = (Triangle *)sp_malloc(sizeof(Triangle) * *nTriangles) ;
    
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
        (*Triangles)[i].triNum = i;
        (*Triangles)[i].d      = d;
        (*Triangles)[i].nd_u   = nd_u;
        (*Triangles)[i].nd_v   = nd_v;
        (*Triangles)[i].k      = k;
        (*Triangles)[i].kbu    = kbu;
        (*Triangles)[i].kbv    = kbv;
        (*Triangles)[i].kbd    = kbd;
        (*Triangles)[i].kcu    = kcu;
        (*Triangles)[i].kcv    = kcv;
        (*Triangles)[i].kcd    = kcd;
        (*Triangles)[i].matInd = tex;
    }
    
//    if(fread(nTextures, sizeof(int), 1, fp)!=1){
//        printf("ERROR: Failed to read nTextures from file %s\n",filename);
//        exit(-1);
//    }
//    *textures = (Texture *)sp_malloc(sizeof(Texture) * *nTextures);
//    
//    int matID;
//    int matNameLen;
//    char *matName = "";
//    for (i=0; i < *nTextures; i++){
//        if(fread(&matID, sizeof(int), 1, fp)!=1){
//            printf("ERROR: Failed to read materialID texture param from file %s\n",filename);
//            exit(-1);
//        }
//        if(fread(&matNameLen, sizeof(int), 1, fp)!=1){
//            printf("ERROR: Failed to read material name length texture param from file %s\n",filename);
//            exit(-1);
//        }
//        matName = (char *)sp_malloc(sizeof (char)*matNameLen);
//        if(fread(matName, sizeof(char), matNameLen, fp)!=matNameLen){
//            printf("ERROR: Failed to read material name length texture param from file %s\n",filename);
//            exit(-1);
//        }
//        (*textures)[i].ka = 0.0 ;
//        (*textures)[i].kd = materialProperties[matID].diffuse ;
//        (*textures)[i].ks = materialProperties[matID].specular ;
//        (*textures)[i].n  = materialProperties[matID].shinyness ;
//    }
    
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
    
    *tricos = (TriCoords *)malloc(sizeof(TriCoords) * *nTriangles);
    if (*tricos == NULL) {
        printf("ERROR: Malloc failed for triangle coordinates. Requested %lu bytes\n",sizeof(tricos) * *nTriangles);
        exit(-1);
    }
    for (i=0; i< *nTriangles; i++){
        if(fread(&Ax, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read A for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&Ay, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read A for triangle %d from file %s\n",i,filename);
            exit(-1);
        }if(fread(&Az, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read A for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&Bx, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read B for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&By, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read B for triangle %d from file %s\n",i,filename);
            exit(-1);
        }if(fread(&Bz, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read B for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&Cx, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read C for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        if(fread(&Cy, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read C for triangle %d from file %s\n",i,filename);
            exit(-1);
        }if(fread(&Cz, sizeof(double), 1, fp)!=1){
            printf("ERROR: Failed to read C for triangle %d from file %s\n",i,filename);
            exit(-1);
        }
        (*tricos)[i].A.x  = Ax ;
        (*tricos)[i].A.y  = Ay ;
        (*tricos)[i].A.z  = Az ;
        (*tricos)[i].B.x  = Bx ;
        (*tricos)[i].B.y  = By ;
        (*tricos)[i].B.z  = Bz ;
        (*tricos)[i].Cc.x = Cx ;
        (*tricos)[i].Cc.y = Cy ;
        (*tricos)[i].Cc.z = Cz ;
    }

    fclose(fp);
    
    return;
    
}

