/***************************************************************************
 *  main.c
 *  triDecomp
 *
 *  Created by Muff Darren on 25/09/2014.
 *  Copyright (c) 2014 [dstl]. All rights reserved.
 *
 *   Description:
 *      Programme to read in a triangle file and search for triangles
 *      made from materials with a diffuse scattering characteristic.
 *      Once found each triangle is decomposed into many smaller Delauney
 *      triangles whose size matches a correlation dimension associated 
 *      with the material type. Once decomposed, the triangle vertices are
 *      randomly aised or lowered above the triangle plane to provide a
 *      surface roughness associated with the material. Finally the new
 *      triangle file is written out.
 *
 *  CLASSIFICATION       :   UNCLASSIFIED
 *  Date of CLASSN       :   25/09/2014
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
 ***************************************************************************/

#include <stdio.h>
#include <SIlib.h>
#include <ctype.h>
#include <math.h>
#include "Clarkson-Delaunay.h"
#include "ranf.h"
#include "boxMullerRandom.h"
#include "materialProperties.h"
#define ROOTPATH "/Users/Darren/Development"
#define SHOWOUTPUT 1

typedef struct triangle{
    int        id ;
    SPVector   AA ;
    SPVector   BB ;
    SPVector   CC ;
    SPVector   NN ;
    char mat[MATBYTES] ;
} triangle;

struct triListNode{
    triangle tri   ;
    struct triListNode *next ;
};

int main(int argc, const char * argv[])
{
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, "triDecomp", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status) ;
    struct triListNode * triList ;
    struct triListNode * triAccessor ;
    triList       = sp_malloc(sizeof(struct triListNode)) ;
    triList->next = NULL ;
    triAccessor   = triList ;
    
    char *instr = input_string((char *)"Input triangle filename", (char *)"infilename",
                       (char *)"The name of a triangle file created with 'colladaToTriFile'",
                       (char *) ROOTPATH"/DATA/triangles.tri");
    
    FILE *fpin ;
    fpin = fopen(instr, "r");
    if (fpin == NULL) {
        printf("Error : failed to open input triangle file %s\n",instr);
        exit(1);
    }
    
    char *oustr = input_string((char *)"Output triangle filename", (char *)"oufilename",
                             (char *)"The name of the triangle file to be created",
                             (char *) ROOTPATH"/DATA/triangles_Delauney.tri");
    
    FILE *fpout ;
    fpout = fopen(oustr, "w");
    if (fpout == NULL) {
        printf("Error : failed to open input triangle file %s\n",oustr);
        exit(1);
    }
    
    int ntri;
    int vertexBytes;
    int materialBytes;
    fread(&ntri, sizeof(int), 1, fpin);
    fread(&vertexBytes,sizeof(int), 1, fpin);
    fread(&materialBytes, sizeof(int), 1, fpin);
    
    if (vertexBytes != sizeof(double) ) {
        printf("Error : Triangle vertices of size %d bytes not yet supported\n",vertexBytes);
        exit(1);
    }
    if (materialBytes != MATBYTES ) {
        printf("Error : Material Bytes in file does not match the program\n");
        exit(1);
    }
    
    triangle * triArray = (triangle *)sp_malloc(sizeof(triangle)*ntri);
    
    for (int itri=0; itri < ntri; itri++ ) {
        fread(&(triArray[itri].id), sizeof(int), 1, fpin) ;
        fread(&(triArray[itri].AA.x), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].AA.y), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].AA.z), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].BB.x), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].BB.y), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].BB.z), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].CC.x), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].CC.y), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].CC.z), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].NN.x), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].NN.y), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].NN.z), sizeof(double), 1, fpin) ;
        fread(&(triArray[itri].mat) , sizeof(char), materialBytes, fpin) ;
    }
    
    if (SHOWOUTPUT) {
        for (int itri=0; itri < ntri; itri++ ) {
            printf("Triangle %d\n", triArray[itri].id);
            printf("  AA       : %3.6f,%3.6f,%3.6f\n", triArray[itri].AA.x,triArray[itri].AA.y,triArray[itri].AA.z);
            printf("  BB       : %3.6f,%3.6f,%3.6f\n", triArray[itri].BB.x,triArray[itri].BB.y,triArray[itri].BB.z);
            printf("  CC       : %3.6f,%3.6f,%3.6f\n", triArray[itri].CC.x,triArray[itri].CC.y,triArray[itri].CC.z);
            printf("  Normal   : %3.6f,%3.6f,%3.6f\n", triArray[itri].NN.x,triArray[itri].NN.y,triArray[itri].NN.z);
            printf("  Material : %s \n",triArray[itri].mat);
        }
    }
    
    fclose(fpin);
    srand((unsigned int)time(NULL));
    
    int triDone ;
    int newTriCnt = 0;

    for (int itri=0; itri < ntri; itri++ ) {
        triangle tri = triArray[itri] ;
        triDone = 0;
        char *mat = tri.mat ;
        int n=0;
        while (mat[n]) {
            mat[n] = toupper(mat[n]);
            n++;
        }
        
        for(int i=0; i<NMATERIALS; i++){
            if( strstr(mat, materialProperties[i].matname)!=NULL && materialProperties[i].diffuse != 0.0){
                triDone = 1 ;
                // workout number of additional points
                //
                double corrPatch = SIPC_pi * materialProperties[i].corlen * materialProperties[i].corlen ;
                SPVector AB, AC, triCross ;
                VECT_SUB(tri.BB, tri.AA, AB);
                VECT_SUB(tri.CC, tri.AA, AC);
                VECT_CROSS(AB, AC, triCross);
                double triArea = VECT_MAG(triCross) / 2 ;
                int np = triArea / corrPatch ;
                
                // generate random points within triangle
                // The following is taken from
                // http://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
                // and from []Osada - "Shape DIstributions"
                // which provides:
                //          P = (1 - sqrt(r1)) * A + (sqrt(r1) * (1 - r2)) * B + (sqrt(r1) * r2) * C
                //
                if(SHOWOUTPUT){
                    printf("Triangle %d requires %d internal random points\n",itri,np);
                }
                float *pointList ;
                pointList = (float *)sp_malloc(sizeof(float) * (np+3) * 3);
                
                for (int n = 0; n<np; n++) {
//                    double r1 =  ((double)(rand() % 65536)) / 65536.0 ;
//                    double r2 =  ((double)(rand() % 65536)) / 65536.0 ;
                    double r1 = Ranf() ;
                    double r2 = Ranf() ;
                    SPVector V1,V2,V3,V ;
                    VECT_SCMULT(tri.AA, (1.0 - sqrt(r1)), V1);
                    VECT_SCMULT(tri.BB, (sqrt(r1) * (1-r2)), V2);
                    VECT_SCMULT(tri.CC, (sqrt(r1) * r2), V3);
                    VECT_ADD(V1, V2, V) ;
                    VECT_ADD(V, V3, V) ;
                    
                    if (SHOWOUTPUT) {
                        printf("  %f,%f,%f\n",V.x,V.y,V.z);
                    }
                    
                    for(int p=0; p<3; p++){
                        pointList[n*3+p] = V.cell[p] ;
                    }
                }
                
                // Add original triangle points to end of array
                //
                pointList[np*3+0] = tri.AA.x ;
                pointList[np*3+1] = tri.AA.y ;
                pointList[np*3+2] = tri.AA.z ;
                pointList[np*3+3] = tri.BB.x ;
                pointList[np*3+4] = tri.BB.y ;
                pointList[np*3+5] = tri.BB.z ;
                pointList[np*3+6] = tri.CC.x ;
                pointList[np*3+7] = tri.CC.y ;
                pointList[np*3+8] = tri.CC.z ;

                // Delauney triangulation of points
                //
                WORD * triangleIndexList ;
                int nTriVerts ;
                triangleIndexList = BuildTriangleIndexList((void *)pointList, 10000.0f, np+3, 3, 0, &nTriVerts) ;
                
                // Adjust interior points by roughness
                // Only do interior points
                //
                for (int n = 0; n<np; n++) {
                    SPVector N = tri.NN ;
                    SPVector V ;
//                    double r3 =  (((double)(rand() % 65536)) / 65536.0) - 0.5 ; // Random number between -0.5 & 0.5
                    float dist = box_muller(0.0, materialProperties[i].roughness) ;
                    
//                    float dist = materialProperties[i].roughness * r3 ;
                    VECT_SCMULT(N, dist, V);
                    pointList[n*3]   = pointList[n*3]   + V.x ;
                    pointList[n*3+1] = pointList[n*3+1] + V.y ;
                    pointList[n*3+2] = pointList[n*3+2] + V.z ;
                }

                int newNtri = nTriVerts / 3 ;
                triangle * newTriArray = (triangle *)sp_malloc(sizeof(triangle) * newNtri);
                for (int v=0; v<newNtri; v++) {
                    SPVector vecA,vecB,vecC ;
                    int pntA,pntB,pntC;
                    double Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz;
                    pntA = triangleIndexList[v*3];
                    pntB = triangleIndexList[v*3+1];
                    pntC = triangleIndexList[v*3+2];
                    Ax = pointList[pntA*3+0];
                    Ay = pointList[pntA*3+1];
                    Az = pointList[pntA*3+2];
                    Bx = pointList[pntB*3+0];
                    By = pointList[pntB*3+1];
                    Bz = pointList[pntB*3+2];
                    Cx = pointList[pntC*3+0];
                    Cy = pointList[pntC*3+1];
                    Cz = pointList[pntC*3+2];
                    VECT_CREATE(Ax, Ay, Az, vecA);
                    VECT_CREATE(Bx, By, Bz, vecB);
                    VECT_CREATE(Cx, Cy, Cz, vecC);
                    newTriArray[v].AA = vecA ;
                    newTriArray[v].BB = vecB ;
                    newTriArray[v].CC = vecC ;
                }
                free(triangleIndexList);
                free(pointList);

                if (SHOWOUTPUT) {
                    printf("Delauney Triangles for tringle %d are :\n",itri);
                    for (int v=0; v<newNtri; v++) {
                        printf("    %8.5f,%10.5f,%10.5f\n",newTriArray[v].AA.x,newTriArray[v].AA.y,newTriArray[v].AA.z);
                        printf("    %10.5f,%10.5f,%10.5f\n",newTriArray[v].BB.x,newTriArray[v].BB.y,newTriArray[v].BB.z);
                        printf("    %10.5f,%10.5f,%10.5f\n",newTriArray[v].CC.x,newTriArray[v].CC.y,newTriArray[v].CC.z);
                        printf("    %10.5f,%10.5f,%10.5f\n",newTriArray[v].AA.x,newTriArray[v].AA.y,newTriArray[v].AA.z);
                    }
                }
                
                
                // Calculate new triangle normals
                //
                for (int v=0; v<newNtri; v++) {
                    SPVector AB,AC, triCross;
                    VECT_SUB(newTriArray[v].BB, newTriArray[v].AA, AB);
                    VECT_SUB(newTriArray[v].CC, newTriArray[v].AA, AC);
                    VECT_CROSS(AB, AC, triCross);
                    VECT_NORM(triCross, newTriArray[v].NN);
                    if(VECT_DOT(tri.NN, newTriArray[v].NN) < 0 ){
                        VECT_SCMULT(newTriArray[v].NN, -1.0, newTriArray[v].NN);
                    }
                    strcpy(newTriArray[v].mat, materialProperties[i].matname) ;
                }
                
                // Write out triangles
                //
                for (int v=0; v<newNtri; v++) {
                    if ( triAccessor != NULL ){
                        while (triAccessor->next != NULL) {
                            triAccessor = triAccessor->next ;
                        }
                    }
                    triAccessor->next = sp_malloc(sizeof(struct triListNode));
                    triAccessor = triAccessor->next ;
                    if( triAccessor == NULL){
                        printf("Error out of memory for linked list\n");
                        exit(1);
                    }
                    triAccessor->next = NULL ;
                    triAccessor->tri  = newTriArray[v] ;
                    newTriCnt++;
                }
                
                // As we have found a material dont search through the material list for
                // any more
                //
                i = NMATERIALS-1;
            }
        }
        
        if (!triDone) {
            if ( triAccessor != NULL ){
                while (triAccessor->next != NULL) {
                    triAccessor = triAccessor->next ;
                }
            }
            triAccessor->next = sp_malloc(sizeof(struct triListNode));
            triAccessor = triAccessor->next ;
            if( triAccessor == NULL){
                printf("Error out of memory for linked list\n");
                exit(1);
            }
            triAccessor->next = NULL ;
            triAccessor->tri  = tri ;
            newTriCnt++;
        }
    }
    
    /*
     Triangle File Structure
     int number_of_triangles
     int sizeof_vertex_information
     int length_of_material_names
     foreach triangle{
     int triangle_number
     <sizeof_vertex_information> vertexA.x
     <sizeof_vertex_information> vertexA.y
     <sizeof_vertex_information> vertexA.z
     <sizeof_vertex_information> vertexB.x
     <sizeof_vertex_information> vertexB.y
     <sizeof_vertex_information> vertexB.z
     <sizeof_vertex_information> vertexC.x
     <sizeof_vertex_information> vertexC.y
     <sizeof_vertex_information> vertexC.z
     <sizeof_vertex_information> normal.x
     <sizeof_vertex_information> normal.y
     <sizeof_vertex_information> normal.z
     char[<length_of_material_names>] Material_name
     }
     */

    fwrite(&newTriCnt, sizeof(int),     1, fpout);
    fwrite(&vertexBytes,sizeof(int),    1, fpout);
    fwrite(&materialBytes, sizeof(int), 1, fpout);
    triAccessor = triList ;

    int itri = 0;
    triangle tri;
    while (triAccessor->next != NULL) {
        triAccessor = triAccessor->next ;
        tri = triAccessor->tri ;
        fwrite(&(itri), sizeof(int), 1, fpout) ;
        fwrite(&(tri.AA.x), sizeof(double), 1, fpout) ;
        fwrite(&(tri.AA.y), sizeof(double), 1, fpout) ;
        fwrite(&(tri.AA.z), sizeof(double), 1, fpout) ;
        fwrite(&(tri.BB.x), sizeof(double), 1, fpout) ;
        fwrite(&(tri.BB.y), sizeof(double), 1, fpout) ;
        fwrite(&(tri.BB.z), sizeof(double), 1, fpout) ;
        fwrite(&(tri.CC.x), sizeof(double), 1, fpout) ;
        fwrite(&(tri.CC.y), sizeof(double), 1, fpout) ;
        fwrite(&(tri.CC.z), sizeof(double), 1, fpout) ;
        fwrite(tri.mat, sizeof(char), materialBytes, fpout);
        itri++;
    }
    printf("Written %d triangles\n",itri);
    
    triAccessor = triList ;
    while( triAccessor->next != NULL ) {
        struct triListNode *node ;
        node = triAccessor ;
        triAccessor = triAccessor->next;
        free(node);
    }
    fclose(fpout);
    
    return 0;
}

