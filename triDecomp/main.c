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
#define ROOTPATH "/Users/Darren/Development"
#define MATBYTES 128
#define DUMP 1

typedef struct triangle{
    int        id ;
    SPVector   AA ;
    SPVector   BB ;
    SPVector   CC ;
    SPVector   NN ;
    char mat[MATBYTES] ;
} triangle;

typedef struct scatProps {
    char   matname[MATBYTES] ;
    float  corlen      ;
    float  roughness   ;
    float  resistivity ;
    float  specular    ;
    float  diffuse     ;
    float  shinyness   ;
} scatProps ;

#define NMATERIALS 9
scatProps materialProperties[NMATERIALS] = {
//   Name          corrLen      Roughness   Resistivity Specular    Diffuse     Shinyness
    {"ASPHALT",     0.5,        0.005,      1.0e18,     0.8,        0.2,        30.0        },
    {"BRICK",       0.1,        0.001,      1.0e18,     0.7,        0.3,        20.0        },
    {"CONCRETE",    0.3,        0.01,       120.0,      0.3,        0.7,        10.0        },
    {"METAL",       100.0,      0.0,        1.0e-8,     1.0,        0.0,        50.0        },
    {"MATERIAL",    100.0,      0.0,        0.0,        1.0,        0.0,        50.0        },
    {"ROOFING",     0.1,        0.1,        1.0e18,     0.6,        0.4,        40.0        },
    {"VEGETATION",  0.01,       0.1,        2000.0,     0.2,        0.8,        5.0         },
    {"WATER",       0.01,       0.1,        2.0e1,      1.0,        0.0,        50.0        },
    {"WOOD",        0.1,        0.001,      1.0e14,     0.6,        0.4,        10.0        }
} ;

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

int main(int argc, const char * argv[])
{
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, "triDecomp", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status) ;

    char *str = input_string((char *)"triangle filename", (char *)"filename",
                       (char *)"The name of a triangle file created with 'colladaToTriFile'",
                       (char *) ROOTPATH"/DATA/triangles.tri");
    
    FILE *fpin ;
    fpin = fopen(str, "r");
    if (fpin == NULL) {
        printf("Error : failed to open input triangle file %s\n",str);
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
    
    if (DUMP) {
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
    
    for (int itri=0; itri < ntri; itri++ ) {
        triangle tri = triArray[itri] ;
        char *mat = tri.mat ;
        int n=0;
        while (mat[n]) {
            mat[n] = toupper(mat[n]);
            n++;
        }
        
        for(int i=0; i<NMATERIALS; i++){
            if( strstr(mat, materialProperties[i].matname)!=NULL){
                printf(" Triangle %d (%s) is of type %s\n", itri, mat, materialProperties[i].matname);
                
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
                printf("%d Random points:\n",np);
                srand((unsigned int)time(NULL));
                float *pointList ;
                pointList = (float *)sp_malloc(sizeof(float) * (np+3) * 3);
                
                for (int n = 0; n<np; n++) {
                    double r1 =  ((double)(rand() % 65536)) / 65536.0 ;
                    double r2 =  ((double)(rand() % 65536)) / 65536.0 ;
                    SPVector V1,V2,V3,V ;
                    VECT_SCMULT(tri.AA, (1 - sqrt(r1)), V1);
                    VECT_SCMULT(tri.BB, (sqrt(r1) * (1-r2)), V2);
                    VECT_SCMULT(tri.CC, (sqrt(r1) * r2), V3);
                    VECT_ADD(V1, V2, V) ;
                    VECT_ADD(V, V3, V) ;
                    
                    printf("  %f,%f,%f\n",V.x,V.y,V.z);
                    
                    for(int p=0; p<3; p++){
                        pointList[n*3+p] = V.cell[p] ;
                    }
                }
                
                // Delauney triangulation of points
                //
                int nTriVerts ;
                pointList[np*3+0] = tri.AA.x ;
                pointList[np*3+1] = tri.AA.y ;
                pointList[np*3+2] = tri.AA.z ;
                pointList[np*3+3] = tri.BB.x ;
                pointList[np*3+4] = tri.BB.y ;
                pointList[np*3+5] = tri.BB.z ;
                pointList[np*3+6] = tri.CC.x ;
                pointList[np*3+7] = tri.CC.y ;
                pointList[np*3+8] = tri.CC.z ;
                
                for(int p=0; p<(np+3)*3; p++){
                    printf("Pointlist[%d] = %f\n",p,pointList[p]);
                }

                WORD * triangleIndexList ;
                
                triangleIndexList = BuildTriangleIndexList((void *)pointList, 10000.0f, np+3, 3, 0, &nTriVerts) ;
                
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
                
                for (int v=0; v<newNtri; v++) {
//                    printf(" Delauney Triangle %d\n",v);
                    printf(" %f,%f,%f\n",newTriArray[v].AA.x,newTriArray[v].AA.y,newTriArray[v].AA.z);
                    printf(" %f,%f,%f\n",newTriArray[v].BB.x,newTriArray[v].BB.y,newTriArray[v].BB.z);
                    printf(" %f,%f,%f\n",newTriArray[v].CC.x,newTriArray[v].CC.y,newTriArray[v].CC.z);
                    printf(" %f,%f,%f\n",newTriArray[v].AA.x,newTriArray[v].AA.y,newTriArray[v].AA.z);

                }
                
                // Adjust interior points by roughness
                
                // Calculate new triangle normals
                
                // Add each triangle to linked list of triangles
                
                // increase triangle count
                
            } else {
                // Add triangle to linked list of triangles
                
                // increase triangle count
            }
        }
    }
    
    // Write out triangles to file
    
    // close files
    
    return 0;
}

