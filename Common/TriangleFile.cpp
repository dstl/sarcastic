/***************************************************************************
 *
 *       Module:    TriangleFile.cpp
 *      Program:    SARCASTIC
 *   Created by:    Darren Muff on 21/07/2015.
 *                  Copyright (c) 2015 Dstl. All rights reserved.
 *
 *   Description:
 *      Class header for TriangleFile class
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  21/07/2015
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

#include "TriangleFile.h"
#include "materialProperties.h"
#include <string.h>

std::vector<SPVector> uniqueVertices(std::vector<Triangle> triangles);
bool SPVectorsAreEqual(SPVector a, SPVector b);
FILE * initialise_ply_file(const char *fname, int nvertices, int ntris);

TriangleFile::TriangleFile( std::vector<Triangle> tris ) {
    triangles = tris ;
}

void TriangleFile::WriteFile(std::string fname){
    FILE *fp;
    fp = fopen(fname.c_str(), "w");
    if (fp==NULL) {
        printf("Error : Could not open file %s\n",fname.c_str());
        exit(1);
    }
    
    int ntri = (int)triangles.size() ;
    int vertexBytes = sizeof(double) ;
    int MAXMATNAME = MATBYTES ; // bytes
    
    fwrite(&ntri, sizeof(int), 1, fp);
    fwrite(&vertexBytes, sizeof(int),1,fp);
    fwrite(&MAXMATNAME, sizeof(int), 1, fp);
    
    for (int i=0; i<ntri; i++){
        
        fwrite(&i, sizeof(int), 1, fp);
        double x,y,z, *p ;
        
        Triangle tri = triangles[i];
        x = tri.AA.x;
        y = tri.AA.y;
        z = tri.AA.z;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&z, sizeof(double), 1, fp);
        
        x = tri.BB.x;
        y = tri.BB.y;
        z = tri.BB.z;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&z, sizeof(double), 1, fp);
        
        x = tri.CC.x;
        y = tri.CC.y;
        z = tri.CC.z;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&z, sizeof(double), 1, fp);
        
        x = tri.NN.x;
        y = tri.NN.y;
        z = tri.NN.z;
        fwrite(&x, sizeof(double), 1, fp);
        fwrite(&y, sizeof(double), 1, fp);
        fwrite(&z, sizeof(double), 1, fp);
        
        x = tri.area ;
        fwrite(&x, sizeof(double), 1, fp);
        
        p = tri.globalToLocalMat;
        for (int j=0; j<9; j++){
            x = p[j] ;
            fwrite(&x, sizeof(double), 1, fp);
        }
        p = tri.localToGlobalMat ;
        for (int j=0; j<9; j++){
            x = p[j] ;
            fwrite(&x, sizeof(double), 1, fp);
        }
        char matname[MAXMATNAME];
        strcpy(matname, materialProperties[tri.matId].matname);
        fwrite(matname, sizeof(char), MAXMATNAME, fp);
        
    }
    fclose(fp);
    return ;
}

void TriangleFile::WritePLYFile( std::string fname, bool binary) {
        
    std::vector<SPVector> v;
    
    v = uniqueVertices(triangles);
    
    FILE * fp = initialise_ply_file(fname.c_str(), (int)v.size(), (int)triangles.size());
    
    // rebuild triangles to use unique vertices
    //
    int t[(int)triangles.size()][3];
    for(int i=0; i<triangles.size(); i++){
        SPVector aa,bb,cc;
        aa = triangles[i].AA ;
        bb = triangles[i].BB ;
        cc = triangles[i].CC ;
        
        for(int j=0; j<v.size(); j++){
            if(SPVectorsAreEqual(v[j], aa))t[i][0] = j ;
            if(SPVectorsAreEqual(v[j], bb))t[i][1] = j ;
            if(SPVectorsAreEqual(v[j], cc))t[i][2] = j ;
        }
    }
    
    // Write out the unique vertices
    //
    for (int i=0; i<v.size(); i++){
        fprintf(fp,"%4.4f %4.4f %4.4f\n",v[i].x,v[i].y,v[i].z);
    }
    // now write out triangles in 'face' element
    //
    for(int i=0; i<triangles.size(); i++){
        fprintf(fp, "3 %d %d %d %d %d %d\n",t[i][0], t[i][1], t[i][2],
                materialColours[triangles[i].matId][0],materialColours[triangles[i].matId][1],materialColours[triangles[i].matId][2]);
    }
    
    fclose(fp);
    return ;


}

FILE * initialise_ply_file(const char *fname, int nvertices, int ntris){
    
    FILE *fp = fopen(fname, "w");
    
    if (fp == NULL) {
        printf("Error : could not open file for writing\n");
        exit(1);
    }
    
    fprintf(fp,"ply\n");
    fprintf(fp,"format ascii 1.0\n");
    fprintf(fp,"comment Triangle PLY File by Sarcastic\n");
    fprintf(fp,"element vertex %d\n",nvertices);
    fprintf(fp,"property float x\n");
    fprintf(fp,"property float y\n");
    fprintf(fp,"property float z\n");
    fprintf(fp,"element face %d\n",ntris);
    fprintf(fp,"property list uchar int vertex_index\n");
    fprintf(fp,"property uchar red\n");
    fprintf(fp,"property uchar green\n");
    fprintf(fp,"property uchar blue\n");
    fprintf(fp,"end_header\n");
    
    return(fp);
}

std::vector<SPVector> uniqueVertices(std::vector<Triangle> triangles){
    int ntris = (int)triangles.size();
    SPVector test;
    bool repeated;
    
    std::vector<SPVector> vertices_orig;
    std::vector<SPVector> vertices_new;
    
    for(int i=0; i<ntris; i++){
        vertices_orig.push_back(triangles[i].AA);
        vertices_orig.push_back(triangles[i].BB);
        vertices_orig.push_back(triangles[i].CC);
    }
    
    for(int i=0; i<vertices_orig.size(); i++){
        test     = vertices_orig[i];
        repeated = false;
        for (int j=0; j<vertices_new.size(); j++){
            if( SPVectorsAreEqual(test, vertices_new[j])){
                repeated = true;
            }
        }
        if (repeated == false) {
            vertices_new.push_back(test);
        }
    }
    return vertices_new ;
}

bool SPVectorsAreEqual(SPVector a, SPVector b){
    return (a.x==b.x && a.y==b.y && a.z==b.z) ;
}
