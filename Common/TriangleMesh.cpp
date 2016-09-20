/***************************************************************************
 *
 *       Module:    TriangleMesh.cpp
 *      Program:    SARCASTIC
 *   Created by:    Darren Muff on 18/09/2016.
 *                  Copyright (c) 2016 Dstl. All rights reserved.
 *
 *   Description:
 *      Class header for TriangleMesh class
 *
 *   CLASSIFICATION        :  UNCLASSIFIED
 *   Date of CLASSN        :  18/09/2016
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

#include "TriangleMesh.hpp"
#include "colladainterface.h"
#include <stdio.h>
#include <fstream>

#define FILE_MAX_LINE_LENGTH (256)

extern "C" {
#include "matrixMultiplication.h"
}

void TriangleMesh::readDAEFile  ( std::string filename ){
    std::vector<ColGeom> geom_vec;    // Vector containing COLLADA meshes
    int num_objects ;
    
    // Geom_vec is a C++ Vector containing an array of 'meshes' from the collada file
    // We are only interested in the triangles in the file so loop though the
    // Geom_vec array, find all the triangles then write them to a triangle vector
    //
    ColladaInterface::readGeometries(&geom_vec, filename.c_str());
    
    std::vector<Triangle> tri_vec;
    num_objects = (int) geom_vec.size();
    for(int i=0; i<num_objects; i++){
        if(geom_vec[i].primitive == GL_TRIANGLES ){
            float *f      = (float *)geom_vec[i].map["POSITION"].data;
            int stride    = geom_vec[i].map["POSITION"].stride ;
            float scaling = geom_vec[i].scaling;
            double localPt[4], globalPt[4];
            
            for (int j=0; j<geom_vec[i].index_count/3; j++){
                short AAi,BBi,CCi;
                AAi = geom_vec[i].indices[j*3+0] ;
                BBi = geom_vec[i].indices[j*3+1] ;
                CCi = geom_vec[i].indices[j*3+2] ;
                float AAx,AAy,AAz, BBx,BBy,BBz, CCx,CCy,CCz;
                
                localPt[0] = f[AAi*stride+0] * scaling;
                localPt[1] = f[AAi*stride+1] * scaling;
                localPt[2] = f[AAi*stride+2] * scaling;
                localPt[3] = scaling ;
                matmul(geom_vec[i].transform, localPt, globalPt, 4, 4, 1, 4);
                AAx = globalPt[0];
                AAy = globalPt[1];
                AAz = globalPt[2];
                
                localPt[0] = f[BBi*stride+0] * scaling;
                localPt[1] = f[BBi*stride+1] * scaling;
                localPt[2] = f[BBi*stride+2] * scaling;
                localPt[3] = scaling ;
                matmul(geom_vec[i].transform, localPt, globalPt, 4, 4, 1, 4);
                BBx = globalPt[0];
                BBy = globalPt[1];
                BBz = globalPt[2];
                
                localPt[0] = f[CCi*stride+0] * scaling;
                localPt[1] = f[CCi*stride+1] * scaling;
                localPt[2] = f[CCi*stride+2] * scaling;
                localPt[3] = scaling ;
                matmul(geom_vec[i].transform, localPt, globalPt, 4, 4, 1, 4);
                CCx = globalPt[0];
                CCy = globalPt[1];
                CCz = globalPt[2];
                
                SPVector AA,BB,CC;
                VECT_CREATE(AAx, AAy, AAz, AA);
                VECT_CREATE(BBx, BBy, BBz, BB);
                VECT_CREATE(CCx, CCy, CCz, CC);
                
                std::string triMaterial ;
                if (geom_vec[i].materialSide1 == std::string("material") || geom_vec[i].materialSide1 == std::string("")) {
                    if (geom_vec[i].materialSide2 != std::string("material") && geom_vec[i].materialSide2 != std::string("")) {
                        triMaterial = geom_vec[i].materialSide2;
                    }else{
                        triMaterial = materialProperties[0].matname;
                    }
                }else if (geom_vec[i].materialSide2 == std::string("material") || geom_vec[i].materialSide2 == std::string("")){
                    if (geom_vec[i].materialSide1 != std::string("material") && geom_vec[i].materialSide1 != std::string("")) {
                        triMaterial = geom_vec[i].materialSide1;
                    }else{
                        triMaterial = materialProperties[0].matname;
                    }
                }
                
                rawTri t = rawTri(AA, BB, CC,triMaterial);
                if(t.A != 0.0){
                    printf("A: %5.2f,%5.2f,%5.2f B: %5.2f,%5.2f,%5.2f C: %5.2f,%5.2f.%5.2f N: %5.2f,%5.2f.%5.2f - %s\n",AA.x,AA.y,AA.z,BB.x,BB.y,BB.z,CC.x,CC.y,CC.z, t.N.x,t.N.y,t.N.z,triMaterial.c_str());
                    rawTriBuff.push_back(t) ;
                }else{
                    printf("ERROR area is %f: A: %f,%f,%f  B: %f,%f,%f  C: %f,%f.%f - %s\n",t.A,AA.x,AA.y,AA.z,BB.x,BB.y,BB.z,CC.x,CC.y,CC.z,triMaterial.c_str());
                    
                }
            }
        }
    }
    
    sortTrianglesAndPoints() ;
    ColladaInterface::freeGeometries(&geom_vec) ;
    
    return ;
}

void TriangleMesh::readPLYFile  ( std::string filename ){
    
    std::ifstream file(filename) ;
    std::string line ;
    std::string name1 ;
    std::string name2 ;
    std::string name3 ;
    
    double x,y,z;
    int t[8];
    int nvertices = 0;
    int ntriangles = 0;

    do {
        std::getline(file, line) ;
        std::stringstream linestream(line) ;
        linestream >> name1 >> name2 >> name3 ;
        
        if(!strcmp(name1.c_str(), "element")){
            if (!strcmp(name2.c_str(), "vertex")) {
                nvertices = std::stoi(name3);
            }else if( !strcmp(name2.c_str(), "face") ){
                ntriangles = std::stoi(name3);
            }
        }
        
    } while (strcmp( name1.c_str(), "end_header") ) ;
    
    if (nvertices == 0 || ntriangles == 0) {
        printf("Error: failed to read number of triangles / vertices from .PLY header\n");
        exit(1);
    }
    
    vertices.reserve(nvertices) ;
    
    for(int i=0; i<nvertices; ++i){
        std::getline(file, line) ;
        std::stringstream linestream(line) ;
        linestream >> x >> y >> z ;
        Triangle3DVec vec(x,y,z);
        vertices.push_back(vec) ;
    }
    
    triangles.reserve(ntriangles) ;
    
    for(int i=0; i<ntriangles; ++i){
        std::getline(file, line) ;
        std::stringstream linestream(line) ;
        linestream >> t[0] >> t[1] >> t[2] >> t[3] >> t[4] >> t[5] >> t[6] >> t[7] ;
        Triangle tri = Triangle(t[1], t[2], t[3], t[7]) ;
        triangles.push_back(tri);
    }
    
    file.close() ;
    return ;
    
}

void TriangleMesh::writePLYFile ( std::string filename ){
    
    if (!sorted) {
        sortTrianglesAndPoints() ;
    }
    
    FILE *fp = fopen(filename.c_str(), "w");
    
    if (fp == NULL) {
        printf("Error : could not open file %s for writing\n",filename.c_str());
        exit(1);
    }
    
    fprintf(fp,"ply\n");
    fprintf(fp,"format ascii 1.0\n");
    fprintf(fp,"comment Triangle PLY File by Sarcastic\n");
    fprintf(fp,"element vertex %d\n",(int)vertices.size());
    fprintf(fp,"property float x\n");
    fprintf(fp,"property float y\n");
    fprintf(fp,"property float z\n");
    fprintf(fp,"element face %d\n",(int)triangles.size());
    fprintf(fp,"property list uchar int vertex_index\n");
    fprintf(fp,"property uchar red\n");
    fprintf(fp,"property uchar green\n");
    fprintf(fp,"property uchar blue\n");
    fprintf(fp,"property uchar mat\n");
    fprintf(fp,"end_header\n");
    
    // Write out the vertices
    //
    for (int i=0; i<vertices.size(); i++){
        fprintf(fp,"%4.4f %4.4f %4.4f\n",vertices[i].x,vertices[i].y,vertices[i].z);
    }
    // now write out triangles in 'face' element
    //
    for(int i=0; i<triangles.size(); i++){
        fprintf(fp, "3 %d %d %d %d %d %d %d\n",triangles[i].a, triangles[i].b, triangles[i].c,
                materialColours[triangles[i].mat][0],materialColours[triangles[i].mat][1],materialColours[triangles[i].mat][2], triangles[i].mat);
    }
    
    fclose(fp);
    return ;
}

void TriangleMesh::sortTrianglesAndPoints(){
    
    Triangle3DVec p;
    Triangle tr;
    int idxA, idxB, idxC;
    
    int npoints = (int)rawTriBuff.size() * 3 ;
    
    vertices.reserve(npoints) ;
    triangles.reserve(rawTriBuff.size()) ;
    
    // Create a vector of triangle vertices
    //
    for(int t=0; t<rawTriBuff.size(); ++t){
        p = Triangle3DVec(rawTriBuff[t].aa.x, rawTriBuff[t].aa.y, rawTriBuff[t].aa.z);
        vertices.push_back(p) ;
        p = Triangle3DVec(rawTriBuff[t].bb.x, rawTriBuff[t].bb.y, rawTriBuff[t].bb.z);
        vertices.push_back(p) ;
        p = Triangle3DVec(rawTriBuff[t].cc.x, rawTriBuff[t].cc.y, rawTriBuff[t].cc.z);
        vertices.push_back(p) ;
    }
    
    // sort the vertices into unique points
    //
    sort(vertices.begin(), vertices.end());
    vertices.erase( unique(vertices.begin(), vertices.end()), vertices.end()) ;
    
    //rebuild triangles using the index of the point in the points vector
    //
    for(int t=0; t<rawTriBuff.size(); ++t){
        p = Triangle3DVec(rawTriBuff[t].aa.x, rawTriBuff[t].aa.y, rawTriBuff[t].aa.z);
        idxA = (int) (lower_bound(vertices.begin(), vertices.end(), p) - vertices.begin()) ;
        p = Triangle3DVec(rawTriBuff[t].bb.x, rawTriBuff[t].bb.y, rawTriBuff[t].bb.z);
        idxB = (int) (lower_bound(vertices.begin(), vertices.end(), p) - vertices.begin()) ;
        p = Triangle3DVec(rawTriBuff[t].cc.x, rawTriBuff[t].cc.y, rawTriBuff[t].cc.z);
        idxC = (int) (lower_bound(vertices.begin(), vertices.end(), p) - vertices.begin()) ;
        
        p = Triangle3DVec(rawTriBuff[t].N.x, rawTriBuff[t].N.y, rawTriBuff[t].N.z);
        tr = Triangle(idxA, idxB, idxC, rawTriBuff[t].mat, rawTriBuff[t].N);
        tr.A = rawTriBuff[t].A ;
        
        triangles.push_back(tr);
    }
    
    // sort the triangle references vector and remove any duplicate triangles
    //
    sort(triangles.begin(), triangles.end());
    triangles.erase( unique(triangles.begin(), triangles.end()), triangles.end()) ;
    
    sorted = true;
    return ;
    
}
