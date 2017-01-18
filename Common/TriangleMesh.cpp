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
        addTriangle(t[1], t[2], t[3], t[7]);
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


void TriangleMesh::buildRawTris(){
    if(rawTriBuff.size() !=0){
        rawTriBuff.clear() ;
    }
    for (int i=0; i<triangles.size(); ++i){
        rawTri rt = rawTri(vertices[triangles[i].a].asSPVector(),
                           vertices[triangles[i].b].asSPVector(),
                           vertices[triangles[i].c].asSPVector(),
                           triangles[i].mat) ;
        rawTriBuff.push_back(rt) ;
    }
    return ;
}

void TriangleMesh::sortTrianglesAndPoints(){
    
    // This only works if there are raw triangles in the triangle buffer.
    // Otherwise it can delete vertices and so the triangle indices are wrong
    //
    
    Triangle3DVec p;
    Triangle tr;
    int idxA, idxB, idxC;
    
    if(sorted){
        return;
    }
    if(rawTriBuff.size() == 0){
        buildRawTris();
    }
    vertices.clear();
    triangles.clear();
    
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
        
        tr = Triangle(idxA, idxB, idxC, rawTriBuff[t].mat, rawTriBuff[t].N, rawTriBuff[t].d);
        tr.Area = rawTriBuff[t].A ;
        
        triangles.push_back(tr);
    }
    
    // sort the triangle references vector and remove any duplicate triangles
    //
    sort(triangles.begin(), triangles.end());
    triangles.erase( unique(triangles.begin(), triangles.end()), triangles.end()) ;
    
    checkIntegrityAndRepair() ;
    
    sorted = true;
    return ;
    
}

void nad(SPVector aa, SPVector bb, SPVector cc, SPVector *normal, float *area, float *distance)
// helper function to take coordinates of a triangle and calculate the
// Normal, Area and Distance (of the triangle plane to the origin) of it
//
{
    SPVector N, ab, bc, x;
    double l;
    VECT_SUB(bb, aa, ab);
    VECT_SUB(cc, bb, bc);
    VECT_CROSS(ab, bc, x);
    l = VECT_MAG(x);
    *area = l * 0.5;
    VECT_SCMULT(x, (1/l), N);
    *distance = VECT_DOT(aa, N) ;
    *normal = N ;
    return ;
}

void TriangleMesh::addTriangle(int a, int b, int c)
{
    if( a > vertices.size() || b > vertices.size() || c > vertices.size()){
        printf("Possible error: Triangle has indices larger than number of vertices in mesh\n");
    }
    addTriangle(a, b, c, 0) ;
    return ;
}

void TriangleMesh::addTriangle(Triangle tri){
    addTriangle(tri.a, tri.b, tri.c, tri.mat);
    return ;
}

void TriangleMesh::addTriangle(int a, int b, int c, int mat)
{
    SPVector aa, bb, cc, NN ;
    float area, distance ;
    if( a > vertices.size() || b > vertices.size() || c > vertices.size()){
        printf("Possible error: Triangle has indices larger than number of vertices in mesh - not adding it\n");
        return ;
    }else{
        VECT_CREATE(vertices[a].x, vertices[a].y, vertices[a].z, aa);
        VECT_CREATE(vertices[b].x, vertices[b].y, vertices[b].z, bb);
        VECT_CREATE(vertices[c].x, vertices[c].y, vertices[c].z, cc);
        nad(aa, bb, cc, &NN, &area, &distance) ;
        Triangle t = Triangle(a, b, c, mat, NN, distance) ;
        t.Area = area ;
        triangles.push_back(t) ;
        sorted = false ;
        return ;
    }
}

void TriangleMesh::addTriangle(rawTri tri){
    SPVector NN;
    float area, distance ;
    int s ;
    
    rawTriBuff.push_back(tri) ;
    s = (int)vertices.size() ;
    vertices.push_back(tri.aa);
    vertices.push_back(tri.bb);
    vertices.push_back(tri.cc);
    nad(tri.aa, tri.bb, tri.cc, &NN, &area, &distance) ;
    Triangle t = Triangle(s, s+1, s+2, tri.mat, NN, distance) ;
    triangles.push_back(t) ;
    
    sorted = false ;
    return ;
}

double pround(double x, int precision)
{
    if (x == 0.)
        return x;
    double a = pow(10,precision);
    double b = x * a ;
    int c = (int)floor(b+0.5);
    double ans = ((double)c) / a ;
    return ans;
}

void TriangleMesh::addTriangle(SPVector AA, SPVector BB, SPVector CC, int mat)
{
    SPVector NN, AAx,BBx,CCx;
    float area, distance ;
    int s ;
    
    AAx.x = pround(AA.x, 6);
    AAx.y = pround(AA.y, 6);
    AAx.z = pround(AA.z, 6);
    BBx.x = pround(BB.x, 6);
    BBx.y = pround(BB.y, 6);
    BBx.z = pround(BB.z, 6);
    CCx.x = pround(CC.x, 6);
    CCx.y = pround(CC.y, 6);
    CCx.z = pround(CC.z, 6);
    
    
    rawTri tri(AAx,BBx,CCx,mat);
    rawTriBuff.push_back(tri);
    s = (int)vertices.size() ;
    vertices.push_back(tri.aa);
    vertices.push_back(tri.bb);
    vertices.push_back(tri.cc);
    nad(tri.aa, tri.bb, tri.cc, &NN, &area, &distance) ;
    Triangle t = Triangle(s, s+1, s+2, tri.mat, NN, distance) ;
    triangles.push_back(t) ;
    
    sorted = false ;
    return ;
}

void TriangleMesh::buildHalfEdges(){
    halfEdge heA, heB, heC;
    
    //    halfedges.reserve(triangles.size()*3) ;
    
    for(int t=0; t<triangles.size(); ++t){
        Triangle tri = triangles[t] ;
        
        heA.vertex = tri.a ;
        heA.face = t ;
        heA.nextVertex = tri.b ;
        halfedges.push_back(heA) ;
        
        heB.vertex = tri.b ;
        heB.face = t ;
        heB.nextVertex = tri.c ;
        halfedges.push_back(heB) ;
        
        heC.vertex = tri.c ;
        heC.face = t ;
        heC.nextVertex = tri.a ;
        halfedges.push_back(heC) ;
    }
    
    for(int it = 0; it < halfedges.size()-1; ++it){
        halfedges[it].nextHalfEdge = it+1 ;
    }
    
    // For each halfedge find its opposite side
    //
    halfEdge h;
    for(int i = 0; i < halfedges.size(); ++i){
        h.vertex = halfedges[i].nextVertex ;
        h.nextVertex = halfedges[i].vertex ;
        auto v = std::find(halfedges.begin(), halfedges.end(), h);
        if(v != halfedges.end()){
            halfedges[i].oppositeHalfEdge = (int)(v-halfedges.begin()) ;
        }
    }
    
    halfEdgesBuilt = true ;
    return ;
    
}

std::vector<halfEdge> TriangleMesh::edges(){
    if(!halfEdgesBuilt){
        buildHalfEdges() ;
    }
    
    std::vector<halfEdge> matches;
    halfEdge h;
    for(int i = 0; i < halfedges.size(); ++i){
        h = halfedges[i] ;
        if (h.oppositeHalfEdge == -1) {
            matches.push_back(h) ;
        }
    }
    
    return(matches);
}

rawTri TriangleMesh::asRawTriangle(long int triangleIndex){
    if(triangleIndex > triangles.size()){
        printf("Error : requested a triangle that does not exist in mesh\n");
        exit(1);
    }
    rawTri t = rawTri(vertices[triangles[triangleIndex].a].asSPVector(),
                      vertices[triangles[triangleIndex].b].asSPVector(),
                      vertices[triangles[triangleIndex].c].asSPVector(),
                      triangles[triangleIndex].mat);
    return t ;
}
void TriangleMesh::printTriangles()
{
    for(auto it=rawTriBuff.begin(); it!=rawTriBuff.end(); ++it){
        it->print();
    }
    return ;
}

void TriangleMesh::checkIntegrityAndRepair(){
    int aa, bb, cc;
    SPVector A,B,Cc,Nclcv;
    Triangle3DVec Ntri, Ncalc;
    float area,distance;
    std::vector<Triangle> newtriangles;
    
    for (int i=0 ; i< triangles.size(); ++i){
        aa = triangles[i].a ;
        bb = triangles[i].b ;
        cc = triangles[i].c ;
        
        A = vertices[aa].asSPVector() ;
        B = vertices[bb].asSPVector() ;
        Cc = vertices[cc].asSPVector() ;
        Ntri = triangles[i].N ;
        
        nad(A, B, Cc, &Nclcv, &area, &distance);
        Ncalc = Triangle3DVec(Nclcv) ;
        
        if (!(Ncalc == Ntri) || area == 0.0 || isnan(distance) ) {
            printf("Mesh integrity check failed for triangle %d\n",i);
            printf("(Naughty triangle has Normal: %f,%f,%f, Planar Distace: %f, Area: %f)\n",
                   Nclcv.x, Nclcv.y, Nclcv.z, distance, area) ;
            printf("Repairing ....\n");
            triangles.erase(triangles.begin()+i);
            --i;
        }
    }
    return ;
}

void TriangleMesh::buildTriangleAABBs()
// function that calculates the axis-aligned bounding box (AABB) for each triangle and
// stores it in a vector
//
{
    SPVector aa, bb, cc;
    SPVector min, max;
    
    AABBs.clear();
    
    for(int i=0; i<triangles.size(); ++i){
        
        VECT_CREATE( 9e9,  9e9,  9e9, min);
        VECT_CREATE(-9e9, -9e9, -9e9, max);
        aa = vertices[triangles[i].a].asSPVector() ;
        bb = vertices[triangles[i].b].asSPVector() ;
        cc = vertices[triangles[i].c].asSPVector() ;
        min.x = (aa.x < min.x) ? aa.x : min.x ;
        min.x = (bb.x < min.x) ? bb.x : min.x ;
        min.x = (cc.x < min.x) ? cc.x : min.x ;
        min.y = (aa.y < min.y) ? aa.y : min.y ;
        min.y = (bb.y < min.y) ? bb.y : min.y ;
        min.y = (cc.y < min.y) ? cc.y : min.y ;
        min.z = (aa.z < min.z) ? aa.z : min.z ;
        min.z = (bb.z < min.z) ? bb.z : min.z ;
        min.z = (cc.z < min.z) ? cc.z : min.z ;
        
        max.x = (aa.x > max.x) ? aa.x : max.x ;
        max.x = (bb.x > max.x) ? bb.x : max.x ;
        max.x = (cc.x > max.x) ? cc.x : max.x ;
        max.y = (aa.y > max.y) ? aa.y : max.y ;
        max.y = (bb.y > max.y) ? bb.y : max.y ;
        max.y = (cc.y > max.y) ? cc.y : max.y ;
        max.z = (aa.z > max.z) ? aa.z : max.z ;
        max.z = (bb.z > max.z) ? bb.z : max.z ;
        max.z = (cc.z > max.z) ? cc.z : max.z ;
        
        AABB aabb = AABB(min,max) ;
        AABBs.push_back(aabb);
    }
    return ;
}

/*void TriangleMesh::buildTriangleAABBs(int dim, float pos)
// function for calculating the AAB for each triangle in the mesh and clips it against
// a split position
//
{
    SPVector aa, bb, cc;
    SPVector min, max;
    
    AABBs.clear() ;
    
    for(int i=0; i<triangles.size(); ++i){
        VECT_CREATE( 9e9,  9e9,  9e9, min);
        VECT_CREATE(-9e9, -9e9, -9e9, max);
        aa = vertices[triangles[i].a].asSPVector() ;
        bb = vertices[triangles[i].b].asSPVector() ;
        cc = vertices[triangles[i].c].asSPVector() ;
        
        std::vector<double> xs;
        std::vector<double> ys;
        std::vector<double> zs;
        
        int k = dim;
        xs.push_back(aa.cell[k]);
        xs.push_back(bb.cell[k]);
        xs.push_back(cc.cell[k]);
        std::sort(xs.begin(), xs.end());
        min.cell[k] = xs[0];
        max.cell[k] = xs[2];
        
        
        // if this triangle straddles the split plane then set the
        // split plane to be the max or min
        // This may then change the AABB on each side of teh split plane
        // so the min and max in the other two dims need to be reconsidered
        //
        if(min.cell[dim] < pos && max.cell[dim] > pos){
            
            
            min.x = (aa.x < min.x) ? aa.x : min.x ;
            min.x = (bb.x < min.x) ? bb.x : min.x ;
            min.x = (cc.x < min.x) ? cc.x : min.x ;
            min.y = (aa.y < min.y) ? aa.y : min.y ;
            min.y = (bb.y < min.y) ? bb.y : min.y ;
            min.y = (cc.y < min.y) ? cc.y : min.y ;
            min.z = (aa.z < min.z) ? aa.z : min.z ;
            min.z = (bb.z < min.z) ? bb.z : min.z ;
            min.z = (cc.z < min.z) ? cc.z : min.z ;
            
            max.x = (aa.x > max.x) ? aa.x : max.x ;
            max.x = (bb.x > max.x) ? bb.x : max.x ;
            max.x = (cc.x > max.x) ? cc.x : max.x ;
            max.y = (aa.y > max.y) ? aa.y : max.y ;
            max.y = (bb.y > max.y) ? bb.y : max.y ;
            max.y = (cc.y > max.y) ? cc.y : max.y ;
            max.z = (aa.z > max.z) ? aa.z : max.z ;
            max.z = (bb.z > max.z) ? bb.z : max.z ;
            max.z = (cc.z > max.z) ? cc.z : max.z ;
        
        }
        
        // If the triangle does not straddle teh split plane then
        // we can just use teh calculated AABB
        //
        
        
        
        
    }
    
}*/

// matrices is a function that calculates the coordinate tranformation matrices
// for converting from a global coordinate system to a local system where the
// 'a' vertex is the origin and the side AB is the ordinate axis.
// This is useful for calculating the surface integral of a triangle
// localToGlobal and globalToLocal are both 9 element (3x3) arrays. The
// memory for them must be allocated before calling this function
//
void Triangle::matrices(double *localToGlobal, double *globalToLocal){
    
    double alpha, beta, T_dash[9], T_dashdash[9] ;
    SPVector zhat = {0,0,1.0} ;
    
    if (localToGlobal == NULL || globalToLocal == NULL) {
        printf("Error matrices must be allocated before calling Triangle::matrices()\n");
        exit(1);
    }
    
    alpha         = atan2(N.y, N.x);
    beta          = acos(VECT_DOT(zhat, N)) ;
    T_dash[0]     = cos(alpha);
    T_dash[1]     = sin(alpha);
    T_dash[2]     = 0;
    T_dash[3]     = -sin(alpha);
    T_dash[4]     = cos(alpha);
    T_dash[5]     = 0;
    T_dash[6]     = 0;
    T_dash[7]     = 0;
    T_dash[8]     = 1;
    T_dashdash[0] = cos(beta);
    T_dashdash[1] = 0;
    T_dashdash[2] = -sin(beta);
    T_dashdash[3] = 0;
    T_dashdash[4] = 1;
    T_dashdash[5] = 0;
    T_dashdash[6] = sin(beta);
    T_dashdash[7] = 0;
    T_dashdash[8] = cos(beta);
    
    matmul(T_dashdash, T_dash, globalToLocal, 3, 3, 3, 3);
    mat3by3inv(globalToLocal, localToGlobal) ;
    
    return ;
}



