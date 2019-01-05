/***************************************************************************
 * 
 *           Module :  TriangleMesh.hpp
 *          Program :  Sarcastic
 *       Created by :  Darren Muff on 18/09/2016
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  03-Nov-2018
 *   Description:
 *      Class header for TriangleMesh class
 * 
 *   (c) Crown Copyright 2018 Defence Science and Technology Laboratory
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software")
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
 ***************************************************************************/

#ifndef TriangleMesh_hpp
#define TriangleMesh_hpp

#include <stdio.h>
#include "materialProperties.h"
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <sarclib/sarclib.h>
#include <set>
#include "AABB.hpp"
#include <algorithm>

enum  TRIVERTEX {
    vertA = 0,
    vertB = 1,
    vertC = 2
} ;

void nad(SPVector aa, SPVector bb, SPVector cc, SPVector *normal, float *area, float *distance) ;

class halfEdge {
public:
    int vertex = -1 ;    // index of vertex owning (starting) this halfEdge
    int nextVertex = -1 ;
    int face = -1 ;      // index of face (triangle) associated with this halfEdge
    int nextHalfEdge     = -1 ;
    int oppositeHalfEdge = -1 ;
    
    halfEdge(){} ;
    
    bool operator<(const halfEdge &o) const {
        if(vertex != o.vertex){
            return vertex < o.vertex;
        }
        return nextVertex < o.nextVertex ;
    }
    
    bool operator==(const halfEdge &o) const {
        if(vertex == o.vertex && nextVertex == o.nextVertex){
            return true;
        }
        return false;
    }
    
};


class rawTri {
public:
    SPVector aa, bb, cc, N;
    int mat=0 ; // material index
    float A=0 ; // area
    float d=0 ; // distance
    
    rawTri(){}
    
    rawTri(SPVector vertA, SPVector vertB, SPVector vertC) : aa(vertA), bb(vertB), cc(vertC)  {
        nad(aa, bb, cc, &N, &A, &d) ;
    } ;
    rawTri(SPVector vertA, SPVector vertB, SPVector vertC, std::string material) : aa(vertA), bb(vertB), cc(vertC) {
        nad(aa, bb, cc, &N, &A, &d) ;
        setMaterial(material) ;
        return ;
    };
    rawTri(SPVector vertA, SPVector vertB, SPVector vertC, int material) : aa(vertA), bb(vertB), cc(vertC), mat(material) {
        nad(aa, bb, cc, &N, &A, &d) ;
        setMaterial(material) ;
        return ;
    };
    
    void setMaterial(int matId) { mat = matId; };
    
    void setMaterial(std::string material){
        std::transform(material.begin(), material.end(),material.begin(), ::toupper);
        for(int i=0; i<globalNumMats; i++){
            if (material.find(std::string(globalMatProps[i].matname)) != std::string::npos ) {
                mat = i ;
                i = globalNumMats ;
            }
        }
        return ;
    }
    
    void print(){
        std::cout << aa.x << " " << aa.y << " " << aa.z << std::endl;
        std::cout << bb.x << " " << bb.y << " " << bb.z << std::endl;
        std::cout << cc.x << " " << cc.y << " " << cc.z << std::endl;
        std::cout << aa.x << " " << aa.y << " " << aa.z << std::endl;
        return ;
    }
};

class Triangle3DVec {
public:
    double x;
    double y;
    double z;
    
    Triangle3DVec() {}
    Triangle3DVec(double x, double y, double z) : x(x), y(y), z(z) {}
    Triangle3DVec(SPVector N) : x(N.x), y(N.y), z(N.z) {}
    SPVector asSPVector(){SPVector v; VECT_CREATE(x, y, z, v); return v;}
    
    bool operator<(const Triangle3DVec &o) const {
        if (x != o.x) {
            return x < o.x;
        }
        if (y != o.y) {
            return y < o.y;
        }
        return z < o.z;
    }
    
    bool operator==(const Triangle3DVec &o) const {
        double epsi = 0.0000001;
        if(fabs(x-o.x)<epsi && fabs(y-o.y)<epsi && fabs(z-o.z)<epsi){
            return true;
        }
        return false;
    }
};

class Triangle {
    
public:
    int a,b,c;          // Index of Vertices in a seperate list of xyz positions (eg TriangleVertex's)
    int mat=0;          // Index of material type
    Triangle3DVec N;    // Triangle Normal
    float dist=0;       // Distance of triangle plane from origin (a.N) Useful for sorting
    float Area=0;       // Area
    double glob2locMatrix[9];
    double loc2GlobMatrix[9];
    
    Triangle(){}
    Triangle(int vertexA, int vertexB, int vertexC): a(vertexA), b(vertexB), c(vertexC) {
    }
    Triangle(int vertexA, int vertexB, int vertexC, int matID): a(vertexA), b(vertexB), c(vertexC), mat(matID) {}
    Triangle(int vertexA, int vertexB, int vertexC, Triangle3DVec N, float dist): a(vertexA), b(vertexB), c(vertexC), N(N), dist(dist) {
        buildTransMats() ;
    }
    Triangle(int vertexA, int vertexB, int vertexC, SPVector N, float dist): a(vertexA), b(vertexB), c(vertexC), N(N), dist(dist) {
        buildTransMats() ;
    }
    Triangle(int vertexA, int vertexB, int vertexC, int matID, Triangle3DVec N, float dist): a(vertexA), b(vertexB), c(vertexC), mat(matID), N(N), dist(dist) {
        buildTransMats() ;
    }
    Triangle(int vertexA, int vertexB, int vertexC, int matID, SPVector N, float dist): a(vertexA), b(vertexB), c(vertexC), mat(matID), N(N), dist(dist) {
        buildTransMats() ;
    }
    Triangle( const Triangle &tri)  // copy constructor
    {
        a=tri.a;
        b=tri.b;
        c=tri.c;
        mat = tri.mat;
        N.x = tri.N.x ;
        N.y = tri.N.y ;
        N.z = tri.N.z ;
        dist = tri.dist ;
        Area = tri.Area ;
        for (int i=0; i<9; ++i){
            glob2locMatrix[i] = tri.glob2locMatrix[i] ;
            loc2GlobMatrix[i] = tri.loc2GlobMatrix[i] ;
        }
    }
    
    // matrices is a function that calculates the coordinate tranformation matrices
    // for converting from a global coordinate system to a local system where the
    // 'a' vertex is the origin and the side AB is the ordinate axis.
    // This is useful for calculating the surface integral of a triangle
    // localToGlobal and globalToLocal are both 9 element (3x3) arrays. The
    // memory for them must be allocated before calling this function
    //
    void matrices(double *localToGlobal, double *globalToLocal) ;
    void buildTransMats(){
        matrices(loc2GlobMatrix, glob2locMatrix) ;
        return ;
    }
    
    // Sort in following order : Material, normal, d, a,b,c
    //
    bool operator<(const Triangle &o) const {
        if( mat != o.mat)   return mat < o.mat ;
        if( ! (N == o.N))   return N < o.N ;
        if( dist != o.dist )      return dist < o.dist ;
        return Area < o.Area;
    }
    
    // Triangles are equal if the have the same material, normal, d and verts in any order
    //
    bool operator==(const Triangle &o) const {
        if (mat != o.mat)   return false;
        if ( dist !=o.dist )      return false ;
        if (!(N == o.N))    return false ;
        std::set<int> s1, s2;
        s1.insert(a); s1.insert(b); s1.insert(c);
        s2.insert(o.a); s2.insert(o.b); s2.insert(o.c);
        if(s1==s2){
            return (s1==s2);
        }else{
            return (s1==s2);
        }
    }
    
    bool coplanar(const Triangle &o) const {
        if ( mat != o.mat ) return false ;
        if (!(N == o.N))    return false ;
        if ( fabs(dist-o.dist) > 1e-8 ) return false ;
        return true;
    }
};

class TriangleMesh {
private:
    std::vector<rawTri> rawTriBuff ;
    bool halfEdgesBuilt = false ;
    bool isPointOnSegment(int pt, int vert1, int vert2) ;
    bool monogamous = false ;

public:
    bool sorted = false ;
    std::vector<Triangle3DVec> vertices ;
    std::vector<Triangle> triangles;
    std::vector<halfEdge> halfedges;
    std::vector<AABB> AABBs ;               // array of AABBs for each triangle
    std::vector<Triangle3DVec> centres;     // centre of each triangle.
    
    TriangleMesh(){}
    
    TriangleMesh(std::vector<rawTri> rawTriangles) : rawTriBuff(rawTriangles){
        sortTrianglesAndPoints() ;
        return;
    }
    
    TriangleMesh(std::vector<Triangle> tris, std::vector<Triangle3DVec> verts) {
        SPVector vertA, vertB, vertC;
        auto it = tris.begin() ;
        while (it != tris.end()){
            VECT_CREATE(verts[it->a].x, verts[it->a].y, verts[it->a].z, vertA);
            VECT_CREATE(verts[it->b].x, verts[it->b].y, verts[it->b].z, vertB);
            VECT_CREATE(verts[it->c].x, verts[it->c].y, verts[it->c].z, vertC);
            rawTri t = rawTri(vertA, vertB, vertC, it->mat) ;
            rawTriBuff.push_back(t) ;
            ++it;
        }
        sortTrianglesAndPoints() ;
        return ;
    }
    
    ~TriangleMesh(){
        rawTriBuff.clear() ;
        vertices.clear() ;
        triangles.clear() ;
        halfedges.clear() ;
    }
    
    void readDAEFile  ( std::string filename );
    void readPLYFile  ( std::string filename );
    void writePLYFile ( std::string filename );
    void readTriFile  ( std::string filename );
    void writeTriFile ( std::string filename );
    void addTriangle  ( int a, int b, int c  );
    void addTriangle  ( int a, int b, int c, int mat );
    void addTriangle(SPVector AA, SPVector BB, SPVector CC, int mat) ;
    void addTriangle  (Triangle tri );
    void addTriangle  (rawTri tri );
    void buildMeshFromRaw (std::vector<rawTri> rawTriangles);
    void buildRawTris();
    void sortTrianglesAndPoints() ;
    void buildHalfEdges() ;
    void printTriangles() ;
    void checkIntegrityAndRepair() ;
    rawTri asRawTriangle(long int triangleIndex);
    std::vector<halfEdge> edges();
    void buildTriangleAABBs();
    TriangleMesh add(const TriangleMesh &newMesh) ;
    SPVector vertAforTri(int triIndx) { return vertices[triangles[triIndx].a].asSPVector() ; }
    SPVector vertBforTri(int triIndx) { return vertices[triangles[triIndx].b].asSPVector() ; }
    SPVector vertCforTri(int triIndx) { return vertices[triangles[triIndx].c].asSPVector() ; }
    long int size() ;
    void buildTrianglelCentres() ;
    void monogamise();
    void removeDegenerates() ;
    bool isMonogamous(){return monogamous;};


    
//    void buildTriangleAABBs(int dim, float pos);
    
};


#endif /* TriangleMesh_hpp */
