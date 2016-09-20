//
//  TriangleMesh.hpp
//  sarcastic
//
//  Created by Darren Muff on 18/09/2016.
//  Copyright © 2016 Dstl. All rights reserved.
//

#ifndef TriangleMesh_hpp
#define TriangleMesh_hpp

#include <stdio.h>
#include "materialProperties.h"
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <SIlib2/SIlib2.h>
#include <set>

class rawTri {
public:
    SPVector aa, bb, cc, N;
    int mat ; // matreial index
    float A ; // area
    float d ; // distance
    
    rawTri(){}
    
    rawTri(SPVector vertA, SPVector vertB, SPVector vertC) : aa(vertA), bb(vertB), cc(vertC)  {
        SPVector ab; VECT_SUB(bb, aa, ab);
        SPVector bc; VECT_SUB(cc, bb, bc);
        SPVector x;  VECT_CROSS(ab, bc, x);
        double l ; l = VECT_MAG(x);
        A = l * 0.5;
        VECT_SCMULT(x, (1/l), N);
        d = VECT_DOT(aa, N) ;
    } ;
    rawTri(SPVector vertA, SPVector vertB, SPVector vertC, std::string material) : aa(vertA), bb(vertB), cc(vertC) {
        *this = rawTri(vertA, vertB, vertC) ;
        setMaterial(material) ;
        return ;
    };
    rawTri(SPVector vertA, SPVector vertB, SPVector vertC, int material) : aa(vertA), bb(vertB), cc(vertC), mat(material) {
        *this = rawTri(vertA, vertB, vertC) ;
        setMaterial(material) ;
        return ;
    };

    void setMaterial(int matId) { mat = matId; };
    
    void setMaterial(std::string material){
        std::transform(material.begin(), material.end(),material.begin(), ::toupper);
        for(int i=0; i<NMATERIALS; i++){
            if (material.find(std::string(materialProperties[i].matname)) != std::string::npos ) {
                mat = i ;
                i = NMATERIALS ;
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
        if(x == o.x && y==o.y && z==o.z){
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
    double d;           // Distance of triangle plane from origin (a.N) Useful for sorting
    double A;           // Area
    
    Triangle(){}
    Triangle(int vertexA, int vertexB, int vertexC): a(vertexA), b(vertexB), c(vertexC) {}
    Triangle(int vertexA, int vertexB, int vertexC, int matID): a(vertexA), b(vertexB), c(vertexC), mat(matID) {}
    Triangle(int vertexA, int vertexB, int vertexC, Triangle3DVec N): a(vertexA), b(vertexB), c(vertexC), N(N) {}
    Triangle(int vertexA, int vertexB, int vertexC, SPVector N): a(vertexA), b(vertexB), c(vertexC), N(N) {}
    Triangle(int vertexA, int vertexB, int vertexC, int matID, Triangle3DVec N): a(vertexA), b(vertexB), c(vertexC), mat(matID), N(N) {}
    Triangle(int vertexA, int vertexB, int vertexC, int matID, SPVector N): a(vertexA), b(vertexB), c(vertexC), mat(matID), N(N) {}

    
    // Sort in following order : Material, normal, d, a,b,c
    //
    bool operator<(const Triangle &o) const {
        if( mat != o.mat)   return mat < o.mat ;
        if( ! (N == o.N))   return N < o.N ;
        if( d != o.d )      return d < o.d ;
        if( a != o.a )      return a < o.a ;
        if( b != o.b )      return b < o.b ;
        return c < o.c ;
    }
    
    // Triangles are equal if the have the same material, normal, d and verts in any order
    //
    bool operator==(const Triangle &o) const {
        if (mat != o.mat)   return false;
        if ( d !=o.d )      return false ;
        if (!(N == o.N))    return false ;
        std::set<int> s1, s2;
        s1.insert(a); s1.insert(b); s1.insert(c);
        s2.insert(o.a); s2.insert(o.b); s2.insert(o.c);
        return (s1==s2);
    }
};

class TriangleMesh {
private:
    std::vector<rawTri> rawTriBuff ;
    void sortTrianglesAndPoints() ;
    bool sorted ;
    
public:
    std::vector<Triangle3DVec> vertices ;
    std::vector<Triangle> triangles;
    
    TriangleMesh(){}
    
    void readDAEFile  ( std::string filename );
    void readPLYFile  ( std::string filename );
    void writePLYFile ( std::string filename );
    void readTriFile  ( std::string filename );
    void writeTriFile ( std::string filename );
    
};


#endif /* TriangleMesh_hpp */
