/***************************************************************************
 *
 *       Module:    trianglecpp.h
 *      Program:    SARCASTIC
 *   Created by:    Darren Muff on 21/07/2015.
 *                  Copyright (c) 2015 Dstl. All rights reserved.
 *
 *   Description:
 *      Class header for trianglecpp class
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

#ifndef __sarcastic__trianglecpp__
#define __sarcastic__trianglecpp__

#include <iostream>
#include <stdio.h>
#include <sarclib/sarclib.h>
#include <set>

class Triangle {
   
public:
    Triangle();
    Triangle(SPVector aa, SPVector bb, SPVector cc) ;
    Triangle(SPVector aa, SPVector bb, SPVector cc, int material) ;
    Triangle(SPVector aa, SPVector bb, SPVector cc, std::string material) ;
    ~Triangle();
    
    void print();
    
    void setMaterial(int matId);
    void setMaterial(std::string material);
    
    SPVector   AA ;
    SPVector   BB ;
    SPVector   CC ;
    SPVector   NN ;
    double     area ;
    int        matId;
    double     d; // distance of triangle plane to origin
    double     globalToLocalMat[9];
    double     localToGlobalMat[9];
};

struct Point {
public:

    float x;
    float y;
    float z;
    Point() {}
    Point(float x, float y, float z) : x(x), y(y), z(z) {}
    
    bool operator<(const Point &o) const {
        if (x != o.x) {
            return x < o.x;
        }
        if (y != o.y) {
            return y < o.y;
        }
        return z < o.z;
    }
    
    bool operator==(const Point &o) const {
        if(x == o.x && y==o.y && z==o.z){
            return true;
        }
        return false;
    }
    
};

class triangleReference{
public:

    int AA;
    int BB;
    int CC;
    int mat ;
    Point NN ;
    
    triangleReference(){}
    triangleReference(int AA, int BB, int CC, int mat, Point NN) : AA(AA), BB(BB), CC(CC), mat(mat), NN(NN) {}
    
    bool operator<(const triangleReference &o) const {
        if( mat != o.mat) {
            return mat < o.mat ;
        }
        if( ! (NN == o.NN) ) {
            return NN < o.NN ;
        }
        if( AA != o.AA ){
            return AA < o.AA ;
        }
        if( BB != o.BB ){
            return BB < o.BB ;
        }
        return CC < o.CC ;
    }
    
    bool operator==(const triangleReference &o) const {
        if (mat != o.mat){
            return false;
        }
        if (!(NN == o.NN)){
            return false ;
        }
        std::set<int> s1, s2;
        s1.insert(AA); s1.insert(BB); s1.insert(CC);
        s2.insert(o.AA); s2.insert(o.BB); s2.insert(o.CC);
        return (s1==s2);
    }
    
} ;

#endif /* defined(__sarcastic__trianglecpp__) */
