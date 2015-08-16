//
//  main.cpp
//  materialise
//
//  Created by Darren on 12/08/2015.
//  Copyright (c) 2015 Dstl. All rights reserved.
//
#define REAL double
#define ANSI_DECLARATORS

#include <iostream>
#include <SIlib/SIlib.h>
#include "trianglecpp.h"
#include "TriangleFile.h"
extern "C" {
#include "matrixMultiplication.h"
#include "JRStriangle.h"
#include "boxMullerRandom.h"
}
#include "materialProperties.h"
#include <sstream>

#define ROOTPATH "/Users/Darren/Development"

bool VECT_EQUAL(SPVector a, SPVector b);
void roughenPoint(SPVector *p, Triangle t);
bool pointsColinear(SPVector p1, SPVector p2, SPVector p3);

int main(int argc, const char * argv[]) {

    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, "materialise", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status) ;
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
                               (char *) ROOTPATH"/DATA/triangles_Delaunay.tri");
    
    FILE *fpout ;
    fpout = fopen(oustr, "w");
    if (fpout == NULL) {
        printf("Error : failed to open input triangle file %s\n",oustr);
        exit(1);
    }
    
    bool plyop=false;
    char *plyFname;
    plyop = input_yesno("Write out triangles as PLY file", "plyop",
                        "The PLY format is a simple-to-understand representation of the triangles",
                        plyop);
    if(plyop){
        plyFname = input_string("name of PLY file to create", "plyFname",
                                "Full path name of the PLY file to create", plyFname);
    }

    
    
    srand((unsigned int)time(NULL));

    TriangleFile triFile = TriangleFile(std::string(instr));
    
    std::vector<Triangle> newTriangleVec;
    SPVector *Points3D ;
    
    for(int t=0; t<triFile.triangles.size(); t++){
        Triangle tri = triFile.triangles[t] ;
        
        // convert 3D coords to 2D plane.
        // AA is the origin
        // Ordinate axis is defined as AABB
        //
        SPVector N = tri.NN;
        SPVector ordinate,abscissa;
        VECT_SUB(tri.BB, tri.AA, ordinate);
        VECT_NORM(ordinate, ordinate);
        VECT_CROSS(N, ordinate, abscissa);
        VECT_NORM(abscissa, abscissa);
        
        SPVector v1,v2;
        double x0,y0,x1,y1,x2,y2;
        
        VECT_SUB(tri.BB, tri.AA, v1);
        VECT_SUB(tri.CC, tri.AA, v2);

        x0 = y0 = 0.0;
        x1 = VECT_DOT(v1, ordinate);
        y1 = VECT_DOT(v1, abscissa);
        x2 = VECT_DOT(v2, ordinate);
        y2 = VECT_DOT(v2, abscissa);
        
        // put 2D point in pointlist
        //
        double points[6];
        points[0] = x0; points[1] = y0;
        points[2] = x1; points[3] = y1;
        points[4] = x2; points[5] = y2;
        
        // find material and from it max area
        //
        double corrlen = materialProperties[tri.matId].corlen ;
        double corrArea = corrlen ;
        std::ostringstream stringstream ;
        stringstream << "zQa" << corrArea ;
        std::string str = stringstream.str();
        const char *triswitches = str.c_str() ;
    
        // delauney triangulate
        //
        struct triangulateio in ;
        struct triangulateio out;
        in.numberofpointattributes = 0;
        in.pointmarkerlist = NULL ;
        in.numberofcorners = 3;
        in.numberoftriangleattributes = 0;
        in.trianglearealist = NULL;
        in.pointlist = points ;
        in.numberofpoints = 3 ;
        out.pointlist = NULL;
        out.pointmarkerlist = NULL ;
        out.trianglelist = NULL ;
        
        triangulate((char *)triswitches, &in, &out, NULL);
        
        // convert all points to 3D points. Do it here so that we can
        // tweek the internal points for surface rouchness.
        
        SPVector OO;
        OO = tri.AA;
        Points3D = (SPVector *)sp_malloc(sizeof(SPVector) * out.numberofpoints);
        
        for (int n=0; n<out.numberofpoints; n++){
            double x2d,y2d ;
            SPVector v,v1,v2 ;
            x2d = out.pointlist[2*n+0] ;
            y2d = out.pointlist[2*n+1] ;
            VECT_SCMULT(ordinate, x2d, v1);
            VECT_SCMULT(abscissa, y2d, v2);
            VECT_ADD(v1, v2, v);
            VECT_ADD(v, OO, v);
            roughenPoint(&v, tri) ;
            Points3D[n] = v ;
        }
        
        // and put into newtri vector
        //
        int ntri = out.numberoftriangles ;
        for (int n=0; n<ntri; n++){
            
            int pointind0 = out.trianglelist[n*3+0] ;
            int pointind1 = out.trianglelist[n*3+1] ;
            int pointind2 = out.trianglelist[n*3+2] ;
            
            SPVector AA,BB,CC ;
            AA = Points3D[pointind0] ;
            BB = Points3D[pointind1] ;
            CC = Points3D[pointind2] ;
            
            Triangle newtri = Triangle(AA, BB, CC, tri.matId);
            newTriangleVec.push_back(newtri) ;
        }
        free(Points3D);
    }
    
    // newTriangleVec now contains new triangles.
    // write it to file
    //
    TriangleFile dtTriFile = TriangleFile(newTriangleVec) ;
    dtTriFile.WriteFile(std::string(oustr));
    if(plyop){
        dtTriFile.WritePLYFile(plyFname, false);
    }
    im_close_lib(&status);
    std::cout << "Done" << std::endl;

    return 0;
}

void roughenPoint(SPVector *p, Triangle t){
    SPVector V ;
    SPVector vin;
    
    vin.x = p->x ;
    vin.y = p->y ;
    vin.z = p->z ;

    if ((!VECT_EQUAL(vin, t.AA)) && (!VECT_EQUAL(vin, t.BB)) && (!VECT_EQUAL(vin, t.CC)) ) {
        
        // Check that vin is not on one of the borders of the original triangle
        //
        if ( !(pointsColinear(vin, t.AA, t.BB) || pointsColinear(vin, t.BB, t.CC) || pointsColinear(vin, t.CC, t.AA))) {
            
            float dist = box_muller(0.0, materialProperties[t.matId].roughness) ;
            VECT_SCMULT(t.NN, dist, V);
            VECT_ADD(V, vin, vin);
            p->x = vin.x ;
            p->y = vin.y ;
            p->z = vin.z ;
        }
    }
    return ;
}

bool pointsColinear(SPVector point, SPVector A, SPVector B){
    double EPSILON = 1.0e-6;
    
    SPVector AB, Ap, Bp, v;
    
    VECT_SUB(B, A, AB) ;
    VECT_SUB(point, A, Ap) ;
    VECT_SUB(point, B, Bp) ;
    
    if (VECT_MAG(Ap) > EPSILON) {
        v=Ap ;
    }else if(VECT_MAG(Bp) > EPSILON){
        v=Bp ;
    }else{
        return true;
    }
    
    VECT_NORM(v, v) ;
    VECT_NORM(AB, AB) ;
    
    double dot = fabs(VECT_DOT(v, AB)) ;
    
    if( fabs(dot-1.0) < EPSILON ){
        return true;
    }else{
        return false;
    }
}


bool VECT_EQUAL(SPVector a, SPVector b){
    return (a.x == b.x && a.y == b.y && a.z == b.z) ;
}
