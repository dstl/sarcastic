//
//  main.cpp
//  materialise
//
//  Created by Darren on 12/08/2015.
//  Copyright (c) 2015 Dstl. All rights reserved.
//
#define REAL double
#define ANSI_DECLARATORS

#include "materialise_version.h"
#include <iostream>
#include <SIlib2/SIlib2.h>
//#include "trianglecpp.h"
//#include "TriangleFile.h"
#include "TriangleMesh.hpp"
extern "C" {
#include "matrixMultiplication.h"
#include "JRStriangle.h"
#include "boxMullerRandom.h"
}
#include "materialProperties.h"
#include "colourCodes.h"
#include <sstream>
#include <forward_list>

#define ROOTPATH "/tmp"
#define REPORTN 100
#define VERBOSE ((int) 0)

using namespace std;

void roughenPoint(SPVector *p, Triangle t);
TriangleMesh growTriangles(TriangleMesh *mesh) ;
void meshTo2D(TriangleMesh *mesh, SPVector *org, SPVector *ord, SPVector *absci,
              int **trinds,int **segments,int *nsegments, double **points, int *nholes, double **holes) ;
void cleanPoints(int ntriangles, double **points, int **trinds, int **segments, int *nPoints, int *nSegs, double **holes, int *nholes);
void PointTo2D(SPVector vec, SPVector origin, SPVector ordinate, SPVector abscissa, double *x, double *y);
bool isPointOnLine(double px, double py, double x1, double y1, double x2, double y2, double approxZero);
char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def);
void printProgress(int percentDone, int nChars) ;
bool pointInPolygon(int npoints, double *points, int nsegments, int *segments, double x, double y);
bool pointInPolygon(int npoints, double *points, int nsegments, int *segments, double x, double y, int *nleft, int *nright);
void testcommonmesh(TriangleMesh *mesh) ;


int main(int argc, const char * argv[]) {
    
    int verbose = VERBOSE ;
    int nSegments, nholes ;
    SPVector origin, abscissa, ordinate, zdir, *Points3D;
    int *trinds;
    int *segments ;
    double *holes ;
    double *points ;

    printf(" \n");
    printf(DARK GREEN "                              Materialise\n" NORMAL);
    printf(BLUE "                          Version :" RED" %s \n", MATERIALISE_FULL_VERSION);
    printf(BLUE "                    Revision: " RED"%s, %s \n",MATERIALISE_REVISION, MATERIALISE_VERSION_DATE);
    printf(BLUE "               Copyright (c) 2013 " WHITE"[" BLUE"Dstl" WHITE"]" BLUE". All rights reserved.\n" RESETCOLOR);
    printf(" \n");

    
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, "materialise", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status) ;
    char *instr = tryReadFile((char *)"Input triangle filename", (char *)"infilename",
                               (char *)"The name of a triangle file created with 'colladaToTriFile'",
                               (char *) ROOTPATH"/triangles.ply");
    
    char *oustr = input_string((char *)"Output triangle filename", (char *)"oufilename",
                               (char *)"The name of the triangle file to be created",
                               (char *) ROOTPATH"/delaunay.ply");
    
    FILE *fpin ;
    fpin = fopen(instr, "r");
    if (fpin == NULL) {
        printf("Error : failed to open input triangle file %s\n",instr);
        exit(1);
    }
    
    FILE *fpout ;
    fpout = fopen(oustr, "w");
    if (fpout == NULL) {
        printf("Error : failed to open input triangle file %s\n",oustr);
        exit(1);
    }
    
    TriangleMesh mesh, newMesh;
    mesh.readPLYFile(string(instr));
    mesh.checkIntegrityAndRepair();
    
    mesh.writePLYFile(std::string("/tmp/junk.ply"));

    mesh.sortTrianglesAndPoints() ;
    long int orgNTris;
    float percentDone;
    orgNTris = mesh.triangles.size() ;
    Timer timer;
    startTimer(&timer, &status);
    int cnt=0;
    while (!mesh.triangles.empty()) {
        
        percentDone = 100.0 * (float)(orgNTris - mesh.triangles.size()) / (float)orgNTris ;
        printProgress(percentDone, 80);
        
        // Take the first triangle in mesh.triangles and calculate all triangles that
        // are common to it (on the same plane and joined togethr
        //
        TriangleMesh commonMesh = growTriangles(&mesh) ;
        
        commonMesh.writePLYFile(std::string("/tmp/dbg")+to_string(cnt)+std::string(".ply"));
        
        // Convert to 2Dpoints
        //
        meshTo2D(&commonMesh, &origin, &ordinate, &abscissa, &trinds, &segments, &nSegments, &points, &nholes, &holes);
        VECT_CROSS(ordinate, abscissa, zdir);
        VECT_NORM(zdir, zdir);
        
        int npoints = (int)commonMesh.vertices.size() ;
        int ntrinds = (int)commonMesh.triangles.size() ;
        
        // find material and from it max area
        //
        double corrlen = materialProperties[commonMesh.triangles[0].mat].corlen ;
        double corrArea = corrlen * corrlen * 0.5 ;
        if(verbose)printf("setting max area to be %f m^2 (correlation length is %f m)\n",corrArea,corrlen);
        
        // Set up switches for Delaunay triangulation
        //
        std::ostringstream stringstream ;
        stringstream << "zpj" ;
        if (!verbose) {
            stringstream << "Q" ;
        }
        stringstream << "a" << corrArea ;
        std::string str = stringstream.str();
        const char *triswitches = str.c_str() ;
        
        // Perform Delaunay triangulation
        //
        struct triangulateio in ;
        struct triangulateio out;
        in.pointlist = points ;
        in.numberofpoints = npoints ;
        in.numberofpointattributes = 0;
        in.pointmarkerlist = NULL ;
        in.pointattributelist = NULL;
        
        in.trianglelist = trinds ;
        in.numberoftriangles = ntrinds ;
        in.trianglearealist = NULL;
        in.numberoftriangleattributes = 0;
        in.triangleattributelist = NULL;
        in.numberofcorners = 3;
        
        in.segmentlist = segments ;
        in.numberofsegments = nSegments ;
        in.segmentmarkerlist = NULL ;
        
        in.numberofholes = nholes;
        in.holelist = holes;
        in.numberofregions = 0;
        in.regionlist = NULL;
        
        out.pointlist = NULL;
        out.pointmarkerlist = NULL ;
        out.trianglelist = NULL ;
        out.neighborlist = NULL ;
        out.segmentlist = NULL;
        out.segmentmarkerlist = NULL;
        out.holelist = NULL;
        
        if(verbose) {
            printf("Triangles\n");
            for (int i=0; i<in.numberoftriangles; ++i) {
                printf(" %f,%f\n %f,%f\n %f,%f\n %f,%f\n",
                       in.pointlist[in.trianglelist[i*3+0]*2+0], in.pointlist[in.trianglelist[i*3+0]*2+1],
                       in.pointlist[in.trianglelist[i*3+1]*2+0], in.pointlist[in.trianglelist[i*3+1]*2+1],
                       in.pointlist[in.trianglelist[i*3+2]*2+0], in.pointlist[in.trianglelist[i*3+2]*2+1],
                       in.pointlist[in.trianglelist[i*3+0]*2+0], in.pointlist[in.trianglelist[i*3+0]*2+1]);

            }
            printf("Done\n");
            printf("segments\n");
            for (int i=0; i<in.numberofsegments; ++i) {
                printf(" %f,%f\n %f,%f\n",
                       in.pointlist[in.segmentlist[i*2+0]*2+0], in.pointlist[in.segmentlist[i*2+0]*2+1],
                       in.pointlist[in.segmentlist[i*2+1]*2+0], in.pointlist[in.segmentlist[i*2+1]*2+1]);
                
            }
            printf("Done\n");
            printf("Holes\n");
            for (int i=0; i<in.numberofholes; ++i) {
                printf(" %f,%f\n",
                       in.holelist[i*2+0], in.holelist[i*2+1]);
            }
            printf("Done\n");
        }
        
        triangulate((char *)triswitches, &in, &out, NULL);
        
        if(verbose) {
            printf("Triangles\n");
            for (int i=0; i<out.numberoftriangles; ++i) {
                printf(" %f,%f\n %f,%f\n %f,%f\n %f,%f\n",
                       out.pointlist[out.trianglelist[i*3+0]*2+0], out.pointlist[out.trianglelist[i*3+0]*2+1],
                       out.pointlist[out.trianglelist[i*3+1]*2+0], out.pointlist[out.trianglelist[i*3+1]*2+1],
                       out.pointlist[out.trianglelist[i*3+2]*2+0], out.pointlist[out.trianglelist[i*3+2]*2+1],
                       out.pointlist[out.trianglelist[i*3+0]*2+0], out.pointlist[out.trianglelist[i*3+0]*2+1]);
                
            }
            printf("Done\n");
            printf("segments\n");
            for (int i=0; i<out.numberofsegments; ++i) {
                printf(" %f,%f\n %f,%f\n",
                       out.pointlist[out.segmentlist[i*2+0]*2+0], out.pointlist[out.segmentlist[i*2+0]*2+1],
                       out.pointlist[out.segmentlist[i*2+1]*2+0], out.pointlist[out.segmentlist[i*2+1]*2+1]);
                
            }
            printf("Done\n");
            printf("Holes\n");
            for (int i=0; i<out.numberofholes; ++i) {
                printf(" %f,%f\n",
                       out.holelist[i*2+0], out.holelist[i*2+1]);
            }
            printf("Done\n");
        }

        if(out.numberoftriangles < in.numberoftriangles){
            if(verbose)printf("error fewer tris buit than input\n");
            out.numberofpoints = in.numberofpoints ;
            out.numberofsegments = in.numberofsegments ;
            out.numberoftriangles = in.numberoftriangles ;
            free(out.pointlist);
            free(out.segmentlist);
            free(out.trianglelist);
            out.pointlist = (double *)sp_malloc(sizeof(double) * out.numberofpoints * 2) ;
            memcpy(out.pointlist, in.pointlist, sizeof(double) * out.numberofpoints * 2) ;
            out.segmentlist = (int *)sp_malloc(sizeof(int) * out.numberofsegments * 2) ;
            memcpy(out.segmentlist,in.segmentlist,sizeof(int) * out.numberofsegments * 2) ;
            out.trianglelist = (int *)sp_malloc(sizeof(int) * out.numberoftriangles * 3) ;
            memcpy(out.trianglelist, in.trianglelist, sizeof(int) * out.numberoftriangles * 3) ;
        }
        
        
        
        // Roughen surface corresponding to material type
        //
        SPVector OO;
        OO = origin;
        Points3D = (SPVector *)sp_malloc(sizeof(SPVector) * out.numberofpoints);
        
        for (int n=0; n<out.numberofpoints; n++){
            double x2d,y2d, Z;
            bool pntOnSegment = false ;
            SPVector v,v1,v2,v3 ;
            x2d = out.pointlist[2*n+0] ;
            y2d = out.pointlist[2*n+1] ;
            
            int s=0;
            Z=0;
            while (s<out.numberofsegments && pntOnSegment==false) {
                double s1x,s1y,s2x,s2y ;
                s1x = out.pointlist[ out.segmentlist[s*2+0]*2+0] ;
                s1y = out.pointlist[ out.segmentlist[s*2+0]*2+1] ;
                s2x = out.pointlist[ out.segmentlist[s*2+1]*2+0] ;
                s2y = out.pointlist[ out.segmentlist[s*2+1]*2+1] ;
                
                if( isPointOnLine(x2d, y2d, s1x,s1y,s2x,s2y, 0.001) ){
                    pntOnSegment=true;
                }
                ++s ;
            }
            if (!pntOnSegment) {
                Z = box_muller(0.0, materialProperties[commonMesh.triangles[0].mat].roughness) ;
            }
            
            VECT_SCMULT(ordinate, x2d, v1);
            VECT_SCMULT(abscissa, y2d, v2);
            VECT_SCMULT(zdir, Z, v3);
            
            VECT_ADD(v1, v2, v);
            VECT_ADD(v3, v, v);
            VECT_ADD(v, OO, v);
            Points3D[n] = v ;
        }
        
        // Add to output triangleMesh
        //
        for (int n=0; n<out.numberoftriangles; n++){
            int pointind0 = out.trianglelist[n*3+0] ;
            int pointind1 = out.trianglelist[n*3+1] ;
            int pointind2 = out.trianglelist[n*3+2] ;
            
            SPVector AA,BB,CC ;
            AA = Points3D[pointind0] ;
            BB = Points3D[pointind1] ;
            CC = Points3D[pointind2] ;

            newMesh.addTriangle(AA, BB, CC, commonMesh.triangles[0].mat);
        }
        
        free(Points3D);
        free(points);
        free(trinds);
        free(segments);
        if(nholes!=0)free(holes);
        free(out.pointlist);
        free(out.pointattributelist);
        free(out.pointmarkerlist);
        free(out.trianglelist);
        free(out.triangleattributelist);
        free(out.trianglearealist);
        free(out.neighborlist);
        free(out.segmentlist);
        free(out.segmentmarkerlist);
        
        // Next triangle
        //
        cnt++;
    }
    

    printf("Writing %lu triangles to file %s...\n",newMesh.triangles.size(),oustr);
    newMesh.writePLYFile(std::string(oustr)) ;
    
    endTimer(&timer, &status);
    im_close_lib(&status);
    std::cout << "Done in " << timeElapsedInSeconds(&timer, &status) << " seconds." << std::endl;
    
    return 0;
}

bool coplanar(const Triangle &a, const Triangle &b)  {
    if ( a.mat != b.mat ) return a.mat < b.mat ;
    if (!(a.N == b.N))    return a.N < b.N ;
    return a.dist<b.dist ;
}

TriangleMesh growTriangles(TriangleMesh *mesh)
// Takes the first triangle from the Triangles list and finds all the coplanar triangles
// inside 'triangles'. These are removed from 'triangles' and returned as a vector
//
{
    Triangle objectTri(mesh->triangles.front());
    auto bounds = std::equal_range(mesh->triangles.begin(), mesh->triangles.end(), objectTri, coplanar) ;
    vector<Triangle> commonTriangles(bounds.first,bounds.second);
    mesh->triangles.erase(bounds.first, bounds.second);
    TriangleMesh newMesh(commonTriangles, mesh->vertices);
    return newMesh ;
}

void meshTo2D(TriangleMesh *mesh, SPVector *org, SPVector *ord, SPVector *absci,
              int **trinds,int **segments,int *nsegments, double **points, int *nholes, double **holes)
{
    
    Triangle tri,tri2;
    SPVector N,ordinate,abscissa, origin;
    double x0,y0;
    
    *points   = (double *)sp_malloc(sizeof(double) * mesh->vertices.size()*2);
    *trinds   = (int *)sp_malloc(sizeof(int) * mesh->triangles.size() * 3);

    tri = *(mesh->triangles.begin()) ;
    N = tri.N.asSPVector();
    origin = mesh->vertices[tri.a].asSPVector() ;
    VECT_SUB(mesh->vertices[tri.b].asSPVector(), mesh->vertices[tri.a].asSPVector(), ordinate);
    VECT_NORM(ordinate, ordinate);
    VECT_CROSS(N, ordinate, abscissa);
    VECT_NORM(abscissa, abscissa);

    // Convert vertices to 2D points
    //
    for (int iv=0; iv < mesh->vertices.size(); ++iv){
        PointTo2D(mesh->vertices[iv].asSPVector(), origin, ordinate, abscissa, &x0, &y0);
        (*points)[2*iv + 0] = x0 ;
        (*points)[2*iv + 1] = y0 ;
    }
    
    // Return triangle indices in mesh as trinds
    //
    for(int it=0; it<mesh->triangles.size(); ++it){
        (*trinds)[it*3+0] = mesh->triangles[it].a ;
        (*trinds)[it*3+1] = mesh->triangles[it].b ;
        (*trinds)[it*3+2] = mesh->triangles[it].c ;
    }
    
    // Return edges in mesh as segments
    //
    std::vector<halfEdge> edges = mesh->edges() ;
    *nsegments = (int)edges.size() ;
    *segments = (int *)sp_malloc(sizeof(int) * edges.size() * 2);
    for(int ie=0; ie<edges.size(); ++ie){
        (*segments)[2*ie+0] = edges[ie].vertex ;
        (*segments)[2*ie+1] = edges[ie].nextVertex ;
    }
    
    SPVector Aa,Bb,dummyHole,xp;
    SPVector AB,outdir,zhat;
    int a,b, nleft,nright;
    bool withinPoly ;
    std::vector<double> holeVect ;
    
    for(int ie=0; ie<edges.size(); ++ie){
        // create a 'dummy' hole just over the halfedge away from its triangle
        //
        a = edges[ie].vertex ;
        b = edges[ie].nextVertex ;
        
        VECT_CREATE((*points)[a*2+0], (*points)[a*2+1], 0.0, Aa);
        VECT_CREATE((*points)[b*2+0], (*points)[b*2+1], 0.0, Bb);
        VECT_CREATE(0,0,1,zhat);
        VECT_SUB(Bb, Aa, AB);
        VECT_SCMULT(AB, 0.5, AB);
        VECT_CROSS(AB,zhat,outdir) ;
        VECT_ADD(Aa,AB,xp);
        VECT_NORM(outdir, outdir);
        VECT_SCMULT(outdir, 0.01, dummyHole);
        VECT_ADD(xp, dummyHole, dummyHole);
        
        // Now see if the 'dummyHole' is outside the closed polygon defined
        // by the halfedges of this mesh.

        withinPoly = pointInPolygon((int)mesh->vertices.size(), *points, *nsegments, *segments, dummyHole.x, dummyHole.y, &nleft, &nright) ;

        if(!withinPoly){
            if(nleft != 0 && nright != 0){
                holeVect.push_back(dummyHole.x);
                holeVect.push_back(dummyHole.y);
            }
        }
    }
    
    *nholes = (int) holeVect.size() / 2 ;
    if(*nholes > 0){
        *holes = (double *)sp_malloc(sizeof(double) * *nholes * 2) ;
        
        for (int i=0; i < *nholes ; ++i) {
            (*holes)[2*i+0] = holeVect[2*i+0] ;
            (*holes)[2*i+1] = holeVect[2*i+1] ;
        }
    }else{
        *holes = NULL ;
    }
    
    *org   = origin ;
    *ord   = ordinate ;
    *absci = abscissa ;
    
    return ;
}

void PointTo2D(SPVector vec, SPVector origin, SPVector ordinate, SPVector abscissa, double *x, double *y){

    SPVector v ;
    
    VECT_SUB(vec, origin, v);
    *x = VECT_DOT(v, ordinate) ;
    *y = VECT_DOT(v, abscissa) ;

    return;
}

bool isPointOnLine(double px, double py, double x1, double y1, double x2, double y2, double approxZero)
// Line is defined as the segment between points x1,y1 and x2,y2
// point px,py is the point being tested
// approx zero is a small number. if the distance of the point p from the line is less than this then
// it is assumed the point is on the line
//
{
    double x,y,u,d ;
    
    if( x1==x2 && y1==y2){
        if(px == x1 && py == y1) {
            return true;
        }else{
            return false;
        }
    }
    if( (px==x1 && py == y1) || (px==x2 && py==y2) ){
        return true;
    }
    
    u = ((px-x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)) / ( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) ) ;
    
    x = x1 + u * (x2 - x1) ;
    y = y1 + u * (y2 - y1) ;
    
    d = sqrt ( ((x-px)*(x-px)) + ((y-py)*(y-py)) ) ;
    
    if (d<approxZero){
        return true ;
    }else{
        return false ;
    }
}

char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def)
///  prompt  : the text to display when asking the user for input (the default value will also be displayed)
///  key     : a unique bit of text (such az "AzPixelSize") which the input system will then use to store parameter values
///  help    : the text to display to the user when they just enter '?'
///  def     : the default, ie the value to take if the user just presses return
///
{
    FILE *fp;
    SPStatus fileStat ;
    char *prmpt, *fname;
    
    size_t len = strlen(def);
    prmpt = (char *)sp_malloc(sizeof(char) * len);
    
    do {
        im_init_status(fileStat, 0) ;
        strcpy(prmpt, def);
        fname = input_string(prompt, key, help, def);
        
        if ( (fp = fopen(fname, "r")) == NULL){
            printf(RED "Cannot access file %s\n" RESETCOLOR,fname);
            fileStat.status = BAD_FILE ;
        }
        
    } while(fileStat.status != NO_ERROR);
    
    fclose(fp);
    free(prmpt) ;
    return(fname) ;
}

void printProgress(int percentDone, int nChars)
{
    int ppc = (float)100 / nChars ; // Percent per character
    putchar('\r');
    for (int i=0; i<nChars; ++i){
        if (i * ppc < percentDone) {
            putchar('>');
        }else{
            putchar('.');
        }
//        putchar('\n');
    }
    std::cout.flush();
}

bool pointInPolygon(int npoints, double *points, int nsegments, int *segments, double x, double y)
{
    bool  oddNodes=false   ;
    
    for(int i=0; i<nsegments; i++){
        int is, ie;
        is = segments[2*i+0] ;
        ie = segments[2*i+1] ;
        double psx, psy, pex, pey;
        psx = points[2*is+0] ;
        psy = points[2*is+1] ;
        pex = points[2*ie+0] ;
        pey = points[2*ie+1] ;
    
        if( ((psy < y && pey >= y) || (pey < y && psy >= y )) && (psx <=x || pex <=x ) ){
            oddNodes ^= ((psx + ((y-psy) / (pey - psy )) * (pex - psx )) < x) ;
        }
        
    }
    return oddNodes;
}

bool pointInPolygon(int npoints, double *points, int nsegments, int *segments, double x, double y, int *nleft, int *nright)
{
    bool  oddNodes=false   ;
    int l=0,r=0;
    
    for(int i=0; i<nsegments; i++){
        int is, ie;
        is = segments[2*i+0] ;
        ie = segments[2*i+1] ;
        double psx, psy, pex, pey;
        psx = points[2*is+0] ;
        psy = points[2*is+1] ;
        pex = points[2*ie+0] ;
        pey = points[2*ie+1] ;
        
        if( ((psy < y && pey >= y) || (pey < y && psy >= y )) && (psx <=x || pex <=x ) ){
            if((psx + ((y-psy) / (pey - psy )) * (pex - psx )) < x){
                ++l;
                oddNodes = !oddNodes ;
            }else{
                ++r;
            }
        }
    }
    *nleft = l;
    *nright = r;
    return oddNodes;
}


void testcommonmesh(TriangleMesh *mesh){
    SPVector N = mesh->triangles[0].N.asSPVector() ;
    double dist = mesh->triangles[0].dist ;
    int mat = mesh->triangles[0].mat ;
    bool pass = true;
    
    for (int i=0 ; i< mesh->triangles.size(); ++i){
        if(!(mesh->triangles[i].N == N) || mesh->triangles[i].dist != dist || mesh->triangles[i].mat != mat){
            printf("Error: Triangles not the same in mesh\n");
            printf("[%d] N: %f,%f,%f (should be (%f,%f,%f), dist: %f (should be %f), mat: %d (should be %d)\n",i
                   ,mesh->triangles[i].N.x,mesh->triangles[i].N.y,mesh->triangles[i].N.z,N.x,N.y,N.z,mesh->triangles[i].dist,dist,
                   mesh->triangles[i].mat, mat);
            pass = false;
        }
    }
    if (!pass )exit(1);
}

