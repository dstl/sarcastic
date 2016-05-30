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
#include <SIlib/SIlib.h>
#include "trianglecpp.h"
#include "TriangleFile.h"
extern "C" {
#include "matrixMultiplication.h"
#include "JRStriangle.h"
#include "boxMullerRandom.h"
}
#include "materialProperties.h"
#include "colourCodes.h"
#include <sstream>
#include <forward_list>

#define ROOTPATH "/Users/Darren/Development"
#define REPORTN 100
#define VERBOSE ((int) 0)

bool VECT_EQUAL(SPVector a, SPVector b);
bool anyVECT_EQUAL(SPVector v, Triangle t);
void roughenPoint(SPVector *p, Triangle t);
bool pointsColinear(SPVector p1, SPVector p2, SPVector p3);
std::forward_list<Triangle> growTriangles(std::forward_list<Triangle> *triangles);
void trianglesToPoints(std::forward_list<Triangle> triangles, int *nTriangles, SPVector *org, SPVector *ord, SPVector *absci, int **trinds, int **segments, double **points);
void cleanPoints(int ntriangles, double **points, int **trinds, int **segments, int *nPoints, int *nSegs, double **holes, int *nholes);
SPVector  createHole(int a, int b, double **points, int c, int nPoints) ;
void PointTo2D(SPVector vec, SPVector origin, SPVector ordinate, SPVector abscissa, double *x, double *y);
bool isPointOnLine(double px, double py, double x1, double y1, double x2, double y2, double approxZero);
char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def);


int main(int argc, const char * argv[]) {
    
    int verbose = VERBOSE ;
    int nCommonTris ;
    int nPoints, nSegments ;
    SPVector origin, abscissa, ordinate,zdir;
    int *trinds;
    int *segments ;
    double *holes ;
    double *points ;
    int nholes ;

    printf(" \n");
    printf(DARK GREEN "                              Materialise\n" NORMAL);
    printf(BLUE "                          Version :" RED" %s \n", FULL_VERSION);
    printf(BLUE "                    Revision: " RED"%s, %s \n",REVISION, VERSION_DATE);
    printf(BLUE "               Copyright (c) 2013 " WHITE"[" BLUE"Dstl" WHITE"]" BLUE". All rights reserved.\n" RESETCOLOR);
    printf(" \n");

    
    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, "materialise", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status) ;
    char *instr = tryReadFile((char *)"Input triangle filename", (char *)"infilename",
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
    char *plyFname=NULL;
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
    
    std::vector<Triangle> triangles = triFile.triangles ;
    size_t nTriangles = triangles.size() ;
    
    std::forward_list<Triangle> triList ;
    for(int tt=0; tt<nTriangles; tt++){
        triList.push_front(triangles[tt]);
    }
    
    while (!triList.empty()) {
        
        std::forward_list<Triangle> commonTris ;
        int i;
        
        if(verbose){
            printf("                             Materialise the following triangles:\n");
            printf("                             ====================================\n");
            i=0;
            for (auto it=triList.begin(); it != triList.end(); ++it) {
                printf("                                          Triangle %03d \n",i++);
                printf("                                          ------------- \n");
                it->print();
            }
        }
        commonTris = growTriangles(&triList) ; // Reduces triList as triangles are added to commonTris
        if(verbose){
            i=0;
            printf("\n");
            printf("                             Triangles grown from Triangle 000 above\n");
            printf("                             ========================================\n");
            
            for (auto it=commonTris.begin(); it != commonTris.end(); ++it) {
                printf("                                          Triangle %03d \n",i++);
                printf("                                          ~~~~~~~~~~~~~ \n");
                it->print();
            }
            printf("\n");
            printf("                                 Remaining Triangles to Process\n");
            printf("                                 ==============================\n");
            
            i=0;
            for (auto it=triList.begin(); it != triList.end(); ++it) {
                printf("                                          Triangle %03d \n",i++);
                printf("                                          +++++++++++++ \n");
                it->print();
            }
            printf("                             ====================================\n\n");
        }

       
        trianglesToPoints(commonTris, &nCommonTris, &origin, &ordinate, &abscissa, &trinds, &segments, &points);
        VECT_CROSS(ordinate, abscissa, zdir);
        VECT_NORM(zdir, zdir);

        cleanPoints(nCommonTris, &points, &trinds, &segments, &nPoints, &nSegments, &holes, &nholes);
        
        if(verbose){
            printf("                     The unique set of points for this growth region are:\n");
            printf("                     ====================================================\n");
            for(int i=0; i< nPoints; i++){
                printf("%f, %f\n",points[i*2+0],points[i*2+1]);
            }
            
            if(nSegments > 0){
                printf("                            The following segments bound this region:\n");
                printf("                            =========================================\n");
                for(int i=0; i<nSegments; i++){
                    printf("                                 %05.2f,%05.2f ----> %05.2f,%05.2f\n",points[segments[i*2+0]*2+0], points[segments[i*2+0]*2+1],
                           points[segments[i*2+1]*2+0],points[segments[i*2+1]*2+1]);
                }
            }else{
                printf("\n");
                printf("                                ++++ This region has no segments +++\n");
                printf("\n");
            }
            
            if(nholes > 0){
                printf("                              The following holes are in this region:\n");
                printf("                              =======================================\n");
                for(int i=0; i<nholes; i++){
                    printf("%f,%f \n",holes[2*i+0],holes[2*i+1]);
                }
            }else{
                printf("\n");
                printf("                                ++++ This region has no holes +++\n");
                printf("\n");
            }
        }
        // find material and from it max area
        //
        double corrlen = materialProperties[commonTris.front().matId].corlen ;
        double corrArea = corrlen * corrlen * 0.5 ;
        printf("setting max area to be %f m^2 (correlation length is %f m)\n",corrArea,corrlen);
        std::ostringstream stringstream ;
        stringstream << "zpj" ;
        if (!verbose) {
            stringstream << "Q" ;
        }
        stringstream << "a" << corrArea ;
        std::string str = stringstream.str();
        const char *triswitches = str.c_str() ;
    
        // delaunay triangulate
        //
        struct triangulateio in ;
        struct triangulateio out;
        in.pointlist = points ;
        in.numberofpoints = nPoints ;
        in.numberofpointattributes = 0;
        in.pointmarkerlist = NULL ;
        in.pointattributelist = NULL;
        
        in.trianglelist = trinds ;
        in.numberoftriangles = nCommonTris ;
        in.trianglearealist = NULL;
        in.numberoftriangleattributes = 0;
        in.triangleattributelist = NULL;
        in.numberofcorners = 3;

        in.segmentlist = segments ;
        in.numberofsegments = nSegments ;
        in.segmentmarkerlist = NULL ;
        
        in.numberofholes = nholes;
        in.holelist = holes; ;
        in.numberofregions = 0;
        in.regionlist = NULL;
        
        out.pointlist = NULL;
        out.pointmarkerlist = NULL ;
        out.trianglelist = NULL ;
        out.neighborlist = NULL ;
        out.segmentlist = NULL;
        out.segmentmarkerlist = NULL;
        
        triangulate((char *)triswitches, &in, &out, NULL);
        
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
                
                if( isPointOnLine(x2d, y2d, s1x,s1y,s2x,s2y, 0.01) ){
                    pntOnSegment=true;
                }
                ++s ;
                
            }
            if (!pntOnSegment) {
                Z = box_muller(0.0, materialProperties[commonTris.front().matId].roughness) ;
            }
            
            VECT_SCMULT(ordinate, x2d, v1);
            VECT_SCMULT(abscissa, y2d, v2);
            VECT_SCMULT(zdir, Z, v3);

            VECT_ADD(v1, v2, v);
            VECT_ADD(v3, v, v);
            VECT_ADD(v, OO, v);
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
            
            Triangle newtri = Triangle(AA, BB, CC, commonTris.front().matId);
            
            if(newtri.area != 0.0){
                newTriangleVec.push_back(newtri) ;
            }
        }
        free(Points3D);
        free(points);
        free(trinds);
        free(segments);
        if(nholes !=0)free(holes);
        free(out.pointlist);
        free(out.pointattributelist);
        free(out.pointmarkerlist);
        free(out.trianglelist);
        free(out.triangleattributelist);
        free(out.trianglearealist);
        free(out.neighborlist);
        free(out.segmentlist);
        free(out.segmentmarkerlist);
        
    }
    
    // newTriangleVec now contains new triangles.
    // write it to file
    //
    printf("Writing %lu triangles to file...\n",newTriangleVec.size());
    TriangleFile dtTriFile = TriangleFile(newTriangleVec) ;
    dtTriFile.WriteFile(std::string(oustr));
    if(plyop){
        printf("Writing %lu triangles to .ply file...\n",newTriangleVec.size());
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
    double eps = 1.0e-6;
    if (fabs(a.x-b.x)<eps && fabs(a.y-b.y)<eps && fabs(a.z-b.z)<eps) {
        return true;
    }else{
        return false;
    }
}

bool anyVECT_EQUAL(SPVector v, Triangle t){
    return ( VECT_EQUAL(v, t.AA) || VECT_EQUAL(v, t.BB) || VECT_EQUAL(v, t.CC) ) ;
}

std::forward_list<Triangle> growTriangles(std::forward_list<Triangle> *triangles)
{
    Triangle subjectTri, tri, cand, commonTri;
    std::forward_list<Triangle> candidates ;
    std::forward_list<Triangle> commonTriangles ;
    std::forward_list<Triangle> moreCommons ;
    std::forward_list<Triangle> beenChecked ;
    std::forward_list<Triangle> newTriangles ;

    int nCommonVerts;
    
    subjectTri = triangles->front() ;
    commonTriangles.push_front(subjectTri) ;
    triangles->pop_front();
    
    while (!triangles->empty()){
        tri=triangles->front() ;
        triangles->pop_front();
        if ( VECT_EQUAL(tri.NN, subjectTri.NN) && (tri.matId == subjectTri.matId)) {
            candidates.push_front(tri);
        }else{
            newTriangles.push_front(tri);
        }
    }
    while (!newTriangles.empty()){
        tri=newTriangles.front();
        newTriangles.pop_front();
        triangles->push_front(tri);
    }
    
    for (auto commonIt=commonTriangles.begin(); commonIt!=commonTriangles.end(); ++commonIt) {
        commonTri = *commonIt ;
        while (!candidates.empty()) {
            cand = candidates.front();
            candidates.pop_front() ;
            nCommonVerts = 0;
            if ( anyVECT_EQUAL(commonTri.AA, cand)) nCommonVerts++ ;
            if ( anyVECT_EQUAL(commonTri.BB, cand)) nCommonVerts++ ;
            if ( anyVECT_EQUAL(commonTri.CC, cand)) nCommonVerts++ ;
            if ( nCommonVerts == 2 ){
                moreCommons.push_front(cand) ;
                commonIt = commonTriangles.before_begin() ;
            } else {
                beenChecked.push_front(cand) ;
            }
        }
        if(!moreCommons.empty()){
            while (!beenChecked.empty()) {
                tri=beenChecked.front() ;
                beenChecked.pop_front() ;
                candidates.push_front(tri) ;
            }
        }
        while (!moreCommons.empty()) {
            tri = moreCommons.front() ;
            moreCommons.pop_front() ;
            commonTriangles.push_front(tri) ;
        }
    }
    while (!beenChecked.empty()) {
        tri = beenChecked.front() ;
        beenChecked.pop_front() ;
        triangles->push_front(tri) ;
    }
    
    return (commonTriangles) ;
}


void trianglesToPoints(std::forward_list<Triangle> triangles, int *nTriangles, SPVector *org, SPVector *ord, SPVector *absci, int **trinds, int **segments, double **points){
    
    Triangle tri,tri2;
    SPVector N,ordinate,abscissa, origin;
    double x0,y0,x1,y1,x2,y2;
    int nt;

    // Work out the number of triangles
    //
    nt = 0 ;
    for( auto it = triangles.begin(); it != triangles.end(); ++it) nt++ ;
    *nTriangles = nt ;

    *points   = (double *)sp_malloc(sizeof(double) * nt * 6);
    *trinds   = (int *)sp_malloc(sizeof(int) * nt * 3);
    *segments = (int *)sp_malloc(sizeof(int) * nt * 3 * 2);

    tri = *(triangles.begin()) ;
    N = tri.NN;
    origin = tri.AA ;
    VECT_SUB(tri.BB, tri.AA, ordinate);
    VECT_NORM(ordinate, ordinate);
    VECT_CROSS(N, ordinate, abscissa);
    VECT_NORM(abscissa, abscissa);

    nt = 0;
    for( auto it = triangles.begin(); it != triangles.end(); ++it){
        tri = *it ;
        
        PointTo2D(tri.AA, origin, ordinate, abscissa, &x0, &y0) ;
        PointTo2D(tri.BB, origin, ordinate, abscissa, &x1, &y1) ;
        PointTo2D(tri.CC, origin, ordinate, abscissa, &x2, &y2) ;
        
        // put 2D point in pointlist
        //
        (*points)[nt*6 + 0] = x0; (*points)[nt*6 + 1] = y0;
        (*points)[nt*6 + 2] = x1; (*points)[nt*6 + 3] = y1;
        (*points)[nt*6 + 4] = x2; (*points)[nt*6 + 5] = y2;
        
        // Create triangle array
        //
        (*trinds)[nt * 3 + 0] = nt*3+0 ;
        (*trinds)[nt * 3 + 1] = nt*3+1 ;
        (*trinds)[nt * 3 + 2] = nt*3+2 ;
        
        // Create array of triangle sides
        //
        (*segments)[nt * 6 + 0 ] = nt*3+0 ;
        (*segments)[nt * 6 + 1 ] = nt*3+1 ;
        (*segments)[nt * 6 + 2 ] = nt*3+1 ;
        (*segments)[nt * 6 + 3 ] = nt*3+2 ;
        (*segments)[nt * 6 + 4 ] = nt*3+2 ;
        (*segments)[nt * 6 + 5 ] = nt*3+0 ;
        
        nt++;

    }
    
    *org   = origin ;
    *ord   = ordinate ;
    *absci = abscissa ;
    
    return ;
    
  
}

void cleanPoints(int ntriangles, double **points, int **trinds, int **segments, int *nPoints, int *nSegs, double **holes, int *nholes)
// Going in
//    points is ntriangles * (3 vertices) * (2 doubles per vertex)
//    trinds are the indices into points of the triangle vertices
//    segments is each side of the triangle (ntriangles * 3 sides * 2 ints per side)
// On output:
//    ntriangles is the same
//    nPoints has been modified to be number of unique points
//    points has been reduced in size to remove duplicates and is now nPoints long (nPoints * 2 * double)
//    trinds are the indices in the new points array
//    segments are now only those sides that have only one triangle connected (ie an edge)
//
{
    
    int np = ntriangles*3 ;
    int nSegments = ntriangles * 3;
    int * vert_remap = (int *)sp_malloc(sizeof(int)*np);
    for( int p=0; p<np; p++)vert_remap[p] = -1;
    
    double *uniquePoints = (double *)sp_malloc(sizeof(double) * np * 2);
    for( int p=0; p<np*2; p++) uniquePoints[p] = -666666.0 ;

    int pntsInList = 0;
    for( int p=0; p<np; p++){
        int q = 0;
        while( !((*points)[2*p+0]==uniquePoints[2*q+0] &&  (*points)[2*p+1]==uniquePoints[2*q+1]) && q<pntsInList ){
            ++q ;
        }
        if(q==pntsInList){
            uniquePoints[2*pntsInList+0] = (*points)[2*p+0] ;
            uniquePoints[2*pntsInList+1] = (*points)[2*p+1] ;
            ++pntsInList ;
        }
        vert_remap[p] = q ;
        
    }
    
    *nPoints = pntsInList ;
    free(*points) ;
    *points = (double *)sp_malloc(sizeof(double) * 2 * *nPoints);
    for(int p=0; p< *nPoints*2; ++p) {
        (*points)[p] = uniquePoints[p] ;
    }
    free(uniquePoints) ;

    for( int p=0; p<ntriangles*3; p++){
        if( vert_remap[(*trinds)[p]] != -1 ){
            (*trinds)[p] = vert_remap[(*trinds)[p]] ;
        }
    }
    for( int p=0; p<nSegments*2; p++){
        if( vert_remap[(*segments)[p]] != -1 ){
            (*segments)[p] = vert_remap[(*segments)[p]] ;
        }
    }
    
    free(vert_remap);

    // Now find segments that only have one triangle connected - these are edges or true segments
    //
    int segA, segB, segAA, segBB ;
    std::forward_list<int> seglist ;
    std::forward_list<double> hole_list ;

    int cnt;
    for( int p=0; p<nSegments; p++){
        cnt = 0;
        segA = (*segments)[2*p+0] ;
        segB = (*segments)[2*p+1] ;
        for( int q=0; q<nSegments; q++){
            if( q!=p){
                segAA = (*segments)[2*q+0] ;
                segBB = (*segments)[2*q+1] ;
                
                if( (segAA == segA && segBB == segB) || (segAA == segB && segBB == segA) ){
                    cnt++;
                }
            }
        }
        if (cnt == 0) {
            // Egde found
            seglist.push_front(segA);
            seglist.push_front(segB);
            
            // find triangle owning this segment
            //
            for (int t=0; t<ntriangles; ++t) {
                int a,b,c ;
                a = (*trinds)[t*3+0] ;
                b = (*trinds)[t*3+1] ;
                c = (*trinds)[t*3+2] ;
                
                SPVector hole ;
                hole.z = -666.0;
                if ( (segA == a && segB == b) ){
                    hole = createHole(a,b,points,c, *nPoints) ;
                    
                }else if ( segA == b && segB == c){
                    hole = createHole(b,c,points,a, *nPoints) ;
                    
                }else if( segA == c && segB == a){
                    hole = createHole(c,a,points,b, *nPoints) ;
                    
                }else if( segA == a && segB == c){
                    hole = createHole(a,c,points,b, *nPoints) ;
                    
                }else if( segA == c && segB == b){
                    hole = createHole(c,b,points,a, *nPoints) ;
                    
                }else if( segA == b && segB == a){
                    hole = createHole(b,a,points,c, *nPoints) ;
                    
                }
                if (hole.z != -666 ) {
                    hole_list.push_front(hole.x);
                    hole_list.push_front(hole.y);
                }
            }
        }
    }
    seglist.reverse() ;
    cnt = 0;
    for( auto it = seglist.begin(); it != seglist.end(); ++it){
        cnt++ ;
    }
    *nSegs = cnt / 2;
    
    free(*segments);
    *segments = (int *)sp_malloc(sizeof(int) * 2 * (*nSegs)) ;
    cnt=0;
    while (!seglist.empty()){
        (*segments)[cnt++] = seglist.front() ;
        seglist.pop_front() ;
    }
    
    hole_list.reverse() ;
    cnt = 0;
    for( auto it = hole_list.begin(); it != hole_list.end(); ++it){
        cnt++;
    }
    *nholes = cnt / 2 ;
    
    if(*nholes != 0){
//        if((*holes) != NULL )free(*holes) ;
        *holes = (double *)sp_malloc(sizeof(double) * 2 * (*nholes));
        cnt=0;
        while (!hole_list.empty()){
            double d = hole_list.front() ;
            (*holes)[cnt++] = d;
            hole_list.pop_front() ;
        }
    }
    
    return ;
    
}
    

SPVector  createHole(int a, int b, double **points, int c, int nPoints)
// Creates a hole coordinate outside of the triangle specified by abc on the
// side of ab
//
{
    SPVector Aa,Bb,Cc,v1,xp,hole;
    SPVector v2,v3,AB,xpc,outdir,indir;

    double dist;
    double epsi = 0.001 ; // small distance to create a hole point over a line segment
   
    
    
    VECT_CREATE((*points)[a*2+0], (*points)[a*2+1], 0.0, Aa);
    VECT_CREATE((*points)[b*2+0], (*points)[b*2+1], 0.0, Bb);
    VECT_CREATE((*points)[c*2+0], (*points)[c*2+1], 0.0, Cc);
    VECT_SUB(Bb, Aa, AB);
    VECT_SCMULT(AB, 0.5, AB);
    VECT_ADD(Aa,AB,xp);
    VECT_SUB(Cc,xp,xpc);
    VECT_CROSS(AB, xpc, v1);
    VECT_CROSS(AB, v1, outdir);
    VECT_NORM(outdir, outdir);
    VECT_MINUS(outdir, indir);
    VECT_SCMULT(outdir, 0.01, v1);
    
    bool pointFurther = false;
    int p=0;
    while (p < nPoints && pointFurther == false){
        // Only make a hole if there is another point further away from cp (central point on segment
        //
        VECT_CREATE((*points)[p*2+0], (*points)[p*2+1], 0.0, v2);
        VECT_SUB(v2, xp, v3);
        dist = VECT_DOT(v3, outdir);
        if (dist > epsi) pointFurther = true ;
        ++p;
    }
    if(pointFurther){
        VECT_ADD(xp, v1, hole);
    }else{
        VECT_CREATE(0,0, -666, hole);
    }
    
    return hole;
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
// approx zero is a small number. if teh distance of the point p from the line is less than this then
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
    
    return(fname) ;
}

