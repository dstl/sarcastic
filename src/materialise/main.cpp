/***************************************************************************
 * 
 *           Module :  main.cpp
 *          Program :  materialise
 *       Created by :  Darren Muff on 11/06/2017
 *   CLASSIFICATION :  Official
 *   Date of CLASSN :  05-Nov-2018
 *   Description:
 *      This programme reads in a .ply file and breaks it down into Delaunay
 *      triangles. The input .ply file can be a mesh or a triangle soup. In the 
 *      latter case (where vertex points can be unique to a triangle) the 
 *      programme searches for them and splits the soup up so that each triangle 
 *      has only onle pertner for each of its edges. The final trianglation has 
 *      a roughness applied that is determined by the input 'material' 
 *      characteristics and the surface roughness parameters specified in the 
 *      materialProperties.h header file. if any of the surfaces of the model
 *      has a surface roughness that is not zero then the triangulation will 
 *      perform a post process pass to randomly 'roughen' the surface such that
 *      the standard deviation is equivalent to the specified surface roughness.
 *
 *      The progremme makes extensive use of teh CGAL library [1] and in 
 *      particular the Delaunay triangulation routines [2]
 *
 *          [1] CGAL, Computational Geometry Algorithms Library, http://www.cgal.org
 *          [2] Laurent Rineau. 2D Conforming Triangulations and Meshes.
 *              In CGAL User and Reference Manual. CGAL Editorial Board, 4.10 edition, 2017
 *
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

#define REAL double
#define ANSI_DECLARATORS

#include "materialise_version.h"
#include <iostream>
#include <sarclib/sarclib.h>
#include "TriangleMesh.hpp"
extern "C" {
#include "matrixMultiplication.h"
#include "boxMullerRandom.h"
#include "printProgress.h"
}
#include <boost/filesystem.hpp>
#include <vector>
#include "readMaterialFile.hpp"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>

typedef int Index ;
typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef CGAL::Triangulation_vertex_base_with_info_2<Index, K>   Vbb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                      Fb;
typedef CGAL::Triangulation_data_structure_2<Vbb,Fb>            Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>      CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                Criteria;
typedef CDT::Vertex_handle                                      Vertex_handle;
typedef CDT::Vertex_iterator                                    Vertex_iterator;
typedef CDT::Vertex                                             Vertex;
typedef CDT::Face                                               Face ;
typedef CDT::Face_circulator                                    Face_circulator ;
typedef CDT::Edge                                               Edge;
typedef CDT::Edge_iterator                                      Edge_iterator;
typedef Face::Face_handle                                       Face_handle ;
typedef K::Point_2                                              Point_2;

//#define DEBUGCOUNT 1529

#define REPORTN 100
#define VERBOSE ((int) 0)

#define RESETCOLOR  "\033[0m"
#define RED         "\033[31m"
#define GREEN       "\033[32m"
#define YELLOW      "\033[33m"
#define BLUE        "\033[34m"
#define WHITE       "\033[37m"
#define NORMAL      "\033[22m"
#define DARK        "\033[1m"

using namespace std;

void roughenPoint(SPVector *p, Triangle t);
TriangleMesh growTriangles(TriangleMesh *mesh) ;
void meshTo2D(TriangleMesh *mesh, SPVector *org, SPVector *absci, SPVector *ord,
              int **trinds,double **segments,int *nsegments, double **points, int *nholes, double **holes) ;
void cleanPoints(int ntriangles, double **points, int **trinds, int **segments, int *nPoints, int *nSegs, double **holes, int *nholes);
void PointTo2D(SPVector vec, SPVector origin, SPVector abscissa, SPVector ordinate, double *x, double *y);
bool isPointOnLine(double px, double py, double x1, double y1, double x2, double y2, double approxZero);
char * tryReadFile(const char *prompt, const char * key, const char * help, const char *def);
bool pointInPolygon(int nsegments, double *segments, double x, double y);
bool pointInPolygon(int nsegments, double *segments, double x, double y, int *nleft, int *nright);
void testcommonmesh(TriangleMesh *mesh) ;
bool edgesIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double &xpt, double &ypt) ;

class Line {
public:
    double p1x ;
    double p1y ;
    double p2x ;
    double p2y ;
    
    Line(){} ;
    Line(double x1, double y1, double x2, double y2): p1x(x1), p1y(y1), p2x(x2), p2y(y2) {} ;

};

int main(int argc, const char * argv[]) {
    
    int verbose = VERBOSE ;
    int nSegments, nholes ;
    SPVector origin, abscissa, ordinate, zdir;
    int *trinds;
    double *segments ;
    double *holes ;
    double *points ;

    printf(" \n");
    printf(DARK GREEN "                              Materialise\n" NORMAL);
    printf(BLUE "                          Version :" RED" %s \n", MATERIALISE_FULL_VERSION);
    printf(BLUE "                    Revision: " RED"%s, %s \n",MATERIALISE_REVISION, MATERIALISE_VERSION_DATE);
    printf(BLUE "               Copyright (c) 2017 " WHITE"[" BLUE"Dstl" WHITE"]" BLUE". All rights reserved.\n" RESETCOLOR);
    printf(" \n");

    SPStatus status ;
    im_init_status(status, 0);
    im_init_lib(&status, "materialise", argc, (char **)argv);
    CHECK_STATUS_NON_PTR(status) ;
    char *instr = tryReadFile((char *)"Input triangle filename", (char *)"infilename",
                               (char *)"The name of a triangle file created with 'colladaToTriFile'",
                               (char *) "triangles.ply");
    
    char *oustr = input_string((char *)"Output triangle filename", (char *)"oufilename",
                               (char *)"The name of the triangle file to be created",
                               (char *) "delaunay.ply");
    
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
    
    // Read in the material properties file if required
    //
    char *matfile = input_string((char *)"Input materialfile filename", (char *)"materialfilename",
                               (char *)"The name of a 'materialfile' or 'none' (defaults used)",
                               (char *) MATERIALPROPS);
    initialiseMaterials(matfile, true);
    
    TriangleMesh mesh, newMesh;
    mesh.readPLYFile(string(instr));
    mesh.checkIntegrityAndRepair();
    mesh.sortTrianglesAndPoints() ;
    long int orgNTris;
    float percentDone;
    orgNTris = mesh.triangles.size() ;
    Timer timer;
    startTimer(&timer, &status);
    int cnt=0;
    
    // set up some timers
    //
    Timer grow_t,meshto2D_t, triangulate_t, roughen_t, addTri_t;
    double growTimer = 0;
    double mesh2DTimer = 0;
    double trianglateTimer = 0;
    double addTriTimer = 0;
    double roughenTimer = 0;
    
    while (!mesh.triangles.empty()) {
        
        if(cnt % 100 == 0){
            percentDone = 100.0 * (float)(orgNTris - mesh.triangles.size()) / (float)orgNTris ;
            printProgress(percentDone, 80);
        }
        
        for(int kk=0; kk<mesh.triangles.size(); ++kk){
            if(mesh.triangles[kk].a == mesh.triangles[kk].b || mesh.triangles[kk].a == mesh.triangles[kk].c || mesh.triangles[kk].b == mesh.triangles[kk].c){
                printf("DEGENERATE TRIANGLE at %d/%d for mesh %d\n",kk,(int)mesh.triangles.size(),cnt);
            }
        }
#ifdef DEBUGCOUNT
        if(cnt == DEBUGCOUNT){
            printf("DEBUG Break \n");
        }
#endif
        
        // Take the first triangle in mesh.triangles and calculate all triangles that
        // are common to it (on the same plane and joined together
        //
        startTimer(&grow_t, &status);
        TriangleMesh commonMesh = growTriangles(&mesh) ;
        endTimer(&grow_t, &status);
        growTimer += timeElapsedInMilliseconds(&grow_t, &status);
      
#ifdef DEBUGCOUNT
        if(cnt == DEBUGCOUNT){
            commonMesh.writePLYFile("debug.ply");
        }
#endif
        // Convert to 2Dpoints
        //
        startTimer(&meshto2D_t, &status);
        meshTo2D(&commonMesh, &origin, &abscissa, &ordinate, &trinds, &segments, &nSegments, &points, &nholes, &holes);
        endTimer(&meshto2D_t, &status);
        mesh2DTimer += timeElapsedInMilliseconds(&meshto2D_t, &status);
        VECT_CROSS(abscissa, ordinate, zdir);
        VECT_NORM(zdir, zdir);

#ifdef DEBUGCOUNT
        if(cnt == DEBUGCOUNT){
            for(auto t=commonMesh.triangles.begin(); t != commonMesh.triangles.end(); t++){
                printf("Triangle %ld\n",t-commonMesh.triangles.begin());
                printf("%f,%f\n",points[2*t->a+0],points[2*t->a+1]);
                printf("%f,%f\n",points[2*t->b+0],points[2*t->b+1]);
                printf("%f,%f\n",points[2*t->c+0],points[2*t->c+1]);

            }
        }
#endif
        
        // find material and from it max area
        //
        double corrlen = globalMatProps[commonMesh.triangles[0].mat].corlen ;
        double corrArea = corrlen * corrlen * 0.5 ;
        if(verbose)printf("setting max area to be %f m^2 (correlation length is %f m)\n",corrArea,corrlen);
        
        // Perform Delaunay triangulation using CGAL
        //
        CDT cdt;
        std::vector<Point_2> dtpoints;
        for(int n=0; n<commonMesh.vertices.size(); ++n){
            dtpoints.push_back(Point_2(points[n*2+0], points[n*2+1])) ;
            
#ifdef DEBUGCOUNT
            if(cnt == DEBUGCOUNT){
                printf("%f, %f\n",points[n*2+0],points[n*2+1]) ;
            }
#endif

        }
        cdt.insert(dtpoints.begin(),dtpoints.end());

        Vertex_handle vha, vhb;
        for (int n=0; n<nSegments; ++n){
            double ax,ay,bx,by;
           
            ax = segments[n*4+0] ;
            ay = segments[n*4+1] ;
            bx = segments[n*4+2] ;
            by = segments[n*4+3] ;
            vha = cdt.insert(Point_2(ax, ay));
            vhb = cdt.insert(Point_2(bx, by));
#ifdef DEBUGCOUNT
            if(cnt == DEBUGCOUNT){
                printf("segment %d\n",n);
                printf("%f,%f\n",ax,ay);
                printf("%f,%f\n",bx,by);
            }
#endif

            cdt.insert_constraint(vha, vhb);
        }
        
        std::list<Point_2> list_of_seeds;
        for(int n=0; n<nholes; ++n){
            list_of_seeds.push_back(Point_2(holes[2*n+0], holes[2*n+1]));
        }
        startTimer(&triangulate_t, &status);
        CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                                     Criteria(0.125, corrlen),false);
        endTimer(&triangulate_t, &status);
        trianglateTimer += timeElapsedInMilliseconds(&triangulate_t, &status);
 
        startTimer(&roughen_t, &status);
        Vertex v ;
        CDT::Locate_type loc = CDT::Locate_type::VERTEX ;
        int li;
        double x2d,y2d, Z;
        SPVector v0,v1,v2,v3,OO ;
        OO = origin;
        int vind = 0;
        SPVector *Points3D = (SPVector *)sp_malloc(sizeof(SPVector) * cdt.number_of_vertices() ) ;
        for(Vertex_iterator vit = cdt.vertices_begin(); vit !=cdt.vertices_end() ; vit++){
            
            v = *vit ;
            vit->info() = vind ;
            Point_2 query = vit->point();
            x2d = query.x() ;
            y2d = query.y() ;
            Z = 0 ;
            Face_handle f =  cdt.locate(query, loc, li);
            
            bool current_face_in_domain ;
            Face_handle start = f ;
            do {
                current_face_in_domain = f->is_in_domain() ;
                Vertex_handle vh = f->vertex(li);
                f = f->neighbor(cdt.ccw(li));
                li = f->index(vh) ;
            }
            while(current_face_in_domain && f != start) ;
            
            if (f == start) {
                // Interior point found - not on the boundary
                //
                Z = box_muller(0.0, globalMatProps[commonMesh.triangles[0].mat].roughness) ;
            }
            
            VECT_SCMULT(abscissa, x2d, v1);
            VECT_SCMULT(ordinate, y2d, v2);
            VECT_SCMULT(zdir, Z, v3);
            
            VECT_ADD(v1, v2, v0);
            VECT_ADD(v3, v0, v0);
            VECT_ADD(v0, OO, v0);
            
            Points3D[vind++] = v0 ;
        }
        endTimer(&roughen_t, &status);
        roughenTimer += timeElapsedInMilliseconds(&roughen_t, &status);

        startTimer(&addTri_t, &status);
        SPVector triVerts[3] ;
        for (CDT::Finite_faces_iterator fit=cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit){
            if(fit->is_in_domain()){
                for (unsigned int vi = 0; vi<3; ++vi) {
                    Vertex_handle vh = fit->vertex(vi) ;
                    int pnt = vh->info() ;
                    triVerts[vi] = Points3D[pnt] ;
                }
            newMesh.addTriangle(triVerts[0],triVerts[1],triVerts[2],commonMesh.triangles[0].mat);
            }
        }

        endTimer(&addTri_t, &status);
        addTriTimer += timeElapsedInMilliseconds(&addTri_t, &status);
        
        free(Points3D);
        free(points);
        free(trinds);
        free(segments);
        if(nholes!=0)free(holes);
        
        // Next triangle
        //
        cnt++;
    }
    
 
    printf("\nWriting %lu triangles to file %s...\n",newMesh.triangles.size(),oustr);
    newMesh.writePLYFile(std::string(oustr)) ;
    endTimer(&timer, &status);
    
    printf("           Timing Summary\n");
    printf("===============================================\n");
    printf("Time spent growing triangles      : %8.2f ms\n",growTimer);
    printf("Time spent converting tris to 2D  : %8.2f ms\n",mesh2DTimer);
    printf("Delaunay triangulation            : %8.2f ms\n",trianglateTimer);
    printf("Time spent roughening surface     : %8.2f ms\n",roughenTimer);
    printf("Time spent adding triangles       : %8.2f ms\n",addTriTimer);
    printf("-----------------------------------------------\n");
    printf("  Total time            : %8.2f ms\n", timeElapsedInMilliseconds(&timer, &status));

    im_close_lib(&status);
    std::cout << "Done in " << timeElapsedInSeconds(&timer, &status) << " seconds." << std::endl;
    
    delete [] globalMatProps ;
    delete [] globalMatColours ;
    free(matfile);
    free(instr);
    free(oustr);

    return 0;
}


TriangleMesh growTriangles(TriangleMesh *mesh)
/// Takes the first triangle from the Triangles vector<> within the input mesh and finds all the coplanar triangles
/// inside the triangles vector<>. These are removed from the input mesh and returned as another mesh
/// The function 'grows' the coplanar region around the first triangle of the input mesh so that
/// other faces in the maesh that have the same normal but are not attached to the first region
/// are not considered. Finally the function checks for 'polygamous' triangles within the input mesh
/// these are triangles that have more that one adjoining partner triangle along one of its edges. If
/// any are found the polygamous triangle is cut into two at the offending vertex to ensure that the output is a true
/// 'mesh' rather than a 'triangle soup'
///
{

    Triangle tri;
    Triangle objectTri(mesh->triangles.back());
    
    // Copy mesh triangles to originals - we will read back the ones that don't match afterwards
    //
    vector<Triangle> originals;
    while(mesh->triangles.size() > 0){
        originals.push_back(mesh->triangles.back());
        mesh->triangles.pop_back();
    }
    
    // Find all the candidates and throw back the ones that are not co-planar with the object tri
    //
    vector<Triangle> candidates;
    for(int i=(int)originals.size()-1; i>=0; --i){
        tri = originals[i] ;
        if ((tri.mat == objectTri.mat) && (fabs(tri.dist-objectTri.dist) < 1e-07) && (tri.N == objectTri.N) ) {
            candidates.push_back(tri) ;
        }else{
            mesh->triangles.push_back(tri) ;
        }
    }
    
    // loop through all the candidates and see if they adjoin any of those in the growing commonMesh
    //
    vector<Triangle> commonTriangles, pending;
    tri = candidates.back() ;
    candidates.pop_back() ;
    commonTriangles.push_back(tri) ;
    bool IGROO ;
    do{
        IGROO = false ;
        do{
            if(candidates.size() > 0){
                tri = candidates.back() ;
                candidates.pop_back() ;
            }else{
                tri.a = tri.b = tri.c = -1 ;
            }
            
            std::vector<int> commonVerts;
            commonVerts.reserve(commonTriangles.size() * 3) ;
            for(int t=0; t < commonTriangles.size(); ++t){
                commonVerts.push_back(commonTriangles[t].a);
                commonVerts.push_back(commonTriangles[t].b);
                commonVerts.push_back(commonTriangles[t].c);
            }
            sort(commonVerts.begin(), commonVerts.end()) ;
            std::vector<int>::iterator afind = std::find(commonVerts.begin(),commonVerts.end(), tri.a);
            std::vector<int>::iterator bfind = std::find(commonVerts.begin(),commonVerts.end(), tri.b);
            std::vector<int>::iterator cfind = std::find(commonVerts.begin(),commonVerts.end(), tri.c);
            if( afind == commonVerts.end() && bfind == commonVerts.end() && cfind == commonVerts.end()){
                // tri doesn't yet touch commonTriangles - next triangle
                //
                if(tri.a != -1 && tri.b != -1 && tri.c != -1)
                    pending.push_back(tri);
            }else{
                // tri is touching commonTriangles
                //
                commonTriangles.push_back(tri) ;
                IGROO = true;
                for(int i = (int)pending.size()-1; i >= 0;  --i){
                    candidates.push_back(pending[i]) ;
                }
                pending.clear() ;
                break ;
            }
            
        }while(candidates.size() > 0);
        
        // Some triangles may be coplanar but never touch the growing region
        // we need to throw these back into the original mesh to be reconsidered
        //
        while(pending.size() > 0){
            mesh->triangles.push_back(pending.back()) ;
            pending.pop_back() ;
        }
        
    }while(IGROO == true );

    TriangleMesh newMesh(commonTriangles, mesh->vertices);
    newMesh.monogamise() ;
    return newMesh ;

}

void meshTo2D(TriangleMesh *mesh, SPVector *org, SPVector *absci, SPVector *ordi,
              int **trinds,double **segments,int *nsegments, double **points, int *nholes, double **holes)
/// Function takes in a triangleMesh ptr. The triangles in the mesh must all be coplanar. ie, each
/// triangle has the same normal. It then converts the mesh into 2D points in terms of x and y and
/// returns the abscissa (x) and ordinate (y) in 'absci' and 'ordi'.
///
/// The function then returns the mesh as 2D pairs in the array 'points' ordered as x0,y0,x1,y1,...xn,yn
/// trinds is anarray that is filled with the indices of the vertex points of each of the input triangles
/// (indexed into the array points) such that triangle #3 with vertices A,B & C has vertices given by:
///  Vertex 'A' ->    x = points[trinds[3*3+0] * 2 + 0]  y = points[trinds[3*3+0] * 2 + 1]
///  Vertex 'B' ->    x = points[trinds[3*3+1] * 2 + 0]  y = points[trinds[3*3+1] * 2 + 1]
///  Vertex 'C' ->    x = points[trinds[3*3+2] * 2 + 0]  y = points[trinds[3*3+2] * 2 + 1]
///
/// The function also returns an array containing the line segments that bound the edge of the mesh
/// This is calculated by finding all the triangle sides in the mesh that do not have a neighbor
/// triangle tht shares the same egde. It is worth noting that for this to work the input mesh muct
/// truly be a 'mesh' where each vertex point shares multiple triangles rather than a 'triangle soup'
/// where a vertex can be unique to a triangle. Segments are stored as 4 doubles with the first two being
/// x0,y0 and the seceond being x1,y1 which are the x,y coordinates of each segment
///
/// Finally any holes within the mesh are claculated and teh x,y coordinates or each hole returned in
/// the array 'holes' with nholes being the number of holes found.
///
/// Is is important that the calling function frees up memory alocated by this function for
///   *point, *trinds, *segments and *holes
///
{
    
    Triangle tri,tri2;
    SPVector N,abscissa,ordinate, origin;
    double x0,y0,x1,y1;
    
    if(!mesh->isMonogamous()){
        printf("Error : this mesh may be a 'triangle soup' with vertciesthat are unique to a triangle\n");
        printf("Error : This can cause problems when calculating the mesh boundary segments.\n");
        printf("Error : Make sure that TriangleMesh.monogamise() is run on the mesh before calling\n");
        printf("Error : this function\n");
        exit(1);
    }
    
    
    *points   = (double *)sp_malloc(sizeof(double) * mesh->vertices.size()*2);
    *trinds   = (int *)sp_malloc(sizeof(int) * mesh->triangles.size() * 3);

    tri = *(mesh->triangles.begin()) ;
    N = tri.N.asSPVector();
    origin = mesh->vertices[tri.a].asSPVector() ;
    VECT_SUB(mesh->vertices[tri.b].asSPVector(), mesh->vertices[tri.a].asSPVector(), abscissa);
    VECT_NORM(abscissa, abscissa);
    VECT_CROSS(N, abscissa, ordinate);
    VECT_NORM(ordinate, ordinate);

    // Convert vertices to 2D points
    //
    for (int iv=0; iv < mesh->vertices.size(); ++iv){
        PointTo2D(mesh->vertices[iv].asSPVector(), origin, abscissa, ordinate, &x0, &y0);
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
    *segments = (double *)sp_malloc(sizeof(double) * edges.size() * 4);
    for(int ie=0; ie<edges.size(); ++ie){
        PointTo2D(mesh->vertices[edges[ie].vertex].asSPVector(), origin, abscissa, ordinate, &x0, &y0);
        PointTo2D(mesh->vertices[edges[ie].nextVertex].asSPVector(), origin, abscissa, ordinate, &x1, &y1);

        (*segments)[4*ie+0] = x0 ;
        (*segments)[4*ie+1] = y0 ;
        (*segments)[4*ie+2] = x1 ;
        (*segments)[4*ie+3] = y1 ;
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

        withinPoly = pointInPolygon(*nsegments, *segments, dummyHole.x, dummyHole.y, &nleft, &nright) ;

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
    *absci = abscissa ;
    *ordi  = ordinate ;
    
    return ;
}

void PointTo2D(SPVector vec, SPVector origin, SPVector abscissa, SPVector ordinate, double *x, double *y){

    SPVector v ;
    
    VECT_SUB(vec, origin, v);
    *x = VECT_DOT(v, abscissa) ;
    *y = VECT_DOT(v, ordinate) ;

    return;
}

bool isPointOnLine(double px, double py, double x1, double y1, double x2, double y2, double approxZero)
/// Line is defined as the segment between points x1,y1 and x2,y2
/// point px,py is the point being tested
/// approx zero is a small number. if the distance of the point p from the line is less than this then
/// it is assumed the point is on the line
///
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
        
        std::string filename(fname);
        string extension = boost::filesystem::extension(filename);
        if(extension != ".ply"){
            printf("Input file does not have \'.ply\' extension %s\n",fname);
            fileStat.status = BAD_FILE ;
        }
        
        if ( (fp = fopen(fname, "r")) == NULL){
            printf("Cannot access file %s\n",fname);
            fileStat.status = BAD_FILE ;
        }
        
    } while(fileStat.status != NO_ERROR);
    
    fclose(fp);
    free(prmpt) ;
    return(fname) ;
}

bool pointInPolygon(int nsegments, double *segments, double x, double y)
/// determines whether point x,y is inside the polygon described by the array
/// of segment edges provided in 'segments'.
/// 'segments' is an array of length 'nsegments' containing for doubles for each
/// segment corresponding to the x0,y0 and x1,y1 coordinates of the start and ends of
/// the segment.
///
{
    bool  oddNodes=false   ;
    
    for(int i=0; i<nsegments; i++){
        double psx, psy, pex, pey;
        psx = segments[4*i+0] ;
        psy = segments[4*i+1] ;
        pex = segments[4*i+2] ;
        pey = segments[4*i+3] ;
    
        if( ((psy < y && pey >= y) || (pey < y && psy >= y )) && (psx <=x || pex <=x ) ){
            oddNodes ^= ((psx + ((y-psy) / (pey - psy )) * (pex - psx )) < x) ;
        }
        
    }
    return oddNodes;
}

/// determines whether point x,y is inside the polygon described by the array
/// of segment edges provided in 'segments'.
/// 'segments' is an array of length 'nsegments' containing for doubles for each
/// segment corresponding to the x0,y0 and x1,y1 coordinates of the start and ends of
/// the segment.
/// This function also returns the number of segments to the left and right of the test point
///
bool pointInPolygon(int nsegments, double *segments, double x, double y, int *nleft, int *nright)
{
    bool  oddNodes=false   ;
    int l=0,r=0;
    
    for(int i=0; i<nsegments; i++){
        double psx, psy, pex, pey;
        psx = segments[4*i+0] ;
        psy = segments[4*i+1] ;
        pex = segments[4*i+2] ;
        pey = segments[4*i+3] ;
        
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

bool edgesIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double &xpt, double &ypt) {
    double A1,A2,B1,B2,C1,C2,det,x,y ;
    A1  = y2-y1 ;
    B1  = x2-x1 ;
    C1  = A1*x1 + B1*y1 ;
    A2  = y4-y3 ;
    B2  = x4-x3 ;
    C2  = A2*x3 + B2*y3 ;
    det = A1*B2 - A2*B1 ;
    if(det == 0){
        // LInes are parallel
        //
        return false ;
    }
    x = (B2*C1 - B1*C2) / det ;
    y = (A1*C2 - A2*C1) / det ;
    
    double line1minx, line1maxx, line2minx, line2maxx ;
    double line1miny, line1maxy, line2miny, line2maxy ;
    double bothminx, bothmaxx, bothminy, bothmaxy ;
    double Ep = 1e-6 ;
    line1minx = std::min(x1,x2); line1maxx = std::max(x1,x2) ;
    line2minx = std::min(x3,x4); line2maxx = std::max(x3,x4) ;
    line1miny = std::min(y1,y2); line1maxy = std::max(y1,y2) ;
    line2miny = std::min(y3,y4); line2maxy = std::max(y3,y4) ;
    
    bothminx = std::max(line1minx,line2minx) ;
    bothmaxx = std::min(line1maxx,line2maxx) ;
    bothminy = std::max(line1miny,line2miny) ;
    bothmaxy = std::min(line1maxy,line2maxy) ;
    
    if( (x > bothminx - Ep) && (x < bothmaxx + Ep) && (y > bothminy - Ep) && (y < bothmaxy + Ep) ){
        xpt = x;
        ypt = y;
        return true;
    }
    
    return false ;
}


