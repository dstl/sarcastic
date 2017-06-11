#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <vector>
#include <CGAL/Delaunay_mesher_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_2                                              Point_2;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                    Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>            Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>      CDT;
typedef CDT::Vertex_handle                                      Vertex_handle;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                Criteria;
typedef CDT::Vertex_iterator                                    Vertex_iterator;
typedef CDT::Vertex                                             Vertex;
typedef CDT::Face                                               Face ;
typedef Face::Face_handle                                       Face_handle ;
typedef CDT::Face_circulator                                    Face_circulator ;


int main(int argc, char* argv[])
{
    
    // Create a vector of the points
    //
    std::vector<Point_2> points2D ;
    points2D.push_back(Point_2(-5.0, -5.0));    // ----------
    points2D.push_back(Point_2( 5.0, -5.0));    // |   OUTER
    points2D.push_back(Point_2( 5.0,  5.0));    // |   SQUARE
    points2D.push_back(Point_2(-5.0,  5.0));    // ----------
    points2D.push_back(Point_2(-2.5, -2.5));    // ----------
    points2D.push_back(Point_2( 2.5, -2.5));    // |   INNER
    points2D.push_back(Point_2( 2.5,  2.5));    // |   SQUARE
    points2D.push_back(Point_2(-2.5,  2.5));    // ----------
    size_t numTestPoints = points2D.size();
    
    // Create a constrained delaunay triangulation and add the points
    //
    CDT cdt;
    std::vector<Vertex_handle> vhs;
    for (unsigned int i=0; i<numTestPoints; ++i){
        vhs.push_back(cdt.insert(points2D[i]));
    }
    
    // Creare constraints of the sides of the mesh
    //
    cdt.insert_constraint(vhs[0],vhs[1]);
    cdt.insert_constraint(vhs[1],vhs[2]);
    cdt.insert_constraint(vhs[2],vhs[3]);
    cdt.insert_constraint(vhs[3],vhs[0]);
    
    cdt.insert_constraint(vhs[4],vhs[5]);
    cdt.insert_constraint(vhs[5],vhs[6]);
    cdt.insert_constraint(vhs[6],vhs[7]);
    cdt.insert_constraint(vhs[7],vhs[4]);
    
    // Create a seed to make sure the inner square is a hole
    //
    std::list<Point_2> list_of_seeds;
    list_of_seeds.push_back(Point_2(0, 0));
    
    // Refine the mesh into triangles of max side length '1' and ensuring seeds are 'holes'
    //
    CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(), list_of_seeds.end(),
                                 Criteria(0.125, 1.5),false);
    
    // The mesh is now created. The next bit swings around each vertex point checking that
    // all faces around it are in the domain. If any are not then the vertex is on the
    // edge of the mesh.
    // thanks to @sloriot for this
    //
    
    Vertex v ;
    std::vector<Point_2> interior_points ;
    std::vector<Point_2> boundary_points ;
    CDT::Locate_type loc = CDT::Locate_type::VERTEX ;
    int li;
    
    Vertex_iterator vit = cdt.vertices_begin(), beyond = cdt.vertices_end() ;
    while (vit != beyond) {
        v = *vit ;
        Point_2 query = vit->point();
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
            interior_points.push_back(query) ;
        }else{
            boundary_points.push_back(query) ;
        }
        ++vit ;
    }
    
    printf("Boundary points:\n");
    for(auto p = boundary_points.begin(); p != boundary_points.end(); ++p){
        printf("%f,%f\n",p->x(), p->y()) ;
    }

    printf("interior points:\n");
    for(auto p = interior_points.begin(); p != interior_points.end(); ++p){
        printf("%f,%f\n",p->x(), p->y()) ;
    }
    
    return 0 ;
}
