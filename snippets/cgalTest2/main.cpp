#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
//#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <vector>
//#include <CGAL/Delaunay_mesher_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef K::Point_2                                              Point_2;
typedef CGAL::Delaunay_mesh_vertex_base_2<K>                    Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K>                      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>            Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>      CDT;
typedef CDT::Vertex_handle                                      Vertex_handle;
//typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                Criteria;
typedef CDT::Face_iterator                                      Face_iterator;
//typedef CDT::Vertex                                             Vertex;
//typedef CDT::Face                                               Face ;
//typedef Face::Face_handle                                       Face_handle ;
//typedef CDT::Face_circulator                                    Face_circulator ;


int main(int argc, char* argv[])
{
    
    // Create a vector of the points
    //
    std::vector<Point_2> points2D ;
    points2D.push_back(Point_2(  1,  1));
    points2D.push_back(Point_2(  1, -1));
    points2D.push_back(Point_2( -1,  1));
    points2D.push_back(Point_2( -1, -1));
    points2D.push_back(Point_2( 0, 0));

    size_t numTestPoints = points2D.size();
    
    // Create a constrained delaunay triangulation and add the points
    //
    CDT cdt;
    std::vector<Vertex_handle> vhs;
    for (unsigned int i=0; i<numTestPoints; ++i){
        vhs.push_back(cdt.insert(points2D[i]));
    }
    
    int i=0;
    for (Face_iterator fit = cdt.faces_begin()  ; fit != cdt.faces_end(); ++fit) {
        printf("Face %d is (%f,%f) -- (%f,%f) -- (%f,%f) \n",i++,
               fit->vertex(0)->point().x(),fit->vertex(0)->point().y(),
               fit->vertex(1)->point().x(),fit->vertex(1)->point().y(),
               fit->vertex(2)->point().x(),fit->vertex(2)->point().y() );
        
    }
    

    return 0 ;
}
