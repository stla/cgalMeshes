#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
namespace SMP = CGAL::Surface_mesh_parameterization;
typedef K::Point_2                                       Point2;

// [[Rcpp::export]]
Rcpp::NumericMatrix testparam() {
  const std::string filename = "torus.off";
  
  Mesh3 sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  
  // A halfedge on the border
  halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;
  
  // The 2D points of the uv parametrisation will be written into this map
  typedef Mesh3::Property_map<vertex_descriptor, Point2>  UV_pmap;
  UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point2>("v:uv").first;
  
  typedef SMP::Discrete_conformal_map_parameterizer_3<Mesh3> Parameterizer;
  
  SMP::Error_code err = SMP::parameterize(sm, Parameterizer(), bhd, uv_map);
  
  if(err != SMP::OK) {
    std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
  }
  
  const size_t nvertices = sm.number_of_vertices();
  Rcpp::NumericMatrix out(nvertices, 2);
  int i = 0;
  for(Mesh3::Vertex_index v : sm.vertices()) {
    Point2 pt = uv_map[v];
    out(i, 0) = pt.x();
    out(i, 1) = pt.y();
    i++;
  }
  
  return out;
}