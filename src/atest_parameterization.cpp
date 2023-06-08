#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>

namespace SMP = CGAL::Surface_mesh_parameterization;
typedef K::Point_2                                       Point2;

// [[Rcpp::export]]
Rcpp::NumericMatrix testparam(std::string filename, int method) {
  
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
  
  SMP::Error_code err;

  typedef SMP::Discrete_conformal_map_parameterizer_3<Mesh3> Parameterizer1;
  
  typedef SMP::Circular_border_arc_length_parameterizer_3<Mesh3>  Border_parameterizer;
  typedef SMP::Discrete_authalic_parameterizer_3<Mesh3, Border_parameterizer> Parameterizer2;
  
  typedef SMP::ARAP_parameterizer_3<Mesh3, Border_parameterizer> Parameterizer3;
  
  if(method == 1) {
    err = SMP::parameterize(sm, Parameterizer1(), bhd, uv_map);
  } else if(method == 2) {
    err = SMP::parameterize(sm, Parameterizer2(), bhd, uv_map);
  } else {
    err = SMP::parameterize(sm, Parameterizer3(10), bhd, uv_map);
  }
  
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