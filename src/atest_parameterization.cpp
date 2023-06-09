#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Discrete_authalic_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>

#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>
#include <CGAL/Unique_hash_map.h>

#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>

namespace SMP = CGAL::Surface_mesh_parameterization;
typedef K::Point_2                                       Point2;

typedef CGAL::Unique_hash_map<vertex_descriptor, Point2>        UV_uhm;
typedef boost::associative_property_map<UV_uhm>                 UV_pmap2;

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
  UV_pmap uv_map;

  typedef Mesh3::Property_map<vertex_descriptor, int>  VI_pmap;
  typedef Mesh3::Property_map<vertex_descriptor, bool>  VB_pmap;
  VI_pmap vi_map = sm.add_property_map<vertex_descriptor, int>("v:vi").first;
  VB_pmap vb_map = sm.add_property_map<vertex_descriptor, bool>("v:vp", false).first;
  
  int k = 0;
  for(vertex_descriptor v : sm.vertices()) {
    vi_map[v] = k++;
  }
  
  UV_uhm uv_uhm;
  UV_pmap2 uv_map2(uv_uhm);
  
  SMP::Error_code err;

  typedef SMP::Discrete_conformal_map_parameterizer_3<Mesh3> Parameterizer1;
  
  typedef SMP::Circular_border_arc_length_parameterizer_3<Mesh3>  Border_parameterizer;
  typedef SMP::Discrete_authalic_parameterizer_3<Mesh3, Border_parameterizer> Parameterizer2;
  
  typedef SMP::ARAP_parameterizer_3<Mesh3, Border_parameterizer> Parameterizer3;
  
  typedef SMP::Iterative_authalic_parameterizer_3<Mesh3, Border_parameterizer> Parameterizer4;

  typedef SMP::Square_border_uniform_parameterizer_3<Mesh3> SquareBorder_parameterizer;
  
  typedef SMP::Discrete_conformal_map_parameterizer_3<Mesh3, SquareBorder_parameterizer> Parameterizer5;
  
  if(method == 1) {
    uv_map = sm.add_property_map<vertex_descriptor, Point2>("v:uv").first;
    err = SMP::parameterize(sm, Parameterizer1(), bhd, uv_map);
  } else if(method == 2) {
    uv_map = sm.add_property_map<vertex_descriptor, Point2>("v:uv").first;
    err = SMP::parameterize(sm, Parameterizer2(), bhd, uv_map);
  } else if(method == 3){
    uv_map = sm.add_property_map<vertex_descriptor, Point2>("v:uv").first;
    // Border_parameterizer border_parameterizer;
    // Parameterizer3 parameterizer;
    // err = parameterizer.parameterize(sm, bhd, uv_map, vi_map, vb_map);
    SMP::ARAP_parameterizer_3<SurfaceMesh> parameterizer;
    err = SMP::parameterize(sm, parameterizer, bhd, uv_map);
  } else if(method == 4) {
    Border_parameterizer border_parameterizer; // the border parameterizer will automatically compute the corner vertices
    Parameterizer4 parameterizer(border_parameterizer);
    const unsigned int iterations = 2;
    err = parameterizer.parameterize(sm, bhd, uv_map2, iterations);
  } else {
    SquareBorder_parameterizer border_param; // the border parameterizer will compute the corner vertices
    err = SMP::parameterize(sm, Parameterizer5(border_param), bhd, uv_map2);    
  }
  
  if(err != SMP::OK) {
    std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
    return EXIT_FAILURE;
  }
  
  const size_t nvertices = sm.number_of_vertices();
  Rcpp::NumericMatrix out(nvertices, 2);
  int i = 0;
  for(Mesh3::Vertex_index v : sm.vertices()) {
    Point2 pt = method < 4 ? uv_map[v] : uv_map2[v];
    out(i, 0) = pt.x();
    out(i, 1) = pt.y();
    i++;
  }
  
  return out;
}