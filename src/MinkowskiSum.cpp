#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/minkowski_sum_3.h>
typedef CGAL::Nef_polyhedron_3<EK> ENef3;

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> MinkowskiSum_cpp(
  Rcpp::XPtr<EMesh3> mesh1XPtr, Rcpp::XPtr<EMesh3> mesh2XPtr
) {
  EMesh3 mesh1 = *(mesh1XPtr.get());
  EMesh3 mesh2 = *(mesh2XPtr.get());
  if(!CGAL::is_triangle_mesh(mesh1)) {
    Rcpp::stop("The first mesh is not triangle.");
  }
  if(!CGAL::is_triangle_mesh(mesh2)) {
    Rcpp::stop("The second mesh is not triangle.");
  }
  ENef3 nef1(mesh1);
  ENef3 nef2(mesh2);
  ENef3 nef = CGAL::minkowski_sum_3(nef1, nef2);
  EMesh3 mesh;
  CGAL::convert_nef_polyhedron_to_polygon_mesh(nef, mesh, false);
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);  
}
