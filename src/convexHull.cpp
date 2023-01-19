#ifndef _HEADER_
#include "cgalMesh.h"
#endif

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> cxhull(Rcpp::NumericMatrix pts) {
  // make points
  int npoints = pts.ncol();
  std::vector<EPoint3> points;
  points.reserve(npoints);
  for(int i = 0; i < npoints; i++) {
    Rcpp::NumericVector pt = pts(Rcpp::_, i);
    points.emplace_back(EPoint3(pt(0), pt(1), pt(2)));
  }
  // convex hull
  EMesh3 mesh;
  CGAL::convex_hull_3(points.begin(), points.end(), mesh);
  //
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}