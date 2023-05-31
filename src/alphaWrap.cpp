#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/alpha_wrap_3.h>

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> alphaWrap_cpp(
    const Rcpp::NumericMatrix pts, const double ralpha, const double roffset
) {
  
  std::vector<Point3> points = matrix_to_points3<Point3>(pts);
  
  CGAL::Bbox_3 bbox = CGAL::bbox_3(std::cbegin(points), std::cend(points));
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double alpha = diag_length / ralpha;
  const double offset = diag_length / roffset;
  
  Mesh3 wrap;
  CGAL::alpha_wrap_3(points, alpha, offset, wrap);
  EMesh3 out;
  CGAL::copy_face_graph(wrap, out);
  
  return Rcpp::XPtr<EMesh3>(new EMesh3(out), false);
}