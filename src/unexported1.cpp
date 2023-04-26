#ifndef _HEADER_
#include "cgalMesh.h"
#endif

// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
void Message(std::string msg) {
  SEXP rmsg = Rcpp::wrap(msg);
  Rcpp::message(rmsg);
}

// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
template <typename PointT>
Rcpp::NumericMatrix points3_to_matrix(std::vector<PointT> points) {
  const size_t npoints = points.size();
  Rcpp::NumericMatrix M(3, npoints);
  for(size_t i = 0; i != npoints; i++) {
    Rcpp::NumericVector col_i(3);
    const PointT point = points[i];
    col_i(0) = CGAL::to_double(point.x());
    col_i(1) = CGAL::to_double(point.y());
    col_i(2) = CGAL::to_double(point.z());
    M(Rcpp::_, i) = col_i;
  }
  return M;
}

template Rcpp::NumericMatrix points3_to_matrix<EPoint3>(std::vector<EPoint3>);

// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix M) {
  const size_t npoints = M.ncol();
  std::vector<PointT> points;
  points.reserve(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt = M(Rcpp::_, i);
    points.emplace_back(PointT(pt(0), pt(1), pt(2)));
  }
  return points;
}

template std::vector<EPoint3> matrix_to_points3<EPoint3>(
  const Rcpp::NumericMatrix
);
template std::vector<Point3> matrix_to_points3<Point3>(
  const Rcpp::NumericMatrix
);

// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
std::vector<std::vector<size_t>> list_to_faces(const Rcpp::List L) {
  const size_t nfaces = L.size();
  std::vector<std::vector<size_t>> faces;
  faces.reserve(nfaces);
  for(size_t i = 0; i < nfaces; i++) {
    Rcpp::IntegerVector face_rcpp = Rcpp::as<Rcpp::IntegerVector>(L(i));
    std::vector<size_t> face(face_rcpp.begin(), face_rcpp.end());
    faces.emplace_back(face);
  }
  return faces;
}
