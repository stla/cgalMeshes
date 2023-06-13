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


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
Mesh3 epeck2epick(EMesh3& emesh) {
  const size_t nvertices = emesh.number_of_vertices();
  const size_t nedges    = emesh.number_of_edges();
  const size_t nfaces    = emesh.number_of_faces();
  Mesh3 mesh;
  mesh.reserve(nvertices, nedges, nfaces);
  for(EMesh3::Vertex_index vd : emesh.vertices()) {
    const EPoint3 vertex = emesh.point(vd);
    const double x = CGAL::to_double<EK::FT>(vertex.x());
    const double y = CGAL::to_double<EK::FT>(vertex.y());
    const double z = CGAL::to_double<EK::FT>(vertex.z());
    mesh.add_vertex(Point3(x, y, z));
  }
  for(EMesh3::Face_index fd : emesh.faces()) {
    std::vector<int> face;
    for(EMesh3::Vertex_index vd : 
          vertices_around_face(emesh.halfedge(fd), emesh)) {
      face.push_back(vd);
    }
    mesh.add_face(
      CGAL::SM_Vertex_index(face[0]), 
      CGAL::SM_Vertex_index(face[1]), 
      CGAL::SM_Vertex_index(face[2])
    );
  }
  return mesh;
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
EMesh3 epick2epeck(Mesh3& mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  const size_t nedges    = mesh.number_of_edges();
  const size_t nfaces    = mesh.number_of_faces();
  EMesh3 emesh;
  emesh.reserve(nvertices, nedges, nfaces);
  for(Mesh3::Vertex_index vd : mesh.vertices()) {
    const Point3 vertex = mesh.point(vd);
    const double x = vertex.x();
    const double y = vertex.y();
    const double z = vertex.z();
    emesh.add_vertex(EPoint3(x, y, z));
  }
  for(Mesh3::Face_index fd : mesh.faces()) {
    std::vector<int> face;
    for(Mesh3::Vertex_index vd : 
          vertices_around_face(mesh.halfedge(fd), mesh)) {
      face.push_back(vd);
    }
    emesh.add_face(
      CGAL::SM_Vertex_index(face[0]), 
      CGAL::SM_Vertex_index(face[1]), 
      CGAL::SM_Vertex_index(face[2])
    );
  }
  return emesh;
}
