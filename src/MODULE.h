#ifndef _HEADER_
#include "cgalMesh.h"
#endif

class CGALmesh {
public:
  EMesh3 mesh;
  Rcpp::XPtr<EMesh3> xptr;
  CGALmesh(const Rcpp::NumericMatrix vertices,
           const Rcpp::List faces,
           const bool clean)
    : mesh(
        csoup2mesh<EMesh3, EPoint3>(
          matrix_to_points3<EPoint3>(vertices), 
          list_to_faces(faces), 
          clean
      )
    ), 
      xptr(Rcpp::XPtr<EMesh3>(&mesh, false)) {}
  CGALmesh(Rcpp::XPtr<EMesh3> xptr_)
    : mesh(*(xptr_.get())), xptr(Rcpp::XPtr<EMesh3>(&mesh, false)) {}
  
  Rcpp::NumericVector centroid() {
    const EPoint3 centroid = PMP::centroid(mesh);
    Rcpp::NumericVector out(3);
    out(0) = CGAL::to_double<EK::FT>(centroid.x());
    out(1) = CGAL::to_double<EK::FT>(centroid.y());
    out(2) = CGAL::to_double<EK::FT>(centroid.z());
    return out;
  }
  
  Rcpp::XPtr<EMesh3> clone() {
    EMesh3 copy;
    CGAL::copy_face_graph(mesh, copy);
    return Rcpp::XPtr<EMesh3>(new EMesh3(copy), false);
  }
  
  bool doesSelfIntersect() {
    return PMP::does_self_intersect(mesh);
  }
  
  Rcpp::List getRmesh(const bool normals) {
    return RSurfEKMesh(mesh, normals);
  }
  
  bool isTriangle() {
    return CGAL::is_triangle_mesh(mesh);
  }
  
  void print() {
    Rcpp::Rcout << "Mesh with " << mesh.number_of_vertices() 
                << " vertices and " << mesh.number_of_faces() << " faces.\n";
  }
  
  void triangulate() {
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  
};