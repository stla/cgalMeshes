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
  CGALmesh(const std::string filename, const bool binary)
    : mesh(readMeshFile(filename)), 
      xptr(Rcpp::XPtr<EMesh3>(&mesh, false)) {}

  double area() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The mesh self-intersects.");
    }
    const EK::FT ar = PMP::area(mesh);
    return CGAL::to_double<EK::FT>(ar);
  }
  
  Rcpp::NumericVector centroid() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    const EPoint3 centroid = PMP::centroid(mesh);
    Rcpp::NumericVector out(3);
    out(0) = CGAL::to_double<EK::FT>(centroid.x());
    out(1) = CGAL::to_double<EK::FT>(centroid.y());
    out(2) = CGAL::to_double<EK::FT>(centroid.z());
    return out;
  }
  
  void clipMesh(Rcpp::XPtr<EMesh3> clipperXPtr, const bool clipVolume) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(clipVolume) {
      if(PMP::does_self_intersect(mesh)) {
        Rcpp::stop("The mesh self-intersects.");
      }
    }
    EMesh3 clipper = *(clipperXPtr.get());
    if(!CGAL::is_triangle_mesh(clipper)) {
      Rcpp::stop("The clipping mesh is not triangle.");
    }
    if(!CGAL::is_closed(mesh)) {
      Rcpp::stop("The clipping mesh is not closed.");
    }
    if(PMP::does_self_intersect(clipper)) {
      Rcpp::stop("The clipping mesh self-intersects.");
    }
    const bool doNotModify = !clipVolume;
    const bool clipping = PMP::clip(
      mesh, clipper, PMP::parameters::clip_volume(clipVolume),
      PMP::parameters::clip_volume(clipVolume).do_not_modify(doNotModify)
    );
    if(!clipping) {
      Rcpp::stop("Clipping has failed.");
    }
    mesh.collect_garbage();
  }
  
  Rcpp::XPtr<EMesh3> clone() {
    EMesh3 copy;
    CGAL::copy_face_graph(mesh, copy);
    return Rcpp::XPtr<EMesh3>(new EMesh3(copy), false);
  }

  bool doesBoundVolume() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    return PMP::does_bound_a_volume(mesh);
  }
  
  bool doesSelfIntersect() {
    return PMP::does_self_intersect(mesh);
  }
  
  Rcpp::DataFrame edges() {
    return getEdges<EK, EMesh3, EPoint3>(mesh);
  }
  
  Rcpp::List getRmesh(const bool normals) {
    if(CGAL::is_triangle_mesh(mesh)) {
      return RSurfEKMesh2(mesh, normals, 3);
    }
    if(CGAL::is_quad_mesh(mesh)) {
      return RSurfEKMesh2(mesh, normals, 4);
    }
    return RSurfEKMesh(mesh, normals);
  }

  bool isClosed() {
    return CGAL::is_closed(mesh);
  }

  bool isOutwardOriented() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    return PMP::is_outward_oriented(mesh);
  }
  
  bool isTriangle() {
    return CGAL::is_triangle_mesh(mesh);
  }

  void orientToBoundVolume() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    PMP::orient_to_bound_a_volume(mesh);
  }
  
  void print() {
    Rcpp::Rcout << "Mesh with " << mesh.number_of_vertices() 
                << " vertices and " << mesh.number_of_faces() << " faces.\n";
  }
  
  void reverseFaceOrientations() {
    PMP::reverse_face_orientations(mesh);  
  }
  
  void triangulate() {
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
  }
  
  Rcpp::NumericMatrix vertices() {
    return getVertices_EK(mesh);
  }
  
  double volume() {
    if(!CGAL::is_closed(mesh)) {
      Rcpp::stop("The mesh is not closed.");
    }
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The mesh self-intersects.");
    }
    const EK::FT vol = PMP::volume(mesh);
    return CGAL::to_double<EK::FT>(vol);
  }
  
  
  void writeFile(
      Rcpp::String filename, const int precision, const bool binary
  ) {
    writeMeshFile(filename, precision, binary, mesh);
  }
  
};