#ifndef _HEADER_
#include "cgalMesh.h"
#endif


Rcpp::NumericMatrix getVertices_EK(EMesh3& mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Vertices(3, nvertices);
  {
    size_t i = 0;
    for(EMesh3::Vertex_index vd : mesh.vertices()) {
      Rcpp::NumericVector col_i(3);
      const EPoint3 vertex = mesh.point(vd);
      col_i(0) = CGAL::to_double<EK::FT>(vertex.x());
      col_i(1) = CGAL::to_double<EK::FT>(vertex.y());
      col_i(2) = CGAL::to_double<EK::FT>(vertex.z());
      Vertices(Rcpp::_, i++) = col_i;
    }
  }
  return Vertices;
}


template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame getEdges(MeshT& mesh) {
  const size_t nedges = mesh.number_of_edges();
  Rcpp::IntegerVector I1(nedges);
  Rcpp::IntegerVector I2(nedges);
  Rcpp::NumericVector Length(nedges);
  Rcpp::LogicalVector Border(nedges);
  Rcpp::NumericVector Angle(nedges);
  {
    size_t i = 0;
    for(typename MeshT::Edge_index ed : mesh.edges()) {
      typename MeshT::Vertex_index s = source(ed, mesh);
      typename MeshT::Vertex_index t = target(ed, mesh);
      I1(i) = (int)s + 1;
      I2(i) = (int)t + 1;
      typename MeshT::Halfedge_index h0 = mesh.halfedge(ed, 0);
      typename KernelT::FT el = PMP::edge_length(h0, mesh);
      Length(i) = CGAL::to_double(el);
      const bool isBorder = mesh.is_border(ed);
      Border(i) = isBorder;
      if(isBorder) {
        Angle(i) = Rcpp::NumericVector::get_na();
      } else {
        std::vector<PointT> points(4);
        points[0] = mesh.point(s);
        points[1] = mesh.point(t);
        points[2] = mesh.point(mesh.target(mesh.next(h0)));
        typename MeshT::Halfedge_index h1 = mesh.halfedge(ed, 1);
        points[3] = mesh.point(mesh.target(mesh.next(h1)));
        typename KernelT::FT angle = CGAL::abs(CGAL::approximate_dihedral_angle(
            points[0], points[1], points[2], points[3]));
        Angle(i) = CGAL::to_double(angle);
      }
      i++;
    }
  }
  Rcpp::DataFrame Edges = Rcpp::DataFrame::create(
    Rcpp::Named("i1")       = I1,
    Rcpp::Named("i2")       = I2,
    Rcpp::Named("length")   = Length,
    Rcpp::Named("border")   = Border,
    Rcpp::Named("angle")    = Angle
  );
  return Edges;
}

template Rcpp::DataFrame getEdges<EK, EMesh3, EPoint3>(EMesh3&);


template <typename MeshT>
Rcpp::List getFaces(MeshT& mesh) {
  const size_t nfaces = mesh.number_of_faces();
  Rcpp::List Faces(nfaces);
  {
    size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      Faces(i++) = col_i;
    }
  }
  return Faces;
}

template Rcpp::List getFaces<EMesh3>(EMesh3&);


template <typename MeshT>
Rcpp::IntegerMatrix getFaces2(MeshT& mesh, const int nsides) {
  const size_t nfaces = mesh.number_of_faces();
  Rcpp::IntegerMatrix Faces(nsides, nfaces);
  {
    size_t i = 0;
    for(typename MeshT::Face_index fd : mesh.faces()) {
      // bool TEST = CGAL::is_triangle(mesh.halfedge(fd), mesh);
      Rcpp::IntegerVector col_i;
      for(typename MeshT::Vertex_index vd :
          vertices_around_face(mesh.halfedge(fd), mesh)) {
        col_i.push_back(vd + 1);
      }
      Faces(Rcpp::_, i++) = col_i;
    }
  }
  return Faces;
}


Rcpp::NumericMatrix getEKNormals(EMesh3& mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Normals(3, nvertices);
  std::pair<CGALnormals_map, bool> vnormals_ = 
    mesh.property_map<vertex_descriptor, EVector3>("v:normals");
  if(vnormals_.second) {
    mesh.remove_property_map(vnormals_.first);
  }
  CGALnormals_map vnormals = 
    mesh.add_property_map<EMesh3::Vertex_index, EVector3>(
                          "v:normals", CGAL::NULL_VECTOR
                        ).first;
  // auto fnormals = mesh.add_property_map<EMesh3::Face_index, EVector3>(
  //                         "f:normals", CGAL::NULL_VECTOR)
  //                     .first;
  PMP::compute_vertex_normals(mesh, vnormals);
  {
    size_t i = 0;
    for(EMesh3::Vertex_index vd : vertices(mesh)) {
      Rcpp::NumericVector col_i(3);
      const EVector3 normal = vnormals[vd];
      col_i(0) = CGAL::to_double<EK::FT>(normal.x());
      col_i(1) = CGAL::to_double<EK::FT>(normal.y());
      col_i(2) = CGAL::to_double<EK::FT>(normal.z());
      Normals(Rcpp::_, i++) = col_i;
    }
  }
  return Normals;
}


Rcpp::List RSurfEKMesh(EMesh3& mesh, const bool normals) {
  Rcpp::NumericMatrix Vertices = getVertices_EK(mesh);
  Rcpp::List Faces = getFaces<EMesh3>(mesh);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getEKNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}


Rcpp::List RSurfEKMesh2(EMesh3& mesh, const bool normals, const int nsides) {
  Rcpp::NumericMatrix Vertices = getVertices_EK(mesh);
  Rcpp::IntegerMatrix Faces = getFaces2<EMesh3>(mesh, nsides);
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("vertices") = Vertices,
                                      Rcpp::Named("faces") = Faces);
  if(normals) {
    Rcpp::NumericMatrix Normals = getEKNormals(mesh);
    out["normals"] = Normals;
  }
  return out;
}
