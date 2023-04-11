#ifndef _HEADER_
#include "cgalMesh.h"
#endif

EMesh3 cxhullMesh(Rcpp::NumericMatrix pts) {
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
  return mesh;
}


// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> cxhull(Rcpp::NumericMatrix pts) {
  EMesh3 mesh = cxhullMesh(pts);
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}


// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> cxhullsIntersection(
  Rcpp::List Pts, Rcpp::Nullable<Rcpp::NumericVector> origin_
) {
  // make convex hulls
  int nmeshes = Pts.size();
  std::vector<EMesh3> meshes;
  meshes.reserve(nmeshes);
  int nfaces = 0;
  for(int i = 0; i < nmeshes; i++) {
    Rcpp::NumericMatrix pts = Rcpp::as<Rcpp::NumericMatrix>(Pts(i));
    EMesh3 cxmesh = cxhullMesh(pts);
    meshes.emplace_back(cxmesh);
    int nf = cxmesh.number_of_faces();
    if(nf < 4) {
      Rcpp::stop("Found a flat convex hull.");
    }
    nfaces += nf;
  }
  // make planes
  std::vector<EPlane3> planes;
  planes.reserve(nfaces);
  for(int i = 0; i < nmeshes; i++) {
    EMesh3 cxmesh = meshes[i];
    for(face_descriptor fd : cxmesh.faces()) {
      auto vs = vertices_around_face(cxmesh.halfedge(fd), cxmesh).begin();
      vertex_descriptor v3 = *(vs++);
      vertex_descriptor v2 = *(vs++);
      vertex_descriptor v1 = *vs;
      EPlane3 plane(cxmesh.point(v1), cxmesh.point(v2), cxmesh.point(v3));
      planes.emplace_back(plane);
    }
  }
  // intersection
  EMesh3 mesh;
  if(origin_.isNotNull()) {
    Rcpp::NumericVector origin(origin_);
    EPoint3 origin_pt(origin(0), origin(1), origin(2));
    CGAL::halfspace_intersection_3(
      planes.begin(), planes.end(), mesh, origin_pt
    );
  } else {
    CGAL::halfspace_intersection_3(
      planes.begin(), planes.end(), mesh
    );
  }
  //
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);  
}
