#ifndef _HEADER_
#include "cgalMesh.h"
#endif


Point3 scalePoint(double lambda, Point3 pt) {
  return Point3(lambda * pt.x(), lambda * pt.y(), lambda * pt.z());
}


Point3 addPoints(Point3 pt1, Point3 pt2) {
  return Point3(pt1.x() + pt2.x(), pt1.y() + pt2.y(), pt1.z() + pt2.z());
}


double sqnormPoint(Point3 pt) {
  return pt.x()*pt.x() + pt.y()*pt.y() + pt.z()*pt.z();
}


double betaF(Point3 A, double s) {
  return s / sqrt(s*s + sqnormPoint(A));
}


double gammaF(Point3 A, double s) {
  return s / sqrt(s*s - sqnormPoint(A));
}


Point3 gyromidpoint(Point3 A, Point3 B, double s) {
  Point3 bA = scalePoint(betaF(A, s), A);
  Point3 bB = scalePoint(betaF(B, s), B);
  double gA = gammaF(bA, s);
  double gB = gammaF(bB, s);
  Point3 M = scalePoint(
    1.0 / (gA + gB), 
    addPoints(scalePoint(gA, bA), scalePoint(gB, bB))
  );
  return scalePoint(gammaF(M, s), M);
}


std::pair<vxdescr, vxdescr> orderedPair(vxdescr vi, vxdescr vj) {
  return int(vi) < int(vj) ? std::make_pair(vi, vj) : std::make_pair(vj, vi);
}


Mesh3 gyroQuadrisection(Mesh3 mesh, double s) {
  Mesh3 newmesh;
  for(Mesh3::Vertex_index vd : mesh.vertices()) {
    newmesh.add_vertex(mesh.point(vd));
  }
  int nverts = mesh.number_of_vertices();
  std::map<std::pair<vxdescr, vxdescr>, vxdescr> middles;
  for(Mesh3::Edge_index ed : mesh.edges()) {
    vxdescr v1 = source(ed, mesh);
    vxdescr v2 = target(ed, mesh);
    Point3 gmidpoint = gyromidpoint(mesh.point(v1), mesh.point(v2), s);
    newmesh.add_vertex(gmidpoint);
    middles.insert(std::make_pair(orderedPair(v1, v2), vxdescr(nverts)));
    nverts++;
  }
  for(Mesh3::Face_index fd: mesh.faces()) {
    auto vs = vertices_around_face(mesh.halfedge(fd), mesh).begin();
    vxdescr v1 = *(vs++);
    vxdescr v2 = *(vs++);
    vxdescr v3 = *vs;
    vxdescr m12 = middles[orderedPair(v1, v2)];
    vxdescr m23 = middles[orderedPair(v2, v3)];
    vxdescr m31 = middles[orderedPair(v3, v1)];
    newmesh.add_face(v1, m12, m31);
    newmesh.add_face(v2, m23, m12);
    newmesh.add_face(v3, m31, m23);
    newmesh.add_face(m12, m23, m31);
  }
  return newmesh;
}


// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> gTriangle(
  Rcpp::NumericVector A, Rcpp::NumericVector B, Rcpp::NumericVector C,
  double s, int iterations
) {
  Point3 pa(A(0), A(1), A(2));
  Point3 pb(B(0), B(1), B(2));
  Point3 pc(C(0), C(1), C(2));
  Mesh3 mesh;
  mesh.add_vertex(pa);
  mesh.add_vertex(pb);
  mesh.add_vertex(pc);
  mesh.add_face(vxdescr(0), vxdescr(1), vxdescr(2));
  std::vector<Mesh3> meshes(iterations);
  meshes[0] = mesh;
  for(int i = 1; i < iterations; i++) {
    meshes[i] = gyroQuadrisection(meshes[i-1], s);
  }
  EMesh3 emesh;
  CGAL::copy_face_graph(meshes[iterations-1], emesh);
  return Rcpp::XPtr<EMesh3>(new EMesh3(emesh), false);
}
