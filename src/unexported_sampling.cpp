#ifndef _HEADER_
#include "cgalMesh.h"
#endif


// ------------------------------------------------------------------------- //
// - volume of tetrahedron p1p2p3p4 ---------------------------------------- //
double volumeTetrahedron(Point3 p1, Point3 p2, Point3 p3, Point3 p4) {
  Vector3 v = p1 - p4;
  Vector3 w = CGAL::cross_product(p2 - p4, p3 - p4);
  double dotprod = CGAL::scalar_product(v, w);
  return fabs(dotprod) / 6.0;
}


// ------------------------------------------------------------------------- //
// - point to vector ------------------------------------------------------- //
Vector3 P3toV3(Point3 p) {
  Vector3 v(p.x(), p.y(), p.z());
  return v;
}


// ------------------------------------------------------------------------- //
// - vector to point ------------------------------------------------------- //
Point3 V3toP3(Vector3 v) {
  Point3 p(v.x(), v.y(), v.z());
  return p;
}


// ------------------------------------------------------------------------- //
// - five tetrahedra partitioning a hexahedron ----------------------------- //
std::array<std::array<Vector3, 4>, 5> hexahedronTetrahedra(
    std::array<Point3, 8> hxh
) {
  Vector3 v1 = P3toV3(hxh[0]);
  Vector3 v2 = P3toV3(hxh[1]);
  Vector3 v3 = P3toV3(hxh[2]);
  Vector3 v4 = P3toV3(hxh[3]);
  Vector3 v5 = P3toV3(hxh[4]);
  Vector3 v6 = P3toV3(hxh[5]);
  Vector3 v7 = P3toV3(hxh[6]);
  Vector3 v8 = P3toV3(hxh[7]);
  std::array<Vector3, 4> th1 = {v1, v5, v3, v7};
  std::array<Vector3, 4> th2 = {v7, v1, v2, v3};
  std::array<Vector3, 4> th3 = {v6, v1, v7, v5};
  std::array<Vector3, 4> th4 = {v8, v7, v3, v5};
  std::array<Vector3, 4> th5 = {v4, v1, v3, v5};
  std::array<std::array<Vector3, 4>, 5> tetrahedra = {th1, th2, th3, th4, th5};
  return tetrahedra;
}


// ------------------------------------------------------------------------- //
// - sample one point in tetrahedron v1v2v3v4  ----------------------------- //
Vector3 sampleTetrahedron(
  Vector3 v1, Vector3 v2, Vector3 v3, Vector3 v4, boost::mt19937 gen
) {
  boost::random::uniform_real_distribution<double> runif(0.0, 1.0);
  double c1 = runif(gen);
  double c2 = runif(gen);
  double c3 = runif(gen);
  double t;
  if(c1 + c2 > 1) {
    c1 = 1 - c1;
    c2 = 1 - c2;
  }
  if(c2 + c3 > 1) {
    t = c3;
    c3 = 1 - c1 - c2;
    c2 = 1 - t;
  } else if(c1 + c2 + c3 > 1) {
    t = c3;
    c3 = c1 + c2 + c3 - 1;
    c1 = 1 - c2 - t;
  }
  double c4 = 1 - (c1 + c2 + c3);
  return c1 * v1 + c2 * v2 + c3 * v3 + c4 * v4;
}
