#ifndef _HEADER_
#include "cgalMesh.h"
#endif

EMesh3 icosphere(EPoint3 center, EK::FT radius, unsigned int iterations) {
  EMesh3 mesh;
  CGAL::make_icosahedron<EMesh3, EPoint3>(mesh, center, radius);
  CGAL::Subdivision_method_3::Loop_subdivision(
    mesh, CGAL::parameters::number_of_iterations(iterations)
  );
  return mesh;
}

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> sTriangle(
  Rcpp::NumericVector A, Rcpp::NumericVector B, Rcpp::NumericVector C,
  Rcpp::NumericVector center, double radius, unsigned int iterations
) {
  EPoint3 pa(A(0), A(1), A(2));
  EPoint3 pb(B(0), B(1), B(2));
  EPoint3 pc(C(0), C(1), C(2));
  EPoint3 O(center(0), center(1), center(2));
  EPlane3 Oab(O, pa, pb);
  if(Oab.has_on_negative_side(pc)) {
    Oab = Oab.opposite();
  }
  EPlane3 Obc(O, pb, pc);
  if(Obc.has_on_negative_side(pa)) {
    Obc = Obc.opposite();
  }
  EPlane3 Oca(O, pc, pa);
  if(Oca.has_on_negative_side(pb)) {
    Oca = Oca.opposite();
  }

  EK::FT r(radius);
  EMesh3 sphereMesh = icosphere(O, r, iterations);

  clippingToPlane(sphereMesh, Oab, false);
  clippingToPlane(sphereMesh, Obc, false);
  clippingToPlane(sphereMesh, Oca, false);

  Normals_map vnormal = 
    sphereMesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
      "v:normal", defaultNormal()
    ).first;
  for(vertex_descriptor vd : sphereMesh.vertices()) {
    EPoint3 vx = sphereMesh.point(vd);
    Rcpp::NumericVector normal = {
      CGAL::to_double<EK::FT>(vx.x()),
      CGAL::to_double<EK::FT>(vx.y()),
      CGAL::to_double<EK::FT>(vx.z())
    };
    vnormal[vd] = normal;
  }

  return Rcpp::XPtr<EMesh3>(new EMesh3(sphereMesh), false);
} 
