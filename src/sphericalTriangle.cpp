#ifndef _HEADER_
#include "cgalMesh.h"
#endif

void clippingToPlane(EMesh3& mesh, EPlane3 plane) {
  const bool clipping = PMP::clip(
    mesh, plane,
    PMP::parameters::clip_volume(false)
  );
  if(!clipping) {
    Rcpp::stop("Clipping has failed.");
  }
  mesh.collect_garbage();
}


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
  double cx = center(0);
  double cy = center(1);
  double cz = center(2);
  EPoint3 O(cx, cy, cz);
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

  clippingToPlane(sphereMesh, Oab);
  clippingToPlane(sphereMesh, Obc);
  clippingToPlane(sphereMesh, Oca);

  Normals_map vnormal = 
    sphereMesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
      "v:normal", defaultNormal()
    ).first;
  for(vertex_descriptor vd : sphereMesh.vertices()) {
    EPoint3 vx = sphereMesh.point(vd);
    Rcpp::NumericVector normal = {
      CGAL::to_double<EK::FT>(vx.x()) - cx,
      CGAL::to_double<EK::FT>(vx.y()) - cy,
      CGAL::to_double<EK::FT>(vx.z()) - cz
    };
    vnormal[vd] = normal;
  }

  return Rcpp::XPtr<EMesh3>(new EMesh3(sphereMesh), false);
} 
