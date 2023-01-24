#ifndef _HEADER_
#include "cgalMesh.h"
#endif

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> Mandelbulb(
  int maxloop,
  Rcpp::NumericVector sphereCenter, double sphereRadius,
  double angle_bound, double radius_bound, double distance_bound, 
  double error_bound
) {
  Tri tr;             // 3D-Delaunay triangulation
  Cplx2 cplx2(tr);    // 2D-complex in 3D-Delaunay triangulation

  // isosurface fun(p)=0
  auto fun = [maxloop](Point_3 p) {
    FT x0 = p.x();
    FT y0 = p.y();
    FT z0 = p.z();
    FT x = x0;
    FT y = y0;
    FT z = z0;
    FT r2, theta, phi, r8;
    for(int i = 0; i < maxloop; i++){
      r2 = x*x + y*y + z*z;
      if(r2 > 4){
        break;
      }
      theta = 8 * atan2(sqrt(x*x + y*y), z);
      phi = 8 * atan2(y, x);
      r8 = r2 * r2 * r2 * r2;
      x = r8 * cos(phi) * sin(theta) + x0;
      y = r8 * sin(phi) * sin(theta) + y0;
      z = r8 * cos(theta) + z0;
    }
    return sqrt(r2) - 2;
  };

  // bounding sphere
  Point_3 bounding_sphere_center(
    sphereCenter(0), sphereCenter(1), sphereCenter(2)
  );
  FT bounding_sphere_squared_radius = sphereRadius * sphereRadius;
  Sphere_3 bounding_sphere(
    bounding_sphere_center, bounding_sphere_squared_radius
  );

  // isosurface
  FT eb(error_bound);
  ImplicitSurface surface(fun, bounding_sphere, eb);

  // defining meshing criteria
  FT ab(angle_bound);
  FT rb(radius_bound);
  FT db(distance_bound); 
  MeshingCriteria criteria(ab, rb, db);

  // meshing surface
  CGAL::make_surface_mesh(
    cplx2, surface, criteria, CGAL::Manifold_with_boundary_tag()
  );
  SurfaceMesh smesh;
  CGAL::facets_in_complex_2_to_triangle_mesh(cplx2, smesh);
  PMP::orient_to_bound_a_volume(smesh);
  if(!PMP::is_outward_oriented(smesh)) {
    PMP::reverse_face_orientations(smesh);
  }

  // convert to EMesh3
  EMesh3 mesh;
  CGAL::copy_face_graph(smesh, mesh);

  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}

