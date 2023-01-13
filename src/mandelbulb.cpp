#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tr;

// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tr> C2t3;

typedef Tr::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;

typedef FT (*Function)(Point_3);

typedef CGAL::Implicit_surface_3<GT, Function> Surface_3;

typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

////

FT mandelbulb_function (Point_3 p) {
  FT x0 = p.x();
  FT y0 = p.y();
  FT z0 = p.z();
  FT x = x0;
  FT y = y0;
  FT z = z0;
  FT r2, theta, phi, r8;
  for(int i = 0; i < 24; i++){
    r2 = x*x + y*y + z*z;
    if(r2 > 4){
      break;
    }
    theta = 8*atan2(sqrt(x*x+y*y),z);
    phi = 8*atan2(y,x);
    r8 = r2*r2*r2*r2;
    x = r8*cos(phi)*sin(theta) + x0;
    y = r8*sin(phi)*sin(theta) + y0;
    z = r8*cos(theta) + z0;
  }
  return sqrt(r2) - 2;
}

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> mandelbulb(double angle_bound, double radius_bound, double distance_bound) {
  Tr tr;            // 3D-Delaunay triangulation
  C2t3 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  Surface_3 surface(mandelbulb_function,             // pointer to function
                    Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
  // Note that "2." above is the *squared* radius of the bounding sphere!

  // defining meshing criteria
  CGAL::Surface_mesh_default_criteria_3<Tr> criteria(35.,  // angular bound
                                                     0.025,  // radius bound
                                                     0.025); // distance bound
  // typedef CGAL::Surface_mesh_default_criteria_3<Tr> Crit; 
  // Crit criteria(tr::FT(angle_bound),  // angular bound
  //                                                    tr::FT(radius_bound),  // radius bound
  //                                                    tr::FT(distance_bound)); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  Surface_mesh sm;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

  EMesh3 mesh;
  CGAL::copy_face_graph(sm, mesh);
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}
