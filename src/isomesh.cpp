#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>

#include <CGAL/Surface_mesh_default_criteria_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Gray_level_image_3.h>

#include <CGAL/Polynomial.h>
#include <CGAL/Polynomial_traits_d.h>
#include <CGAL/Polynomial_type_generator.h>
#include <CGAL/polynomial_utils.h>

// default triangulation for Surface_mesher
typedef CGAL::Surface_mesh_default_triangulation_3 Tri;
typedef CGAL::Surface_mesh_default_criteria_3<Tri> MeshingCriteria; 
// c2t3
typedef CGAL::Complex_2_in_triangulation_3<Tri> Cplx2;
typedef Tri::Geom_traits GT;
typedef GT::Sphere_3 Sphere_3;
typedef GT::Point_3 Point_3;
typedef GT::FT FT;
typedef FT (*Function)(Point_3);
typedef CGAL::Implicit_surface_3<GT, Function> ImplicitSurface;
typedef CGAL::Surface_mesh<Point_3> SurfaceMesh;

typedef CGAL::Gray_level_image_3<FT, Point_3> Gray_level_image;
typedef CGAL::Implicit_surface_3<GT, Gray_level_image> Surface_gray;

typedef CGAL::Polynomial_type_generator<FT, 3>::Type     Poly3;
typedef CGAL::Polynomial_traits_d<Poly3>                 PT3;
typedef PT3::Innermost_coefficient_type                  Real;

////

Poly3 Polynomial(Rcpp::IntegerMatrix powers, Rcpp::NumericVector coeffs) {

  PT3::Construct_polynomial construct_polynomial;

  std::list<std::pair<CGAL::Exponent_vector, Real>> innermost_coeffs;
  for(int i = 0; i < coeffs.size(); i++) {
    Rcpp::IntegerVector pows = powers(i, Rcpp::_);
    FT coeff(coeffs(i));
    innermost_coeffs.push_back(
      std::make_pair(CGAL::Exponent_vector(pows(0), pows(1), pows(2)), coeff)
    );
  }
  return 
    construct_polynomial(innermost_coeffs.begin(), innermost_coeffs.end());
}


// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> algebraicMesh(
  Rcpp::IntegerMatrix powers, Rcpp::NumericVector coeffs, double isolevel,
  Rcpp::NumericVector sphereCenter, double sphereRadius,
  double angle_bound, double radius_bound, double distance_bound
) {
  Tri tr;            // 3D-Delaunay triangulation
  Cplx2 cplx2(tr);    // 2D-complex in 3D-Delaunay triangulation

  // isosurface fun(p)=0
  Poly3 P = Polynomial(powers, coeffs);
  PT3::Substitute substitute;
  FT isoval(isolevel);

  auto fun = [P, substitute, isoval](Point_3 p) {
    FT x = p.x();
    FT y = p.y();
    FT z = p.z();
    std::list<Real> xyz;
    xyz.push_back(x); 
    xyz.push_back(y);
    xyz.push_back(z);
    FT val = substitute(P, xyz.begin(), xyz.end());
    return val - isoval;
  };

  // bounding sphere
  Point_3 bounding_sphere_center(
    sphereCenter(0), sphereCenter(1), sphereCenter(2)
  );
  FT bounding_sphere_squared_radius = sphereRadius * sphereRadius;
  Sphere_3 bounding_sphere(
    bounding_sphere_center, bounding_sphere_squared_radius
  );

  // check
  {
    FT val = fun(bounding_sphere_center);
    if(val >= 0) {
      Rcpp::stop(
        "The value of the polynomial at the center of the bounding sphere ",
        "must be less than the isovalue."
      );
    }
  }

  // isosurface
  ImplicitSurface surface(fun, bounding_sphere);

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

  // convert to EMesh3
  EMesh3 mesh;
  SurfaceMesh::Property_map<SurfaceMesh::Vertex_index, vertex_descriptor> 
    v2vmap = 
      smesh.add_property_map<SurfaceMesh::Vertex_index, vertex_descriptor>(
        "v:v"
      ).first;
  CGAL::copy_face_graph(
    smesh, mesh, CGAL::parameters::vertex_to_vertex_map(v2vmap)
  );

  // normals
  PT3::Differentiate differentiate;
  Poly3 dPx = differentiate(P, 0);
  Poly3 dPy = differentiate(P, 1);
  Poly3 dPz = differentiate(P, 2);
  Normals_map vnormal = 
    mesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
      "v:normal", defaultNormal()
    ).first;
  for(SurfaceMesh::Vertex_index vi : smesh.vertices()) {
    Point_3 p = smesh.point(vi);
    std::list<Real> xyz;
    xyz.push_back(p.x()); 
    xyz.push_back(p.y()); 
    xyz.push_back(p.z());
    FT nx = substitute(dPx, xyz.begin(), xyz.end());
    FT ny = substitute(dPy, xyz.begin(), xyz.end());
    FT nz = substitute(dPz, xyz.begin(), xyz.end());
    Rcpp::NumericVector normal = {-nx, -ny, -nz};
    vnormal[v2vmap[vi]] = normal;
  }

  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}


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

FT spikes_function (Point_3 p) {
  FT x = p.x();
  FT y = p.y();
  FT z = p.z();
  FT x2 = x*x;
  FT y2 = y*y;
  FT z2 = z*z;
  FT val = - exp(-(x2+y2+z2)*2)*1 + exp(-(x2*0.03+y2*0.03+z2)*10) +
    exp(-(x2*0.03+y2+z2*0.03)*10) + exp(-(x2+y2*0.03+z2*0.03)*10);
  return val - 1;
}

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> mandelbulb(double angle_bound, double radius_bound, double distance_bound) {
  Tri tr;            // 3D-Delaunay triangulation
  Cplx2 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  auto fun = [](Point_3 p) {
    FT x = p.x();
    FT y = p.y();
    FT z = p.z();
    FT x2 = x*x;
    FT y2 = y*y;
    FT z2 = z*z;
    return x2 + y2 + z2 - 1;
  };
  ImplicitSurface surface(fun,      // pointer to function
                    Sphere_3(CGAL::ORIGIN, 2.)); // bounding sphere
  // Note that "2." above is the *squared* radius of the bounding sphere!

  // defining meshing criteria
  // CGAL::Surface_mesh_default_criteria_3<Tri> criteria(35.,  // angular bound
  //                                                    0.025,  // radius bound
  //                                                    0.025); // distance bound
  FT ab(angle_bound);
  FT rb(radius_bound);
  FT db(distance_bound); 
  typedef CGAL::Surface_mesh_default_criteria_3<Tri> Crit; 
  Crit criteria(ab,  // angular bound
                                                     rb,  // radius bound
                                                     db); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  SurfaceMesh sm;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

  EMesh3 mesh;
  CGAL::copy_face_graph(sm, mesh);
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> spikes(double angle_bound, double radius_bound, double distance_bound) {
  Tri tr;            // 3D-Delaunay triangulation
  Cplx2 c2t3 (tr);   // 2D-complex in 3D-Delaunay triangulation

  // defining the surface
  ImplicitSurface surface(spikes_function,             // pointer to function
                    Sphere_3(CGAL::ORIGIN, 16.)); // bounding sphere
  // Note that "2." above is the *squared* radius of the bounding sphere!

  // defining meshing criteria
  // CGAL::Surface_mesh_default_criteria_3<Tri> criteria(35.,  // angular bound
  //                                                    0.025,  // radius bound
  //                                                    0.025); // distance bound
  FT ab(angle_bound);
  FT rb(radius_bound);
  FT db(distance_bound); 
  typedef CGAL::Surface_mesh_default_criteria_3<Tri> Crit; 
  Crit criteria(ab,  // angular bound
                                                     rb,  // radius bound
                                                     db); // distance bound
  // meshing surface
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());

  SurfaceMesh sm;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

  EMesh3 mesh;
  CGAL::copy_face_graph(sm, mesh);
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}


// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> Isomesh(
  std::string filename, double isovalue, 
  Rcpp::NumericVector center, double radius,
  double angle_bound, double radius_bound, double distance_bound) {

  Tri tr;            // 3D-Delaunay triangulation
  Cplx2 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation

  // the 'function' is a 3D gray level image
  FT isoval(isovalue);
  Gray_level_image image(CGAL::data_file_path(filename), isoval);

  // Carefully chosen bounding sphere: the center must be inside the
  // surface defined by 'image' and the radius must be high enough so that
  // the sphere actually bounds the whole image.
  Point_3 bounding_sphere_center(center(0), center(1), center(2));
  FT bounding_sphere_squared_radius = radius * radius;
  Sphere_3 bounding_sphere(
    bounding_sphere_center, bounding_sphere_squared_radius
  );

  // definition of the surface, with 10^-5 as relative precision
  Surface_gray surface(image, bounding_sphere, 1e-5);

  // defining meshing criteria
  FT ab(angle_bound);
  FT rb(radius_bound);
  FT db(distance_bound); 
  MeshingCriteria criteria(ab, rb, db);

  // meshing surface, with the "manifold without boundary" algorithm
  CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Manifold_tag());
  SurfaceMesh sm;
  CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, sm);

  EMesh3 mesh;
  CGAL::copy_face_graph(sm, mesh);
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}
