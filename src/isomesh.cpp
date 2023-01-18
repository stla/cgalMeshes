#ifndef _HEADER_
#include "cgalMesh.h"
#endif

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
Rcpp::XPtr<EMesh3> AlgebraicMesh(
  Rcpp::IntegerMatrix powers, Rcpp::NumericVector coeffs, double isolevel,
  Rcpp::NumericVector sphereCenter, double sphereRadius,
  double angle_bound, double radius_bound, double distance_bound, 
  double error_bound
) {
  Tri tr;             // 3D-Delaunay triangulation
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
      std::string msg = "The value of the polynomial at the center of the "; 
      msg += "bounding sphere must be less than the isovalue.";
      Rcpp::stop(msg);
    }
  }

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


// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> VoxelToMesh(
  std::string filename, double isovalue, 
  Rcpp::NumericVector center, double radius,
  double angle_bound, double radius_bound, double distance_bound, 
  double error_bound
) {

  Tri tr;            // 3D-Delaunay triangulation
  Cplx2 cplx2(tr);   // 2D-complex in 3D-Delaunay triangulation

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

  // definition of the surface
  FT eb(error_bound);
  Surface_gray surface(image, bounding_sphere, eb);

  // defining meshing criteria
  FT ab(angle_bound);
  FT rb(radius_bound);
  FT db(distance_bound); 
  MeshingCriteria criteria(ab, rb, db);

  // meshing surface, with the "manifold with boundary" algorithm
  CGAL::make_surface_mesh(
    cplx2, surface, criteria, CGAL::Manifold_with_boundary_tag()
  );
  SurfaceMesh smesh;
  CGAL::facets_in_complex_2_to_triangle_mesh(cplx2, smesh);

  EMesh3 mesh;
  CGAL::copy_face_graph(smesh, mesh);
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}
