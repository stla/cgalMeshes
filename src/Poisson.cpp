#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/jet_estimate_normals.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/pca_estimate_normals.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
typedef K::Vector_3                                             Vector3;
typedef std::pair<Point3, Vector3>                              P3wn;
typedef CGAL::Polyhedron_3<K, CGAL::Polyhedron_items_with_id_3> Polyhedron;
typedef CGAL::Parallel_if_available_tag Concurrency_tag;


// [[Rcpp::export]]
Rcpp::NumericMatrix jet_normals(
  const Rcpp::NumericMatrix pts, const unsigned nb_neighbors
) {
  const int npoints = pts.ncol();
  std::vector<P3wn> points;
  points.reserve(npoints);
  for(int i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i); 
    points.emplace_back(
      std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)), Vector3(0.0, 0.0, 0.0))
    );
  }

  CGAL::jet_estimate_normals<Concurrency_tag>(
    points, nb_neighbors,
    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
                     .normal_map(CGAL::Second_of_pair_property_map<P3wn>())
  );

  CGAL::mst_orient_normals(
    points, nb_neighbors,
    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
                     .normal_map(CGAL::Second_of_pair_property_map<P3wn>())
  );

  Rcpp::NumericMatrix Normals(3, npoints);
  for(int i = 0; i < npoints; i++) {
    Rcpp::NumericVector normal_i(3);
    const Vector3 normal = points[i].second;
    normal_i(0) = normal.x();
    normal_i(1) = normal.y();
    normal_i(2) = normal.z();
    Normals(Rcpp::_, i) = normal_i;
  }

  return Normals;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix pca_normals(
  const Rcpp::NumericMatrix pts, const unsigned nb_neighbors
) {
  const int npoints = pts.ncol();
  std::vector<P3wn> points;
  points.reserve(npoints);
  for(int i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i); 
    points.emplace_back(
      std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)), Vector3(0.0, 0.0, 0.0))
    );
  }
  
  CGAL::pca_estimate_normals<Concurrency_tag>(
    points, nb_neighbors,
    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
                     .normal_map(CGAL::Second_of_pair_property_map<P3wn>())
  );

  CGAL::mst_orient_normals(
    points, nb_neighbors,
    CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
                     .normal_map(CGAL::Second_of_pair_property_map<P3wn>())
  );

  Rcpp::NumericMatrix Normals(3, npoints);
  for(int i = 0; i < npoints; i++) {
    Rcpp::NumericVector normal_i(3);
    const Vector3 normal = points[i].second;
    normal_i(0) = normal.x();
    normal_i(1) = normal.y();
    normal_i(2) = normal.z();
    Normals(Rcpp::_, i) = normal_i;
  }

  return Normals;
}


// [[Rcpp::export]]
Rcpp::List Poisson_reconstruction_cpp(
  const Rcpp::NumericMatrix pts, const Rcpp::NumericMatrix normals,
  double spacing, 
  const double sm_angle, const double sm_radius, const double sm_distance
) {
  const int npoints = pts.ncol();
  std::vector<P3wn> points;
  points.reserve(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt_i   = pts(Rcpp::_, i); 
    const Rcpp::NumericVector nrml_i = normals(Rcpp::_, i); 
    points.emplace_back(
      std::make_pair(
        Point3(pt_i(0), pt_i(1), pt_i(2)),
        Vector3(nrml_i(0), nrml_i(1), nrml_i(2))
      )
    );
  }
  
  if(spacing == -1.0) {
    spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(
      points, 6, // ??
      CGAL::parameters::point_map(CGAL::First_of_pair_property_map<P3wn>())
    );
  }

  Polyhedron pmesh;
  const bool psr = CGAL::poisson_surface_reconstruction_delaunay(
    points.begin(), points.end(), CGAL::First_of_pair_property_map<P3wn>(),
    CGAL::Second_of_pair_property_map<P3wn>(), pmesh, spacing, sm_angle,
    sm_radius, sm_distance
  );
  
  if(!psr) {
    throw Rcpp::exception("Poisson surface reconstruction has failed.");
  }
  
  Mesh3 mesh;
  int id = 0;
  for(
    Polyhedron::Vertex_iterator vit = pmesh.vertices_begin();
                                vit != pmesh.vertices_end(); ++vit
  ) {
    mesh.add_vertex(vit->point());
    vit->id() = id;
    id++;
  }
  for(
    Polyhedron::Facet_iterator fit = pmesh.facets_begin();
                               fit != pmesh.facets_end(); fit++
  ) {
    int i1 = fit->halfedge()->vertex()->id();
    int i2 = fit->halfedge()->next()->vertex()->id();
    int i3 = fit->halfedge()->opposite()->vertex()->id();
    mesh.add_face(i1, i2, i3); 
  }
  
  EMesh3 emesh;
  CGAL::copy_face_graph(mesh, emesh);
  
  return Rcpp::List::create(
    Rcpp::Named("xptr")    = Rcpp::XPtr<EMesh3>(new EMesh3(emesh), false),
    Rcpp::Named("spacing") = spacing
);
}
