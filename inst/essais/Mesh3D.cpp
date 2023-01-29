#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
//#include <CGAL/refine_mesh_3.h>
#include <CGAL/IO/facets_in_complex_3_to_triangle_mesh.h>


// Domain
//typedef CGAL::Surface_mesh<K::Point_3> Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_3<Mesh3, K> MeshDomain;

#ifdef CGAL_CONCURRENT_MESH_3
typedef CGAL::Parallel_tag   Concurrency_tag;
#else
typedef CGAL::Sequential_tag Concurrency_tag;
#endif

// Triangulation
typedef CGAL::Mesh_triangulation_3<MeshDomain, CGAL::Default, Concurrency_tag>::type Tr3;

typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr3> Cplx3;

// Criteria
typedef CGAL::Mesh_criteria_3<Tr3> MeshCriteria3;


// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> mmesh3d(double timeLimit, size_t maxIterations) {

  Mesh3 mesh;
  std::ifstream input("/home/stla/Documents/R/MyPackages/cgalMeshes/inst/trash/mesh.off");
  input >> mesh;
  if(input.fail()){
    std::cerr << "Error: Cannot read file " << std::endl;
    Rcpp::stop("");
  }
  Rcpp::Rcout << "file readed\n";
  input.close();
  Rcpp::Rcout << "file closed\n";

  if (!CGAL::is_triangle_mesh(mesh)){
    std::cerr << "Input geometry is not triangulated." << std::endl;
    Rcpp::stop("");
  }

  // Create domain
  MeshDomain domain(mesh);
    Rcpp::Rcout << "domain created\n";


  // Mesh criteria (no cell_size set)
  MeshCriteria3 criteria(
    CGAL::parameters::facet_angle= 5.0,
    CGAL::parameters::facet_size = 1.15,
    CGAL::parameters::facet_distance= 1.8,
    CGAL::parameters::cell_radius_edge_ratio = 3.0);
  Rcpp::Rcout << "criterai defined\n";

  // Mesh generation
  Cplx3 cplx3 = CGAL::make_mesh_3<Cplx3>(
    domain, criteria, CGAL::parameters::lloyd(
      CGAL::parameters::time_limit = timeLimit, 
      CGAL::parameters::max_iteration_number = maxIterations
    )
  );
  Rcpp::Rcout << "make_mesh_3 done\n";

  Mesh3 outmesh;
  CGAL::facets_in_complex_3_to_triangle_mesh(cplx3, outmesh);

  // Output
  EMesh3 emesh;
  CGAL::copy_face_graph(outmesh, emesh);

  return Rcpp::XPtr<EMesh3>(new EMesh3(emesh), false);
}

