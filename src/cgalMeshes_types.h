#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef EK::Point_3                                       EPoint3;
typedef CGAL::Surface_mesh<EPoint3>                       EMesh3;

struct MyMesh {
  EMesh3 mesh; 
  Rcpp::Nullable<Rcpp::NumericMatrix> normals = R_NilValue;
  Rcpp::Nullable<Rcpp::StringVector> vcolors = R_NilValue;
  Rcpp::Nullable<Rcpp::StringVector> fcolors = R_NilValue;
};
#define MYMESH_OPEN [&] {MyMesh xxx = {};  
#define MYMESH_CLOSE ; return xxx; }()
#define MYMESH(x) MYMESH_OPEN x MYMESH_CLOSE
