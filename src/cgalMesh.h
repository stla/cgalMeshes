#ifndef _HEADER_
#define _HEADER_
#endif

#include <Rcpp.h>

#define CGAL_EIGEN3_ENABLED 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/number_utils.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
// #include <CGAL/Surface_mesh/IO/PLY.h>
// #include <locale>  // tolower
// #include <CGAL/IO/io.h>

// -------------------------------------------------------------------------- //
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
//typedef K::Point_3 Point3;
typedef EK::Point_3 EPoint3;
typedef CGAL::Surface_mesh<EPoint3> EMesh3;
typedef EK::Vector_3 EVector3;

// -------------------------------------------------------------------------- //
namespace PMP = CGAL::Polygon_mesh_processing;

// -------------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT csoup2mesh(std::vector<PointT>, std::vector<std::vector<int>>, const bool);

std::vector<std::vector<int>> list_to_faces(const Rcpp::List);




template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix);

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame getEdges(MeshT);

//Rcpp::NumericMatrix getEKNormals(EMesh3);

Rcpp::List RSurfEKMesh(EMesh3, const bool);
Rcpp::List RSurfEKMesh2(EMesh3, const bool, const int);

void Message(std::string);

EMesh3 readMeshFile(const std::string, const bool);
void writeMeshFile(const std::string, const int precision, EMesh3);