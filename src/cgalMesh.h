#ifndef _HEADER_
#define _HEADER_
#endif

#include <Rcpp.h>

#define CGAL_EIGEN3_ENABLED 1

#include "cgalMeshes_types.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/number_utils.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
//#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <boost/graph/connected_components.hpp>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
// #include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <locale>  // tolower
#include <CGAL/IO/io.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
//////
//#include <CGAL/Surface_mesh_parameterization/Discrete_conformal_map_parameterizer_3.h>
//#include <CGAL/Surface_mesh_parameterization/Error_code.h>
//#include <CGAL/Surface_mesh_parameterization/parameterize.h>
//#include <CGAL/boost/graph/Seam_mesh.h>

//#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>

//#include <CGAL/Unique_hash_map.h>

//#include <CGAL/Surface_mesh_parameterization/Orbifold_Tutte_parameterizer_3.h>
//#include <CGAL/boost/graph/properties.h>


//namespace SMP = CGAL::Surface_mesh_parameterization;



// -------------------------------------------------------------------------- //
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point3;
typedef CGAL::Surface_mesh<Point3> Mesh3;
typedef EK::Vector_3 EVector3;
typedef CGAL::Nef_polyhedron_3<EK, CGAL::SNC_indexed_items> NefPol;
typedef CGAL::Polyhedron_3<EK> EPolyhedron;

typedef CGAL::Triangle_3<EK> Triangle;
typedef CGAL::Tetrahedron_3<EK> Tetrahedron;

//typedef CGAL::SM_index_pmap<EPoint3, boost::graph_traits<EMesh3>::vertex_descriptor> boostmap;
typedef boost::graph_traits<EMesh3>::vertex_descriptor vertex_descriptor;
typedef EMesh3::Property_map<vertex_descriptor, double> Vertex_distance_map;
typedef EMesh3::Property_map<vertex_descriptor, std::size_t> Vertex_index_map;
typedef boost::graph_traits<EMesh3>::face_descriptor face_descriptor;
typedef EMesh3::Property_map<face_descriptor, std::size_t> Face_index_map;
typedef boost::graph_traits<Mesh3>::halfedge_descriptor halfedge_descriptor;
typedef EMesh3::Property_map<halfedge_descriptor, std::size_t> Halfedge_index_map;
typedef EMesh3::Property_map<vertex_descriptor, Rcpp::NumericVector> Normals_map;
typedef EMesh3::Property_map<vertex_descriptor, std::string> Vcolors_map;
typedef EMesh3::Property_map<face_descriptor, std::string> Fcolors_map;

typedef CGAL::Advancing_front_surface_reconstruction<> AFS_reconstruction;
typedef AFS_reconstruction::Triangulation_3 AFS_triangulation3;
typedef AFS_reconstruction::Triangulation_data_structure_2 AFS_Tds2;

//typedef CGAL::Face_filtered_graph<EMesh3, Face_index_map, Vertex_index_map, Halfedge_index_map> Filtered_graph;
typedef CGAL::Face_filtered_graph<EMesh3> Filtered_graph;

///////////////
// typedef boost::graph_traits<Mesh3>::halfedge_descriptor halfedge_descriptor;
// typedef K::Point_2 Point2;
// typedef boost::graph_traits<Mesh3>::vertex_descriptor v_descriptor;
// typedef boost::graph_traits<Mesh3>::edge_descriptor edge_descriptor;
// 
// typedef Mesh3 PolyMesh;
// typedef boost::graph_traits<PolyMesh>::edge_descriptor SM_edge_descriptor;
// typedef boost::graph_traits<PolyMesh>::halfedge_descriptor SM_halfedge_descriptor;
// typedef boost::graph_traits<PolyMesh>::vertex_descriptor SM_vertex_descriptor;
// 
// typedef CGAL::Unique_hash_map<SM_halfedge_descriptor, Point2> UV_uhm;
// typedef CGAL::Unique_hash_map<SM_edge_descriptor, bool> Seam_edge_uhm;
// typedef CGAL::Unique_hash_map<SM_vertex_descriptor, bool> Seam_vertex_uhm;
// 
// typedef boost::associative_property_map<UV_uhm> Seam_UV_pmap;
// typedef boost::associative_property_map<Seam_edge_uhm> Seam_edge_pmap;
// typedef boost::associative_property_map<Seam_vertex_uhm> Seam_vertex_pmap;
// 
// typedef Mesh3::Property_map<SM_edge_descriptor, bool>           edge_pmap;
// typedef Mesh3::Property_map<SM_vertex_descriptor, bool>         vertex_pmap;
// 
// typedef CGAL::Seam_mesh<PolyMesh, edge_pmap, vertex_pmap> SMesh;
// 
// typedef CGAL::Seam_mesh<PolyMesh, Seam_edge_pmap, Seam_vertex_pmap> Mesh;
// 
// typedef boost::graph_traits<Mesh>::vertex_descriptor seam_vertex_descriptor;
// typedef boost::graph_traits<Mesh>::halfedge_descriptor seam_halfedge_descriptor;
// typedef boost::graph_traits<Mesh>::face_descriptor seam_face_descriptor;
// 
// typedef Mesh3::Property_map<halfedge_descriptor, Point2>      UV_pmap;


// -------------------------------------------------------------------------- //
namespace PMP = CGAL::Polygon_mesh_processing;

// -------------------------------------------------------------------------- //
template <typename MeshT, typename PointT>
MeshT csoup2mesh(std::vector<PointT>, std::vector<std::vector<int>>, const bool);

std::vector<std::vector<int>> list_to_faces(const Rcpp::List);

template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix);

template <typename KernelT, typename MeshT, typename PointT>
Rcpp::DataFrame getEdges(MeshT&);

//Rcpp::NumericMatrix getEKNormals(EMesh3);

Rcpp::NumericMatrix getVertices_EK(EMesh3&);
Rcpp::List RSurfEKMesh(EMesh3&, const bool);
Rcpp::List RSurfEKMesh2(EMesh3&, const bool, const int);

template <typename MeshT>
Rcpp::List getFaces(MeshT&);

void Message(std::string);

EMesh3 readMeshFile(const std::string);
void writeMeshFile(const std::string, const int, const bool, EMesh3&);

EMesh3 dualMesh(EMesh3&);

EMesh3 makeMesh(const Rcpp::NumericMatrix,
                const Rcpp::List,
                const bool,
                const Rcpp::Nullable<Rcpp::NumericMatrix>&,
                const Rcpp::Nullable<Rcpp::StringVector>&,
                const Rcpp::Nullable<Rcpp::StringVector>&);

EMesh3 cloneMesh(EMesh3&, const bool, const bool, const bool);
void removeProperties(EMesh3&, const bool, const bool, const bool);
std::pair<std::map<vertex_descriptor, Rcpp::NumericVector>, bool> copy_vnormal(EMesh3&);
std::pair<std::map<vertex_descriptor, std::string>, bool> copy_vcolor(EMesh3&);
std::pair<std::map<face_descriptor, std::string>, bool> copy_fcolor(EMesh3&);

void clipping(EMesh3&, EMesh3&, const bool);
//////////////////////////////////////////
void new_vertex_added(std::size_t, vertex_descriptor, const EMesh3&);

struct MyVisitor : 
  public PMP::Corefinement::Default_visitor<EMesh3>
{
  // void new_vertex_added(std::size_t i_id, vertex_descriptor v, const EMesh3 & tm) {
  //   Rcpp::Rcout << v << "\n";
  // }
  void before_subface_creations(face_descriptor fsplit, const EMesh3 & tm) {
    *ofaceindex = fsplit;
    // Rcpp::Rcout << j++ << "\n";
    // Rcpp::Rcout << tm.has_garbage() << "\n";
    // Rcpp::Rcout << "\n";
    //Rcpp::Rcout << tm.number_of_faces() << "\n";
  }
  void after_subface_created(face_descriptor fnew, const EMesh3 & tm) {
    (*fmap).insert(std::make_pair(fnew, *ofaceindex));
    // Rcpp::Rcout << fnew << "\n";
    // Rcpp::Rcout << tm.number_of_faces() << "\n";
    // Rcpp::Rcout << tm.has_garbage() << "\n";
    // Rcpp::Rcout << "\n";
  }
  // void in_place_operation(PMP::Corefinement::Boolean_operation_type t) {
  //   Rcpp::Rcout << t << "\n";
  // }
  // void after_face_copy(face_descriptor fsrc, const EMesh3 & tmsrc, face_descriptor ftgt, const EMesh3 & tmtgt) {
  //   (*vmap).insert(std::make_pair(fsrc, ftgt));
  // }
  // void after_edge_duplicated(halfedge_descriptor hsrc, halfedge_descriptor hnew, const EMesh3 & tm) {
  //   (*vmap).insert(std::make_pair(hsrc, hnew));
  // }
  // void intersection_edge_copy(halfedge_descriptor hsrc1, const EMesh3 & tmsrc1, halfedge_descriptor hsrc2, const EMesh3 & tmsrc2, halfedge_descriptor htgt, const EMesh3 & tmtgt) {
  //   (*vmap).insert(std::make_pair(hsrc1, htgt));
  // }
  // void intersection_point_detected(std::size_t i_id, int sdim, halfedge_descriptor h_f, halfedge_descriptor h_e, const EMesh3 & tm_f, const EMesh3 & tm_e, bool is_target_coplanar, bool is_source_coplanar) {
  //   (*vmap).insert(std::make_pair(h_f, sdim));
  //   (*i).push_back(tm_f.number_of_faces());
  // }
  
  MyVisitor()
    : fmap(new std::map<face_descriptor, face_descriptor>()),
      ofaceindex(new face_descriptor())
  {}
  
  std::shared_ptr<std::map<face_descriptor, face_descriptor>> fmap;
  std::shared_ptr<face_descriptor> ofaceindex;
};

struct UnionVisitor : 
  public PMP::Corefinement::Default_visitor<EMesh3>
{
  void before_face_copy(face_descriptor fsrc, const EMesh3 & tmsrc, const EMesh3 & tmtgt) {
    Rcpp::Rcout << fsrc;
  }

  void after_face_copy(face_descriptor fsrc, const EMesh3 & tmsrc, face_descriptor ftgt, const EMesh3 & tmtgt) {
    (*fmap).insert(std::make_pair(fsrc, ftgt));
  }
  
  UnionVisitor()
    : fmap(new std::map<face_descriptor, face_descriptor>()),
      ofaceindex(new face_descriptor())
  {}
  
  std::shared_ptr<std::map<face_descriptor, face_descriptor>> fmap;
  std::shared_ptr<face_descriptor> ofaceindex;
};

struct TriangulateVisitor : 
  public PMP::Triangulate_faces::Default_visitor<EMesh3>
{
  void before_subface_creations(face_descriptor fsplit) {
    *ofaceindex = fsplit;
  }
  void after_subface_created(face_descriptor fnew) {
    (*fmap).insert(std::make_pair(fnew, *ofaceindex));
  }
  
  TriangulateVisitor()
    : fmap(new std::map<face_descriptor, face_descriptor>()),
      ofaceindex(new face_descriptor())
  {}
  
  std::shared_ptr<std::map<face_descriptor, face_descriptor>> fmap;
  std::shared_ptr<face_descriptor> ofaceindex;
};

// template <class MeshT>
//using Visitor = PMP::PMPCorefinementVisitor;

//typedef struct PMPCorefinementVisitor Visitor;

//template <typename>
//void Visitor::new_vertex_added(std::size_t, vertex_descriptor, const EMesh3&);

//void Visitor<EMesh3>::new_vertex_added(std::size_t, vertex_descriptor, const EMesh3&);
