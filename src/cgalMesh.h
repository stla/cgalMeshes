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
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/merge_border_vertices.h>
#include <CGAL/number_utils.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/helpers.h>
//#include <CGAL/boost/graph/Euler_operations.h>
// #include <CGAL/boost/graph/iterator.h>
// #include <CGAL/Iterator_range.h>
//#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/copy.hpp>
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
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
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
typedef boost::graph_traits<EMesh3>::halfedge_descriptor halfedge_descriptor;
typedef EMesh3::Property_map<halfedge_descriptor, std::size_t> Halfedge_index_map;
typedef EMesh3::Property_map<vertex_descriptor, Rcpp::NumericVector> Normals_map;
typedef std::pair<std::map<vertex_descriptor, Rcpp::NumericVector>, bool> MaybeNormalMap;
typedef EMesh3::Property_map<vertex_descriptor, EVector3> CGALnormals_map;
typedef EMesh3::Property_map<vertex_descriptor, std::string> Vcolors_map;
typedef std::pair<std::map<vertex_descriptor, std::string>, bool> MaybeVcolorMap;
typedef EMesh3::Property_map<face_descriptor, std::string> Fcolors_map;
typedef std::pair<std::map<face_descriptor, std::string>, bool> MaybeFcolorMap;
typedef EMesh3::Property_map<vertex_descriptor, double> Vscalars_map;
typedef std::pair<std::map<vertex_descriptor, double>, bool> MaybeVscalarMap;
typedef EMesh3::Property_map<face_descriptor, double> Fscalars_map;
typedef std::pair<std::map<face_descriptor, double>, bool> MaybeFscalarMap;
typedef std::map<face_descriptor, face_descriptor> MapBetweenFaces;
typedef boost::graph_traits<EMesh3>::edge_descriptor edge_descriptor;


typedef CGAL::Advancing_front_surface_reconstruction<> AFS_reconstruction;
typedef AFS_reconstruction::Triangulation_3 AFS_triangulation3;
typedef AFS_reconstruction::Triangulation_data_structure_2 AFS_Tds2;

//typedef CGAL::Face_filtered_graph<EMesh3, Face_index_map, Vertex_index_map, Halfedge_index_map> Filtered_mesh;
typedef CGAL::Face_filtered_graph<EMesh3> Filtered_graph;

typedef boost::graph_traits<Filtered_graph>::vertex_descriptor 
  ffg_vertex_descriptor;
typedef std::map<ffg_vertex_descriptor, vertex_descriptor> 
  MapBetweenVertexDescriptors;
typedef boost::graph_traits<Filtered_graph>::face_descriptor 
  ffg_face_descriptor;
typedef std::map<ffg_face_descriptor, face_descriptor> 
  MapBetweenFaceDescriptors;

struct trivial_edge_predicate {
  trivial_edge_predicate() { }
  bool operator()(const edge_descriptor& e) const {
    return true;
  }
};

struct positive_vertex_scalar {
  positive_vertex_scalar() { }
  positive_vertex_scalar(Vscalars_map vmap) : vscalar(vmap) { }
  bool operator()(const vertex_descriptor& v) const {
    return 0 < vscalar[v];
  }
  Vscalars_map vscalar;
};

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

EMesh3 cloneMesh(EMesh3&, std::vector<std::string>);
void removeProperties(EMesh3&, std::vector<std::string>);
MaybeNormalMap copy_vnormal(EMesh3&);
MaybeVcolorMap copy_vcolor(EMesh3&);
MaybeFcolorMap copy_fcolor(EMesh3&);
void triangulateMesh(EMesh3&);
Rcpp::NumericVector defaultNormal();

template <typename Keytype, typename Valuetype>
std::pair<std::map<Keytype, Valuetype>, bool> copy_prop(
  EMesh3&, std::string
);

Rcpp::List clipping(EMesh3&, EMesh3&, const bool);

template <typename SourceDescriptor, typename TargetDescriptor, typename Valuetype>
void copy_property(
  EMesh3&, EMesh3&, std::map<SourceDescriptor, TargetDescriptor>, std::string 
);
//////////////////////////////////////////

struct ClipVisitor : 
  public PMP::Corefinement::Default_visitor<EMesh3>
{
  // void new_vertex_added(std::size_t i_id, vertex_descriptor v, const EMesh3 & tm) {
  //   Rcpp::Rcout << v << "\n";
  // }
  void before_subface_creations(face_descriptor fsplit, const EMesh3 & tm) {
    *ofaceindex = fsplit;
    if(*is_tm) {
      size_t nf = tm.number_of_faces();
      if((*nfaces).size() >= 1 && nf < (*nfaces).back()) {
        *is_tm = false;
      } else {
        (*nfaces).push_back(nf);
      }
    }
    (*action).push_back("before_subface_creations");
  }
  void after_subface_created(face_descriptor fnew, const EMesh3 & tm) {
    if(*is_tm) {
      (*fmap_tm).insert(std::make_pair(fnew, *ofaceindex));
    } else {
      (*fmap_clipper).insert(std::make_pair(fnew, *ofaceindex));
    }
//    (*pairs).push_back(std::make_pair(fnew, *is_tm));
    (*nfaces2).push_back(tm.number_of_faces());
    (*action).push_back("after_subface_created");
    // Rcpp::Rcout << fnew << "\n";
    // Rcpp::Rcout << tm.number_of_faces() << "\n";
    // Rcpp::Rcout << tm.has_garbage() << "\n";
    // Rcpp::Rcout << "\n";
  }
  // void in_place_operation(PMP::Corefinement::Boolean_operation_type t) {
  //   Rcpp::Rcout << t << "\n";
  // }
  void after_face_copy(
    face_descriptor fsrc, const EMesh3 & tmsrc, 
    face_descriptor ftgt, const EMesh3 & tmtgt
  ) {
    (*ftargets).insert(std::make_pair(ftgt, fsrc));
    (*action).push_back("after_face_copy");
  }

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
  
  ClipVisitor()
    : fmap_tm(new MapBetweenFaces()),
      fmap_clipper(new MapBetweenFaces()),
      ofaceindex(new face_descriptor()),
      nfaces(new std::vector<size_t>()),
      nfaces2(new std::vector<size_t>()),
      ftargets(new MapBetweenFaces()),
//      pairs(new std::vector<std::pair<face_descriptor, bool>>()),
      is_tm(new bool(true)),
      action(new std::vector<std::string>())
  {}
  
  std::shared_ptr<MapBetweenFaces> fmap_tm;
  std::shared_ptr<MapBetweenFaces> fmap_clipper;
  std::shared_ptr<MapBetweenFaces> ftargets;
//  std::shared_ptr<std::vector<std::pair<face_descriptor, bool>>> pairs;
  std::shared_ptr<face_descriptor> ofaceindex;
  std::shared_ptr<std::vector<size_t>> nfaces;
  std::shared_ptr<std::vector<size_t>> nfaces2;
  std::shared_ptr<bool> is_tm;
  std::shared_ptr<std::vector<std::string>> action;
};

struct UnionVisitor : 
  public PMP::Corefinement::Default_visitor<EMesh3>
{
  void before_subface_creations(face_descriptor fsplit, const EMesh3 & tm) {
    *ofaceindex = fsplit;
    (*action).push_back("before_subface_creations");
  }
  void after_subface_created(face_descriptor fnew, const EMesh3 & tm) {
    if(*is_tm) {
      if(int(fnew) - 2 > *fprev) {
        *is_tm = false;
        (*fmap_mesh2).insert(std::make_pair(fnew, *ofaceindex));
      } else {
        *fprev = int(fnew) - 1;
        (*fmap_mesh1).insert(std::make_pair(fnew, *ofaceindex));
      }
    } else {
      (*fmap_mesh2).insert(std::make_pair(fnew, *ofaceindex));
    }
    (*nfaces2).push_back(tm.number_of_faces());
    (*action).push_back("after_subface_created");
  }
  void after_face_copy(
    face_descriptor fsrc, const EMesh3 & tmsrc, 
    face_descriptor ftgt, const EMesh3 & tmtgt
  ) {
    (*fmap_union).insert(std::make_pair(ftgt, fsrc));
    if(*nfaces_umesh1 == -1 && int(fsrc) == *fprev) {
      *nfaces_umesh1 = int(ftgt) + 1;
      Rcpp::Rcout << "FTGT: " << ftgt << "\n";
    }
    (*action).push_back("after_face_copy");
  }
  void after_vertex_copy(
    vertex_descriptor vsrc, const EMesh3 & tmsrc, 
    vertex_descriptor vtgt, const EMesh3 & tmtgt
  ) {
    (*vmap_union).insert(std::make_pair(vtgt, vsrc));
  }
  
  UnionVisitor()
    : fmap_mesh1(new MapBetweenFaces()),
      fmap_mesh2(new MapBetweenFaces()),
      fprev(new int(INT_MAX)),
      nfaces_umesh1(new int(-1)),
      ofaceindex(new face_descriptor()),
      nfaces(new std::vector<size_t>()),
      nfaces2(new std::vector<size_t>()),
      fmap_union(new MapBetweenFaces()),
      vmap_union(new std::map<vertex_descriptor, vertex_descriptor>()),
      is_tm(new bool(true)),
      action(new std::vector<std::string>())
  {}
  
  std::shared_ptr<MapBetweenFaces> fmap_mesh1;
  std::shared_ptr<MapBetweenFaces> fmap_mesh2;
  std::shared_ptr<int> fprev;
  std::shared_ptr<int> nfaces_umesh1;
  std::shared_ptr<MapBetweenFaces> fmap_union;
  std::shared_ptr<std::map<vertex_descriptor, vertex_descriptor>> vmap_union;
  std::shared_ptr<face_descriptor> ofaceindex;
  std::shared_ptr<std::vector<size_t>> nfaces;
  std::shared_ptr<std::vector<size_t>> nfaces2;
  std::shared_ptr<bool> is_tm;
  std::shared_ptr<std::vector<std::string>> action;
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
    : fmap(new MapBetweenFaces()),
      ofaceindex(new face_descriptor())
  {}
  
  std::shared_ptr<MapBetweenFaces> fmap;
  std::shared_ptr<face_descriptor> ofaceindex;
};

struct SoupVisitor : 
  public PMP::Default_orientation_visitor
{
  void non_manifold_vertex(std::size_t vid, std::size_t nb_link_ccs) {
    Message("Detected and duplicated a non-manifold vertex.");
  }

  void non_manifold_edge(std::size_t id1, std::size_t id2, std::size_t nb_polygons) {
    Message("Detected a non-manifold edge.");
  }

  // void polygon_orientation_reversed(std::size_t id) {
  //   Rcpp::Rcout << "ORIENTATION REVERSED: " << id << "\n";
  // }

  SoupVisitor() {} 
};

// template <class MeshT>
//using Visitor = PMP::PMPCorefinementVisitor;

//typedef struct PMPCorefinementVisitor Visitor;

//template <typename>
//void Visitor::new_vertex_added(std::size_t, vertex_descriptor, const EMesh3&);

//void Visitor<EMesh3>::new_vertex_added(std::size_t, vertex_descriptor, const EMesh3&);
