#ifndef _HEADER_
#define _HEADER_
#endif

#include <Rcpp.h>
#include <RcppColors.h>

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
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/number_utils.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/boost/graph/helpers.h>
#include <boost/graph/connected_components.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/copy.hpp>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <locale>  // tolower
#include <CGAL/IO/io.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
//#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Triangle_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Plane_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Iso_cuboid_3.h>
// #include <CGAL/Surface_mesh_default_triangulation_3.h>
// #include <CGAL/Complex_2_in_triangulation_3.h>
// #include <CGAL/make_surface_mesh.h>
// #include <CGAL/Implicit_surface_3.h>
// #include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
// #include <CGAL/Surface_mesh_default_criteria_3.h>
// #include <CGAL/Complex_2_in_triangulation_3.h>
// #include <CGAL/Polynomial.h>
// #include <CGAL/Polynomial_traits_d.h>
// #include <CGAL/Polynomial_type_generator.h>
// #include <CGAL/polynomial_utils.h>
#include <CGAL/Subdivision_method_3/subdivision_methods_3.h>
#include <CGAL/boost/graph/generators.h>
// #include <CGAL/convex_hull_3.h>
// #include <CGAL/Convex_hull_3/dual/halfspace_intersection_3.h>
#include <CGAL/Side_of_triangle_mesh.h>

// -------------------------------------------------------------------------- //
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3                                          Point3;
typedef CGAL::Surface_mesh<Point3>                          Mesh3;
typedef EK::Vector_3                                        EVector3;
typedef EK::Plane_3                                         EPlane3;
typedef CGAL::Bbox_3                                        Bbox3;
typedef CGAL::Iso_cuboid_3<EK>                              IsoCuboid3;
typedef CGAL::Nef_polyhedron_3<EK, CGAL::SNC_indexed_items> NefPol;
typedef CGAL::Polyhedron_3<EK>                              EPolyhedron;

typedef CGAL::Triangle_3<EK>    Triangle;
typedef CGAL::Tetrahedron_3<EK> Tetrahedron;

typedef boost::graph_traits<Mesh3>::vertex_descriptor                     vxdescr;
typedef boost::graph_traits<EMesh3>::vertex_descriptor                    vertex_descriptor;
typedef EMesh3::Property_map<vertex_descriptor, double>                   Vertex_distance_map;
typedef EMesh3::Property_map<vertex_descriptor, std::size_t>              Vertex_index_map;
typedef boost::graph_traits<EMesh3>::face_descriptor                      face_descriptor;
typedef EMesh3::Property_map<face_descriptor, std::size_t>                Face_index_map;
typedef boost::graph_traits<EMesh3>::halfedge_descriptor                  halfedge_descriptor;
typedef EMesh3::Property_map<halfedge_descriptor, std::size_t>            Halfedge_index_map;
typedef EMesh3::Property_map<vertex_descriptor, Rcpp::NumericVector>      Normals_map;
typedef std::pair<std::map<vertex_descriptor, Rcpp::NumericVector>, bool> MaybeNormalMap;
typedef EMesh3::Property_map<vertex_descriptor, EVector3>                 CGALnormals_map;
typedef EMesh3::Property_map<vertex_descriptor, std::string>              Vcolors_map;
typedef std::pair<std::map<vertex_descriptor, std::string>, bool>         MaybeVcolorMap;
typedef EMesh3::Property_map<face_descriptor, std::string>                Fcolors_map;
typedef std::pair<std::map<face_descriptor, std::string>, bool>           MaybeFcolorMap;
typedef EMesh3::Property_map<vertex_descriptor, double>                   Vscalars_map;
typedef std::pair<std::map<vertex_descriptor, double>, bool>              MaybeVscalarMap;
typedef EMesh3::Property_map<face_descriptor, double>                     Fscalars_map;
typedef std::pair<std::map<face_descriptor, double>, bool>                MaybeFscalarMap;
typedef std::map<face_descriptor, face_descriptor>                        MapBetweenFaces;
typedef std::map<vertex_descriptor, vertex_descriptor>                    MapBetweenVertices;
typedef boost::graph_traits<EMesh3>::edge_descriptor                      edge_descriptor;

typedef CGAL::Advancing_front_surface_reconstruction<>     AFS_reconstruction;
typedef AFS_reconstruction::Triangulation_3                AFS_triangulation3;
typedef AFS_reconstruction::Triangulation_data_structure_2 AFS_Tds2;

typedef CGAL::Scale_space_surface_reconstruction_3<K>                SSS_reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother<K> SSS_smoother;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher<K>    SSS_mesher;
typedef SSS_reconstruction::Facet_const_iterator                      SSS_facet_iterator;
typedef SSS_reconstruction::Point_const_iterator                      SSS_point_iterator;

typedef CGAL::Face_filtered_graph<EMesh3>                      Filtered_graph;
typedef boost::graph_traits<Filtered_graph>::vertex_descriptor ffg_vertex_descriptor;
typedef std::map<ffg_vertex_descriptor, vertex_descriptor>     MapBetweenVertexDescriptors;
typedef boost::graph_traits<Filtered_graph>::face_descriptor   ffg_face_descriptor;
typedef std::map<ffg_face_descriptor, face_descriptor>         MapBetweenFaceDescriptors;

typedef CGAL::IO::Color Color;

// typedef CGAL::Surface_mesh_default_triangulation_3 Tri;
// typedef CGAL::Surface_mesh_default_criteria_3<Tri> MeshingCriteria; 
// typedef CGAL::Complex_2_in_triangulation_3<Tri>    Cplx2;
// typedef Tri::Geom_traits                           GT;
// typedef GT::Sphere_3                               Sphere_3;
// typedef GT::Point_3                                Point_3;
// typedef GT::FT                                     FT;
// typedef                                            FT (*Function)(Point_3);
// typedef CGAL::Implicit_surface_3<GT, Function>     ImplicitSurface;
// typedef CGAL::Surface_mesh<Point_3>                SurfaceMesh;
// 
// typedef CGAL::Polynomial_type_generator<FT, 3>::Type     Poly3;
// typedef CGAL::Polynomial_traits_d<Poly3>                 PT3;
// typedef PT3::Innermost_coefficient_type                  Real;
// 
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

Rcpp::NumericMatrix getVertices_EK(EMesh3&);
Rcpp::List RSurfEKMesh(EMesh3&, const bool);
Rcpp::List RSurfEKMesh2(EMesh3&, const bool, const int);

template <typename MeshT>
Rcpp::List getFaces(MeshT&);

void Message(std::string);

EMesh3 readMeshFile(const std::string, bool);
void writeMeshFile(
  const std::string, const int, const bool, std::string, EMesh3&
);

EMesh3 dualMesh(EMesh3&);

EMesh3 makeMesh(const Rcpp::NumericMatrix,
                const Rcpp::List,
                const bool,
                const Rcpp::Nullable<Rcpp::NumericMatrix>&,
                const Rcpp::Nullable<Rcpp::StringVector>&,
                const Rcpp::Nullable<Rcpp::StringVector>&);

EMesh3 cloneMesh(EMesh3&, std::vector<std::string>);
void removeProperties(EMesh3&, std::vector<std::string>);
void triangulateMesh(EMesh3&);
Rcpp::NumericVector defaultNormal();

template <typename Keytype, typename Valuetype>
std::pair<std::map<Keytype, Valuetype>, bool> copy_prop(
  EMesh3&, std::string
);

template <typename Keytype, typename Valuetype>
void removeProperty(EMesh3&, std::string);

template <typename SourceDescriptor, typename TargetDescriptor, typename Valuetype>
void copy_property(
  EMesh3&, EMesh3&, std::map<SourceDescriptor, TargetDescriptor>, std::string 
);

//////////////////////////////////////////

struct ClipVisitor : 
  public PMP::Corefinement::Default_visitor<EMesh3>
{
  void before_subface_creations(face_descriptor fsplit, const EMesh3 & tm) {
    *ofaceindex = fsplit;
  }

  void after_subface_created(face_descriptor fnew, const EMesh3 & tm) {
    if(*is_tm) {
      if(tm.property_map<face_descriptor, std::size_t>("f:i").second) {
        (*fmap_tm).insert(std::make_pair(fnew, *ofaceindex));
      } else {
        *is_tm = false;
        (*fmap_clipper).insert(std::make_pair(fnew, *ofaceindex));
      }
    } else {
      (*fmap_clipper).insert(std::make_pair(fnew, *ofaceindex));
    }
  }

  void after_face_copy(
    face_descriptor fsrc, const EMesh3 & tmsrc, 
    face_descriptor ftgt, const EMesh3 & tmtgt
  ) {
    (*ftargets).insert(std::make_pair(ftgt, fsrc));
  }
  
  ClipVisitor()
    : ofaceindex(new face_descriptor()),
      fmap_tm(new MapBetweenFaces()),
      fmap_clipper(new MapBetweenFaces()),
      ftargets(new MapBetweenFaces()),
      is_tm(new bool(true))
  {}
  
  std::shared_ptr<face_descriptor> ofaceindex;
  std::shared_ptr<MapBetweenFaces> fmap_tm;
  std::shared_ptr<MapBetweenFaces> fmap_clipper;
  std::shared_ptr<MapBetweenFaces> ftargets;
  std::shared_ptr<bool> is_tm;
};


struct DifferenceVisitor : 
  public PMP::Corefinement::Default_visitor<EMesh3>
{
  void before_subface_creations(face_descriptor fsplit, const EMesh3 & tm) {
    *ofaceindex = fsplit;
  }

  void after_subface_created(face_descriptor fnew, const EMesh3 & tm) {
    if(*is_mesh1) {
      if(tm.property_map<face_descriptor, std::size_t>("f:i").second) {
        (*fmap_mesh1).insert(std::make_pair(fnew, *ofaceindex));
      } else {
        *is_mesh1 = false;
        (*fmap_mesh2).insert(std::make_pair(fnew, *ofaceindex));
      }
    } else {
      (*fmap_mesh2).insert(std::make_pair(fnew, *ofaceindex));
    }
  }

  void after_face_copy(
    face_descriptor fsrc, const EMesh3 & tmsrc, 
    face_descriptor ftgt, const EMesh3 & tmtgt
  ) {
    if(*is_mesh1src) {
      *is_mesh1src = 
        tmsrc.property_map<face_descriptor, std::size_t>("f:i").second;
      (*nfaces_dmesh1)++;
    }
    (*fmap_difference).insert(std::make_pair(ftgt, fsrc));
  }

  void after_vertex_copy(
    vertex_descriptor vsrc, const EMesh3 & tmsrc, 
    vertex_descriptor vtgt, const EMesh3 & tmtgt
  ) {
    (*vmap_difference).insert(std::make_pair(vtgt, vsrc));
  }
  
  DifferenceVisitor()
    : fmap_mesh1(new MapBetweenFaces()),
      fmap_mesh2(new MapBetweenFaces()),
      ofaceindex(new face_descriptor()),
      fmap_difference(new MapBetweenFaces()),
      nfaces_dmesh1(new int(-1)),
      vmap_difference(new MapBetweenVertices()),
      is_mesh1(new bool(true)),
      is_mesh1src(new bool(true))
  {}
  
  std::shared_ptr<MapBetweenFaces> fmap_mesh1;
  std::shared_ptr<MapBetweenFaces> fmap_mesh2;
  std::shared_ptr<face_descriptor> ofaceindex;
  std::shared_ptr<MapBetweenFaces> fmap_difference;
  std::shared_ptr<int> nfaces_dmesh1;
  std::shared_ptr<MapBetweenVertices> vmap_difference;
  std::shared_ptr<bool> is_mesh1;
  std::shared_ptr<bool> is_mesh1src;
};


struct UnionVisitor : 
  public PMP::Corefinement::Default_visitor<EMesh3>
{
  void before_subface_creations(face_descriptor fsplit, const EMesh3 & tm) {
    *ofaceindex = fsplit;
  }

  void after_subface_created(face_descriptor fnew, const EMesh3 & tm) {
    if(*is_mesh1) {
      if(tm.property_map<face_descriptor, std::size_t>("f:i").second) {
        (*fmap_mesh1).insert(std::make_pair(fnew, *ofaceindex));
        *fprev = int(fnew) - 1;
      } else {
        *is_mesh1 = false;
        (*fmap_mesh2).insert(std::make_pair(fnew, *ofaceindex));
      }
    } else {
      (*fmap_mesh2).insert(std::make_pair(fnew, *ofaceindex));
    }
  }

  void after_face_copy(
    face_descriptor fsrc, const EMesh3 & tmsrc, 
    face_descriptor ftgt, const EMesh3 & tmtgt
  ) {
    (*fmap_union).insert(std::make_pair(ftgt, fsrc));
    if(*nfaces_umesh1 == -1 && int(fsrc) == *fprev) {
      *nfaces_umesh1 = int(ftgt) + 1;
    }
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
      fprev(new int(0)),
      nfaces_umesh1(new int(-1)),
      ofaceindex(new face_descriptor()),
      fmap_union(new MapBetweenFaces()),
      vmap_union(new std::map<vertex_descriptor, vertex_descriptor>()),
      is_mesh1(new bool(true))
  {}
  
  std::shared_ptr<MapBetweenFaces> fmap_mesh1;
  std::shared_ptr<MapBetweenFaces> fmap_mesh2;
  std::shared_ptr<int> fprev;
  std::shared_ptr<int> nfaces_umesh1;
  std::shared_ptr<face_descriptor> ofaceindex;
  std::shared_ptr<MapBetweenFaces> fmap_union;
  std::shared_ptr<std::map<vertex_descriptor, vertex_descriptor>> vmap_union;
  std::shared_ptr<bool> is_mesh1;
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

  SoupVisitor() {} 
};
