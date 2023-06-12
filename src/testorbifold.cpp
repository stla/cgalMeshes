#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Surface_mesh_parameterization/Orbifold_Tutte_parameterizer_3.h>
//#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/properties.h>
// #include <CGAL/Timer.h>
// #include <unordered_map>
// #include <fstream>
// #include <iostream>
// #include <list>
// #include <string>
// #include <utility>
// #include <vector>

// typedef CGAL::Simple_cartesian<double>            Kernel;
// typedef Kernel::Point_2                           Point_2;
// typedef Kernel::Point_3                           Point_3;
// typedef CGAL::Surface_mesh<Kernel::Point_3>       SurfaceMesh;
// typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     SM_vertex_descriptor;
// typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   SM_halfedge_descriptor;
// typedef boost::graph_traits<SurfaceMesh>::edge_descriptor       SM_edge_descriptor;
typedef Mesh3::Property_map<edescr, bool>           Seam_edge_pmap;
typedef Mesh3::Property_map<vxdescr, bool>         Seam_vertex_pmap;
typedef CGAL::Seam_mesh<Mesh3, Seam_edge_pmap, Seam_vertex_pmap>  Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor                    svertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor                  shalfedge_descriptor;
typedef Mesh3::Property_map<hgdescr, Point2>      hUV_pmap;

// [[Rcpp::export]]
Rcpp::NumericMatrix testo(std::string filename)
{
  Mesh3 sm;
  if(!CGAL::IO::read_polygon_mesh(filename, sm))
  {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  // Selection file that contains the cones and possibly the path between cones
  // -- the first line for the cones indices
  // -- the second line must be empty
  // -- the third line optionally provides the seam edges indices as 'e11 e12 e21 e22 e31 e32' etc.
  const char* cone_filename = "cones.selection.txt";
  // Read the cones and compute their corresponding vertex_descriptor in the underlying mesh 'sm'
  std::vector<vxdescr> cone_sm_vds;
  SMP::read_cones<Mesh3>(sm, cone_filename, std::back_inserter(cone_sm_vds));
  // Two property maps to store the seam edges and vertices
  Seam_edge_pmap seam_edge_pm = sm.add_property_map<edescr, bool>("e:on_seam", false).first;
  Seam_vertex_pmap seam_vertex_pm = sm.add_property_map<vxdescr, bool>("v:on_seam",false).first;
  // The seam mesh
  Mesh mesh(sm, seam_edge_pm, seam_vertex_pm);
  // If provided, use the path between cones to create a seam mesh
  hgdescr smhd = mesh.add_seams(cone_filename);
  // If not provided, compute the paths using shortest paths
  if(smhd == hgdescr() ) {
    std::cout << "No seams given in input, computing the shortest paths between consecutive cones" << std::endl;
    std::list<edescr> seam_edges;
    SMP::compute_shortest_paths_between_cones(sm, cone_sm_vds.begin(), cone_sm_vds.end(), seam_edges);
    // Add the seams to the seam mesh
    for(edescr e : seam_edges) {
      mesh.add_seam(source(e, sm), target(e, sm));
    }
  }
  std::cout << mesh.number_of_seam_edges() << " seam edges in input" << std::endl;
  // Index map of the seam mesh (assuming a single connected component so far)
  typedef std::unordered_map<svertex_descriptor, int> Indices;
  Indices indices;
  boost::associative_property_map<Indices> vimap(indices);
  int counter = 0;
  for(svertex_descriptor vd : vertices(mesh)) {
    put(vimap, vd, counter++);
  }
  // Mark the cones in the seam mesh
  std::unordered_map<svertex_descriptor, SMP::Cone_type> cmap;
  SMP::locate_cones(mesh, cone_sm_vds.begin(), cone_sm_vds.end(), cmap);
  // The 2D points of the uv parametrisation will be written into this map
  // Note that this is a halfedge property map, and that uv values
  // are only stored for the canonical halfedges representing a vertex
  hUV_pmap uvmap = sm.add_property_map<hgdescr, Point2>("h:uv").first;
  // Parameterizer
  typedef SMP::Orbifold_Tutte_parameterizer_3<Mesh>         Parameterizer;
  Parameterizer parameterizer(SMP::Triangle, SMP::Cotangent);
  // a halfedge on the (possibly virtual) border
  // only used in output (will also be used to handle multiple connected components in the future)
  shalfedge_descriptor bhd = PMP::longest_border(mesh).first;
  parameterizer.parameterize(mesh, bhd, cmap, uvmap, vimap);
  
  const size_t nvertices = sm.number_of_vertices();
  Rcpp::NumericMatrix out(nvertices, 2);
  int i = 0;
  for(auto v : sm.vertices()) {
    Point2 pt = uvmap[sm.halfedge(v)];
    out(i, 0) = pt.x();
    out(i, 1) = pt.y();
    i++;
  }
  return out;
}
