#ifndef _HEADER_
#include "cgalMesh.h"
#endif

// centroid of a face of the mesh
EPoint3 centroid(EMesh3& mesh, EMesh3::Face_index f) {
  EK::FT x = 0;
  EK::FT y = 0;
  EK::FT z = 0;
  int nvertices = 0;
  for(EMesh3::Vertex_index v :
           vertices_around_face(mesh.halfedge(f), mesh)) {
    nvertices++;
    EPoint3 p = mesh.point(v);
    x += p.x();
    y += p.y();
    z += p.z();
  }
  return EPoint3(x / nvertices, y / nvertices, z / nvertices);
}

EMesh3 dualMesh(EMesh3& mesh) {
  EMesh3 dualmesh;
  // property to remember new vertices per face
  auto fvertex = 
    mesh.add_property_map<EMesh3::Face_index, EMesh3::Vertex_index>(
      "f:vertex", EMesh3::null_vertex()).first;
  // for each face add its centroid to the dual mesh
  for(EMesh3::Face_index f : mesh.faces()) {
    fvertex[f] = dualmesh.add_vertex(centroid(mesh, f));
  }
  // add new face for each vertex
  for(EMesh3::Vertex_index v : mesh.vertices()) {
    std::vector<EMesh3::Vertex_index> vs;
    for(EMesh3::Face_index f : faces_around_target(mesh.halfedge(v), mesh)) {
      vs.push_back(fvertex[f]);
    }
    dualmesh.add_face(vs);
  }
  return dualmesh;
}
