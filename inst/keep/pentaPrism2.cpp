#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef EK::Point_3                                       Point;
typedef CGAL::Surface_mesh<Point>                         EMesh;
typedef boost::graph_traits<EMesh>::face_descriptor       face_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

int main() {
  // vertices -----------------------------------------------------------------
  const double rho = sqrt((5 - sqrt(5))/10);
  const double R = sqrt((25 - 11*sqrt(5))/10);
  const double pi = M_PI;
  std::vector<Point> vertices;
  for(int i = 0; i < 5; i++) {
    vertices.push_back(Point(rho*cos(2*i*pi/5), rho*sin(2*i*pi/5), 0.1))
  }
  for(int i = 0; i < 5; i++) {
    vertices.push_back(Point(rho*cos(2*i*pi/5), rho*sin(2*i*pi/5), -0.1))
  }
  for(int i = 0; i < 5; i++) {
    vertices.push_back(Point(R*cos(2*i*pi/5 + pi/5), R*sin(2*i*pi/5 + pi/5), 0.1))
  }
  for(int i = 0; i < 5; i++) {
    vertices.push_back(Point(R*cos(2*i*pi/5 + pi/5), R*sin(2*i*pi/5 + pi/5), -0.1))
  }
  // faces --------------------------------------------------------------------
  std::vector<std::vector<int>> triangles = 
    {
      {14, 0, 10}, 
      {10, 1, 11}, 
      {11, 2, 12},
      {12, 3, 13},
      {13, 4, 14},
      {15, 5, 19},
      {16, 6, 15},
      {17, 7, 16},
      {18, 8, 17},
      {19, 9, 18}
    };
  std::vector<std::vector<int>> rectangles =
  {
    {0, 5, 15, 10},
    {1, 6, 16, 11},
    {2, 7, 17, 12},
    {3, 8, 18, 13},
    {4, 9, 19, 14},
    {6, 1, 10, 15},
    {7, 2, 11, 16},
    {8, 3, 12, 17},
    {9, 4, 13, 18},
    {5, 0, 14, 19}
  };
  std::vector<std::vector<int>> pentagons =
  {
    {10, 11, 12, 13, 14},
    {19, 18, 17, 16, 15}
  };
  // mesh ---------------------------------------------------------------------
  EMesh mesh;
  const int nverts = vertices.size();
  for(int i = 0; i < nverts; i++) {
    mesh.add_vertex(vertices[i]);
  }
  const int ntriangles = triangles.size();
  for(int i = 0; i < ntriangles; i++) {
    std::vector<int> intface = triangles[i];
    std::vector<EMesh::Vertex_index> face;
    face.reserve(3);
    for(int k = 0; k < 3; k++) {
      face.emplace_back(CGAL::SM_Vertex_index(intface[k]));
    }
    face_descriptor fd = mesh.add_face(face);
    if(fd == EMesh::null_face()) {
      std::cout << "Cannot add triangle " + std::to_string(i) + ".";
      return 1;
    }
  }
  const int nrectangles = rectangles.size();
  for(int i = 0; i < nrectangles; i++) {
    std::vector<int> intface = rectangles[i];
    std::vector<EMesh::Vertex_index> face;
    face.reserve(4);
    for(int k = 0; k < 4; k++) {
      face.emplace_back(CGAL::SM_Vertex_index(intface[k]));
    }
    face_descriptor fd = mesh.add_face(face);
    if(fd == EMesh::null_face()) {
      std::cout << "Cannot add rectangle " + std::to_string(i) + ".";
      return 1;
    }
  }
  const int npentagons = pentagons.size();
  for(int i = 0; i < npentagons; i++) {
    std::vector<int> intface = pentagons[i];
    std::vector<EMesh::Vertex_index> face;
    face.reserve(5);
    for(int k = 0; k < 5; k++) {
      face.emplace_back(CGAL::SM_Vertex_index(intface[k]));
    }
    face_descriptor fd = mesh.add_face(face);
    if(fd == EMesh::null_face()) {
      std::cout << "Cannot add pentagon " + std::to_string(i) + ".";
      return 1;
    }
  }
  // triangulation
  const bool success = PMP::triangulate_faces(mesh);
  if(!success) {
    std::cout << "Triangulation has failed.";
    return 1;
  }
  return 0;
}