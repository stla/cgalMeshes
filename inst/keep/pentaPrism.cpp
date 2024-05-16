#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Surface_mesh<EK::Point_3> EMesh;
namespace PMP = CGAL::Polygon_mesh_processing;

int main() {
  EMesh mesh;
  const bool readoff = CGAL::IO::read_OFF("pentagrammicPrism.off", mesh);
  if(!readoff) {
    std::cout << "read_OFF has failed.";
    return 1;
  }
  const bool success = PMP::triangulate_faces(mesh);
  if(!success) {
    std::cout << "Triangulation has failed.";
    return 1;
  }
  return 0;
}