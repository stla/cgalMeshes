#ifndef _HEADER_
#include "cgalMesh.h"
#endif

std::string toLower(std::string s) {
  for(char& c : s) {
    c = std::tolower(c);
  }
  return s;
}

EMesh3 readMeshFile(const std::string filename) {
  EMesh3 mesh;
  const bool ok = PMP::IO::read_polygon_mesh(
    filename, mesh, CGAL::parameters::verbose(true)
  );
  if(!ok) {
    Rcpp::stop("Reading failure.");
  }
  const bool valid = mesh.is_valid(false);
  if(!valid) {
    Rcpp::warning("The mesh is not valid.");
  }
  return mesh;
}

// void writeMeshFile(const std::string filename,
//                    const int precision,
//                    EMesh3 mesh) {
//   const bool ok = CGAL::IO::write_polygon_mesh(
//     filename, mesh,
//     CGAL::parameters::stream_precision(precision).verbose(true)
//   );
//   if(!ok) {
//     Rcpp::stop("Failed to write file.");
//   }
// }

void writeMeshFile(const std::string filename,
                   const int precision,
                   const bool binary,
                   EMesh3& mesh) {
  const std::string ext = toLower(filename.substr(filename.length() - 4, 4));
  bool ok = false;
  std::ofstream outfile;
  if(binary) {
    outfile.open(filename, std::ios::binary);
    CGAL::IO::set_binary_mode(outfile);
  } else {
    outfile.open(filename);
  }
  if(ext == ".ply") {
    ok = CGAL::IO::write_PLY(
      outfile, mesh,
      CGAL::parameters::stream_precision(precision)
    );
  } else if(ext == ".off") {
    ok = CGAL::IO::write_OFF(
      outfile, mesh,
      CGAL::parameters::stream_precision(precision)
    );
  } else {
    outfile.close();
    Rcpp::stop("Unknown file extension.");
  }
  outfile.close();
  if(!ok) {
    Rcpp::stop("Failed to write file.");
  }
}
