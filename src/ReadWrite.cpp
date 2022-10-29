#ifndef _HEADER_
#include "cgalMesh.h"
#endif

std::string toLower(std::string s) {
  for(char& c : s) {
    c = std::tolower(c);
  }
  return s;
}

EMesh3 readMeshFile(const std::string filename, const bool binary) {
  const std::string ext = toLower(filename.substr(filename.length() - 3, 3));
  EMesh3 mesh;
  bool ok = false;
  std::ifstream infile;
  if(binary) {
    infile.open(filename, std::ios::binary);
  } else {
    infile.open(filename);
  }
  if(ext == "ply") {
    ok = CGAL::IO::read_PLY(
      infile, mesh
    );
  } else if(ext == "off") {
    ok = CGAL::IO::read_OFF(infile, mesh);
  } else {
    infile.close();
    Rcpp::stop("Unknown file extension.");
  }
  infile.close();
  if(!ok) {
    Rcpp::stop("Reading failure.");
  }
  return mesh;
}

void writeMeshFile(const std::string filename,
                   const bool binary,
                   const int precision,
                   EMesh3 mesh) {
  const std::string ext = toLower(filename.substr(filename.length() - 3, 3));
  bool ok = false;
  std::ofstream outfile; 
  if(binary) {
    outfile.open(filename, std::ios::binary);
    CGAL::IO::set_binary_mode(outfile);
  } else {
    outfile.open(filename);
  }
  if(ext == "ply") {
    ok = CGAL::IO::write_PLY(
      outfile, mesh,
      CGAL::parameters::stream_precision(precision)
    );
  } else if(ext == "off") {
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
