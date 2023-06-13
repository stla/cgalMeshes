#ifndef _HEADER_
#include "cgalMesh.h"
#endif


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
std::string toLower(std::string s) {
  for(char& c : s) {
    c = std::tolower(c);
  }
  return s;
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
EMesh3 readValidMesh(const std::string filename, bool binary) {
  Mesh3 mesh;
  const std::string ext = toLower(filename.substr(filename.length() - 4, 4));
  bool ok = false;
  std::ifstream infile;
  if(binary) {
    infile.open(filename, std::ios::binary);
  } else {
    infile.open(filename);
  }
  std::string comments;
  if(ext == ".ply") {
    ok = CGAL::IO::read_PLY(infile, mesh, comments);
  } else if(ext == ".off") {
    ok = CGAL::IO::read_OFF(infile, mesh);
  } else {
    ok = PMP::IO::read_polygon_mesh(
      filename, mesh, CGAL::parameters::verbose(true)
    );
  }
  infile.close();
  if(!ok) {
    Rcpp::stop("Reading failure.");
  }
  if(!comments.empty()) {
    Message("Comments found in " + filename + ":");
    Message("------------------");
    Message(comments);
  }
  const bool valid = mesh.is_valid(false);
  if(!valid) {
    Rcpp::warning("The mesh is not valid.");
  }
  
  EMesh3 emesh = epick2epeck(mesh);

  std::pair<std::map<fdescr, Color>, bool> fcolor_ =
    copy_prop<fdescr, Color, K>(mesh, "f:color");
  if(fcolor_.second) {
    Fcolors_map fcolor_map =
      emesh.add_property_map<face_descriptor, std::string>(
                          "f:color", ""
                        ).first;
    for(Mesh3::Face_index fi : mesh.faces()) {
      const Color color = fcolor_.first[fi];
      const double red   = color.red();
      const double green = color.green();
      const double blue  = color.blue();
      fcolor_map[fi] = RcppColors::rgb2hex(red, green, blue);
    }
  }
  
  std::pair<std::map<vxdescr, Color>, bool> vcolor_ =
    copy_prop<vxdescr, Color, K>(mesh, "v:color");
  if(vcolor_.second) {
    Vcolors_map vcolor_map =
      emesh.add_property_map<vertex_descriptor, std::string>(
                          "v:color", ""
                        ).first;
    for(Mesh3::Vertex_index vi : mesh.vertices()) {
      const Color color = vcolor_.first[vi];
      const double red   = color.red();
      const double green = color.green();
      const double blue  = color.blue();
      vcolor_map[vi] = RcppColors::rgb2hex(red, green, blue);
    }
  }
  
  std::pair<std::map<vxdescr, Vector3>, bool> vnormal_ =
    copy_prop<vxdescr, Vector3, K>(mesh, "v:normal");
  if(vnormal_.second) {
    Normals_map vnormal_map =
      emesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
                          "v:normal", defaultNormal()
                        ).first;
    for(Mesh3::Vertex_index vi : mesh.vertices()) {
      Rcpp::NumericVector rcppnormal(3);
      const Vector3 normal = vnormal_.first[vi];
      rcppnormal(0) = normal.x();
      rcppnormal(1) = normal.y();
      rcppnormal(2) = normal.z();
      vnormal_map[vi] = rcppnormal;
    }
  }
  
  return emesh;
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
EMesh3 readPolygonSoup(const std::string filename, bool binary) {
  std::vector<EPoint3> points;
  std::vector<std::vector<size_t>> faces;
  const std::string ext = toLower(filename.substr(filename.length() - 4, 4));
  bool ok = false;
  std::ifstream infile;
  if(binary) {
    infile.open(filename, std::ios::binary);
  } else {
    infile.open(filename);
  }
  if(ext == ".ply") {
    ok = CGAL::IO::read_PLY(infile, points, faces);
  } else if(ext == ".off") {
    ok = CGAL::IO::read_OFF(infile, points, faces);
  } else {
    ok = CGAL::IO::read_polygon_soup(
      filename, points, faces, CGAL::parameters::verbose(true)
    );
  }
  infile.close();
  if(!ok) {
    Rcpp::stop("Reading failure.");
  }
  return csoup2mesh<EMesh3, EPoint3>(points, faces, true);
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
EMesh3 readMeshFile(const std::string filename, bool binary, bool soup) {
  if(soup) {
    return readPolygonSoup(filename, binary);
  } else {
    return readValidMesh(filename, binary);
  }
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
void writeMeshFile(const std::string filename,
                   const int precision,
                   const bool binary,
                   std::string comments,
                   Mesh3& mesh) {
  const std::string ext = toLower(filename.substr(filename.length() - 4, 4));
  bool ok = false;
  std::ofstream outfile;
  if(binary && ext == ".ply") {
    outfile.open(filename, std::ios::binary);
    CGAL::IO::set_binary_mode(outfile);
  } else {
    outfile.open(filename);
  }
  if(ext == ".ply") {
    ok = CGAL::IO::write_PLY(
      outfile, mesh, comments,
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
