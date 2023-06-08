void writeFile(
    Rcpp::String filename, const int precision, 
    const bool binary, std::string comments,
    Rcpp::Nullable<Rcpp::NumericMatrix> normals_,
    Rcpp::Nullable<Rcpp::IntegerMatrix> fcolors_,
    Rcpp::Nullable<Rcpp::IntegerMatrix> vcolors_
) {
  EMesh3 meshcopy;
  MapBetweenFaces f2fmap_;
  boost::associative_property_map<MapBetweenFaces> f2fmap(f2fmap_);
  MapBetweenVertices v2vmap_;
  boost::associative_property_map<MapBetweenVertices> v2vmap(v2vmap_);
  CGAL::copy_face_graph(
    mesh, meshcopy, 
    CGAL::parameters::face_to_face_map(f2fmap).vertex_to_vertex_map(v2vmap)
  );
  if(normals_.isNotNull()) {
    Rcpp::NumericMatrix normals(normals_);
    CGALnormals_map vnormal = 
      meshcopy.add_property_map<vertex_descriptor, EVector3>(
          "v:normal", CGAL::NULL_VECTOR
      ).first;
    for(const auto& [vsource, vtarget] : v2vmap_) {
      Rcpp::NumericVector normal = normals(Rcpp::_, int(vsource));
      if(!Rcpp::NumericVector::is_na(normal(0))) {
        vnormal[vtarget] = EVector3(normal(0), normal(1), normal(2));
      }
    }
  }
  if(fcolors_.isNotNull()) {
    Rcpp::IntegerMatrix fcolors(fcolors_);
    EMesh3::Property_map<face_descriptor, Color> 
      fcolor = meshcopy.add_property_map<face_descriptor, Color>(
        "f:color", Color()
      ).first;      
    for(const auto& [fsource, ftarget] : f2fmap_) {
      Rcpp::IntegerVector color = fcolors(Rcpp::_, int(fsource));
      unsigned char red   = color(0);
      unsigned char green = color(1);
      unsigned char blue  = color(2);
      fcolor[ftarget] = Color(red, green, blue);
    }
  }
  if(vcolors_.isNotNull()) {
    Rcpp::IntegerMatrix vcolors(vcolors_);
    EMesh3::Property_map<vertex_descriptor, Color> 
      vcolor = meshcopy.add_property_map<vertex_descriptor, Color>(
        "v:color", Color()
      ).first;      
    for(const auto& [vsource, vtarget] : v2vmap_) {
      Rcpp::IntegerVector color = vcolors(Rcpp::_, int(vsource));
      unsigned char red   = color(0);
      unsigned char green = color(1);
      unsigned char blue  = color(2);
      vcolor[vtarget] = Color(red, green, blue);
    }
  }
  writeMeshFile(filename, precision, binary, comments, meshcopy);
}
