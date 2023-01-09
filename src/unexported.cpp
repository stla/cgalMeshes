#ifndef _HEADER_
#include "cgalMesh.h"
#endif

template <typename MeshT, typename PointT>
MeshT csoup2mesh(std::vector<PointT> points,
                 std::vector<std::vector<int>> faces,
                 const bool clean) {
  if(clean) {
    PMP::repair_polygon_soup(points, faces);
  }
  const bool success = PMP::orient_polygon_soup(points, faces);
  if(!success) {
    Rcpp::warning("Polygon orientation failed.");
  }
  MeshT mesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
  const bool valid = mesh.is_valid(false);
  if(!valid) {
    Rcpp::warning("The mesh is not valid.");
  }
  return mesh;
}

template EMesh3 csoup2mesh<EMesh3, EPoint3>(
  std::vector<EPoint3>, std::vector<std::vector<int>>, const bool
);

EMesh3 vf2mesh(const Rcpp::NumericMatrix vertices,
               const Rcpp::List faces) {
  EMesh3 mesh;
  const int nv = vertices.ncol();
  for(int j = 0; j < nv ; j++) {
    Rcpp::NumericVector vertex = vertices(Rcpp::_, j);
    EPoint3 pt = EPoint3(vertex(0), vertex(1), vertex(2));
    mesh.add_vertex(pt);
  }
  const int nf = faces.size();
  Rcpp::Rcout << "nf: " << nf << "\n";
  for(int i = 0; i < nf; i++) {
    Rcpp::IntegerVector intface = Rcpp::as<Rcpp::IntegerVector>(faces(i));
    const int sf = intface.size();
    Rcpp::Rcout << "sf: " << sf << "\n";
    std::vector<EMesh3::Vertex_index> face;
    face.reserve(sf);
    for(int k = 0; k < sf; k++) {
      Rcpp::Rcout << intface(k) << " - ";
      face.emplace_back(CGAL::SM_Vertex_index(intface(k)));
    }
    Rcpp::Rcout << "\n";
    face_descriptor fd = mesh.add_face(face);
    Rcpp::Rcout << "added face: " << fd << "\n";
  }
  Rcpp::Rcout << "nfaces: " << mesh.number_of_faces() << "\n";
  return mesh;
}

Rcpp::NumericVector defaultNormal() {
  Rcpp::NumericVector def = 
    {
      Rcpp::NumericVector::get_na(), 
      Rcpp::NumericVector::get_na(), 
      Rcpp::NumericVector::get_na()
    };
  return def;
}

EMesh3 makeMesh(const Rcpp::NumericMatrix vertices,
                const Rcpp::List faces,
                bool soup,
                const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_,
                const Rcpp::Nullable<Rcpp::StringVector> &vcolors_,
                const Rcpp::Nullable<Rcpp::StringVector> &fcolors_) {
  if(soup) {
    return csoup2mesh<EMesh3, EPoint3>(
      matrix_to_points3<EPoint3>(vertices), 
      list_to_faces(faces), 
      true
    );
  }
  EMesh3 mesh = vf2mesh(vertices, faces);
  if(normals_.isNotNull()) {
    Rcpp::NumericMatrix normals(normals_);
    if(mesh.number_of_vertices() != normals.ncol()) {
      Rcpp::stop(
        "The number of normals does not match the number of vertices."
      );
    }
    Rcpp::NumericVector def = defaultNormal();
    Normals_map normalsmap = 
      mesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
        "v:normal", def
      ).first;
    for(int j = 0; j < normals.ncol(); j++) {
      Rcpp::NumericVector normal = normals(Rcpp::_, j);
      normalsmap[CGAL::SM_Vertex_index(j)] = normal;
    }
  }
  if(vcolors_.isNotNull()) {
    Rcpp::StringVector vcolors(vcolors_);
    if(mesh.number_of_vertices() != vcolors.size()) {
      Rcpp::stop(
        "The number of vertex colors does not match the number of vertices."
      );
    }
    Vcolors_map vcolorsmap = 
      mesh.add_property_map<vertex_descriptor, std::string>("v:color", "").first;
    for(int i = 0; i < vcolors.size(); i++) {
      vcolorsmap[CGAL::SM_Vertex_index(i)] = vcolors(i);
    }
  }
  if(fcolors_.isNotNull()) {
    Rcpp::StringVector fcolors(fcolors_);
    if(mesh.number_of_faces() != fcolors.size()) {
      Rcpp::stop(
        "The number of face colors does not match the number of vertices."
      );
    }
    Fcolors_map fcolorsmap = 
      mesh.add_property_map<face_descriptor, std::string>("f:color", "").first;
    for(int i = 0; i < fcolors.size(); i++) {
      fcolorsmap[CGAL::SM_Face_index(i)] = fcolors(i);
    }
  }
  return mesh;
}

EMesh3 cloneMesh(
  EMesh3& mesh, std::vector<std::string> props
) {
  EMesh3::Property_map<vertex_descriptor, vertex_descriptor> v2vmap = 
    mesh.add_property_map<vertex_descriptor, vertex_descriptor>("v:v").first;
  EMesh3::Property_map<face_descriptor, face_descriptor> f2fmap = 
    mesh.add_property_map<face_descriptor, face_descriptor>("f:f").first;
  EMesh3 out;
  CGAL::copy_face_graph(
    mesh, out, 
    CGAL::parameters::vertex_to_vertex_map(v2vmap).face_to_face_map(f2fmap)
  );
  for(int i = 0; i < props.size(); i++) {
    std::string prop = props[i];
    if(prop == "f:color") {
      std::pair<Fcolors_map, bool> pmap_ = 
        mesh.property_map<face_descriptor, std::string>(prop);
      if(pmap_.second) {
        Fcolors_map pmap = 
          out.add_property_map<face_descriptor, std::string>(
            prop, ""
          ).first;
        for(EMesh3::Face_index fi : mesh.faces()) {
          pmap[f2fmap[fi]] = pmap_.first[fi];
        }
      }
    } else if(prop == "v:color") {
      std::pair<Vcolors_map, bool> pmap_ = 
        mesh.property_map<vertex_descriptor, std::string>(prop);
      if(pmap_.second) {
        Vcolors_map pmap = 
          out.add_property_map<vertex_descriptor, std::string>(
            prop, ""
          ).first;
        for(EMesh3::Vertex_index vi : mesh.vertices()) {
          pmap[v2vmap[vi]] = pmap_.first[vi];
        }
      }
    } else if(prop == "v:normal") {
      std::pair<Normals_map, bool> pmap_ = 
        mesh.property_map<vertex_descriptor, Rcpp::NumericVector>(prop);
      if(pmap_.second) {
        Normals_map pmap = 
          out.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
            prop, defaultNormal()
          ).first;
        for(EMesh3::Vertex_index vi : mesh.vertices()) {
          pmap[v2vmap[vi]] = pmap_.first[vi];
        }
      }
    } else if(prop == "v:scalar") {
      std::pair<Vscalars_map, bool> pmap_ = 
        mesh.property_map<vertex_descriptor, double>(prop);
      if(pmap_.second) {
        Vscalars_map pmap = 
          out.add_property_map<vertex_descriptor, double>(
            prop, nan("")
          ).first;
        for(EMesh3::Vertex_index vi : mesh.vertices()) {
          pmap[v2vmap[vi]] = pmap_.first[vi];
        }
      }
    } else if(prop == "f:scalar") {
      std::pair<Fscalars_map, bool> pmap_ = 
        mesh.property_map<face_descriptor, double>(prop);
      if(pmap_.second) {
        Fscalars_map pmap = 
          out.add_property_map<face_descriptor, double>(
            prop, nan("")
          ).first;
        for(EMesh3::Face_index fi : mesh.faces()) {
          pmap[f2fmap[fi]] = pmap_.first[fi];
        }
      }
    }
  }
  mesh.remove_property_map(v2vmap);
  mesh.remove_property_map(f2fmap);
  return out;
}

void removeProperties(
  EMesh3& mesh, std::vector<std::string> props
) {
  for(int i = 0; i < props.size(); i++) {
    std::string prop = props[i];
    if(prop == "f:color") {
      std::pair<Fcolors_map, bool> pmap_ = 
        mesh.property_map<face_descriptor, std::string>("f:color");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "v:color") {
      std::pair<Vcolors_map, bool> pmap_ = 
        mesh.property_map<vertex_descriptor, std::string>("v:color");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "v:normal") {
      std::pair<Normals_map, bool> pmap_ = 
        mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "v:scalar") {
      std::pair<Vscalars_map, bool> pmap_ = 
        mesh.property_map<vertex_descriptor, double>("v:scalar");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    } else if(prop == "f:scalar") {
      std::pair<Fscalars_map, bool> pmap_ = 
        mesh.property_map<face_descriptor, double>("f:scalar");
      if(pmap_.second) {
        mesh.remove_property_map(pmap_.first);
      }
    }
  }
}

/* MaybeFcolorMap copy_fcolor(EMesh3& mesh) {
  std::pair<Fcolors_map, bool> fcolors_ = 
    mesh.property_map<face_descriptor, std::string>("f:color");
  bool has_fcolor = fcolors_.second;
  std::map<face_descriptor, std::string> fcolorsmap;
  if(has_fcolor) {
    for(EMesh3::Face_index fi: mesh.faces()) {
      fcolorsmap[fi] = fcolors_.first[fi];
    }
    mesh.remove_property_map(fcolors_.first);
  }
  return std::make_pair(fcolorsmap, has_fcolor);
}

MaybeVcolorMap copy_vcolor(EMesh3& mesh) {
  std::pair<Vcolors_map, bool> vcolors_ = 
    mesh.property_map<vertex_descriptor, std::string>("v:color");
  bool has_vcolor = vcolors_.second;
  std::map<vertex_descriptor, std::string> vcolorsmap;
  if(has_vcolor) {
    for(EMesh3::Vertex_index vi: mesh.vertices()) {
      vcolorsmap[vi] = vcolors_.first[vi];
    }
    mesh.remove_property_map(vcolors_.first);
  }
  return std::make_pair(vcolorsmap, has_vcolor);
}

MaybeNormalMap copy_vnormal(EMesh3& mesh) {
  std::pair<Normals_map, bool> vnormals_ = 
    mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
  bool has_vnormal = vnormals_.second;
  std::map<vertex_descriptor, Rcpp::NumericVector> vnormalsmap;
  if(has_vnormal) {
    for(EMesh3::Vertex_index vi: mesh.vertices()) {
      vnormalsmap[vi] = vnormals_.first[vi];
    }
    mesh.remove_property_map(vnormals_.first);
  }
  return std::make_pair(vnormalsmap, has_vnormal);
}
 */
template <typename Keytype, typename Valuetype>
std::pair<std::map<Keytype, Valuetype>, bool> copy_prop(
  EMesh3& mesh, std::string propname
) {
  std::pair<EMesh3::Property_map<Keytype, Valuetype>, bool> pmap_ = 
    mesh.property_map<Keytype, Valuetype>(propname);
  bool has_prop = pmap_.second;
  std::map<Keytype, Valuetype> pmap;
  if(has_prop) {
    std::string descriptor = propname.substr(0, 1);
    std::size_t n = 
      descriptor == "v" ? mesh.number_of_vertices() : mesh.number_of_faces();
    for(std::size_t idx = 0; idx < n; idx++) {
      pmap[Keytype(idx)] = pmap_.first[Keytype(idx)];
    }
    mesh.remove_property_map(pmap_.first);
  }
  return std::make_pair(pmap, has_prop);
}

template MaybeFcolorMap copy_prop<face_descriptor, std::string>(EMesh3&, std::string);
template MaybeFscalarMap copy_prop<face_descriptor, double>(EMesh3&, std::string);
template MaybeVcolorMap copy_prop<vertex_descriptor, std::string>(EMesh3&, std::string);
template MaybeVscalarMap copy_prop<vertex_descriptor, double>(EMesh3&, std::string);
template MaybeNormalMap copy_prop<vertex_descriptor, Rcpp::NumericVector>(EMesh3&, std::string);

template <typename SourceDescriptor, typename TargetDescriptor, typename Valuetype>
void copy_property(
  EMesh3& mesh, EMesh3& fmesh, std::map<SourceDescriptor, TargetDescriptor> dmap, std::string propname
) {
  std::pair<EMesh3::Property_map<TargetDescriptor, Valuetype>, bool> pmap_ = 
    mesh.property_map<TargetDescriptor, Valuetype>(propname);
  bool has_prop = pmap_.second;
  if(has_prop) {
    EMesh3::Property_map<TargetDescriptor, Valuetype> pmap = 
      fmesh.add_property_map<TargetDescriptor, Valuetype>(
        propname
      ).first;
    for(const auto& [source_decriptor, target_decriptor] : dmap) {
      pmap[target_decriptor] = pmap_.first[TargetDescriptor(int(source_decriptor))];
    }
  }
}

template void copy_property<ffg_vertex_descriptor, vertex_descriptor, Rcpp::NumericVector>(EMesh3&, EMesh3&, MapBetweenVertexDescriptors, std::string);
template void copy_property<ffg_vertex_descriptor, vertex_descriptor, std::string>(EMesh3&, EMesh3&, MapBetweenVertexDescriptors, std::string);
template void copy_property<ffg_vertex_descriptor, vertex_descriptor, double>(EMesh3&, EMesh3&, MapBetweenVertexDescriptors, std::string);
template void copy_property<ffg_face_descriptor, face_descriptor, std::string>(EMesh3&, EMesh3&, MapBetweenFaceDescriptors, std::string);
template void copy_property<ffg_face_descriptor, face_descriptor, double>(EMesh3&, EMesh3&, MapBetweenFaceDescriptors, std::string);


void triangulateMesh(EMesh3& mesh) {
  MaybeFcolorMap fcolormap_ = 
    copy_prop<face_descriptor, std::string>(mesh, "f:color");
  const bool hasFcolors = fcolormap_.second;
  MaybeFscalarMap fscalarmap_ = 
    copy_prop<face_descriptor, double>(mesh, "f:scalar");
  const bool hasFscalars = fscalarmap_.second;
  removeProperties(mesh, {"v:normal"});
  TriangulateVisitor vis;
  const bool success = 
    PMP::triangulate_faces(mesh, CGAL::parameters::visitor(vis));
  if(!success) {
    Rcpp::stop("Triangulation has failed.");
  }
  if(hasFcolors) {
    MapBetweenFaces fmap = *(vis.fmap);
    Fcolors_map fcolors = 
      mesh.add_property_map<face_descriptor, std::string>(
        "f:color", ""
      ).first;
    for(EMesh3::Face_index fi: mesh.faces()) {
      fcolors[fi] = fcolormap_.first[fmap[fi]];
    }
  }
  if(hasFscalars) {
    MapBetweenFaces fmap = *(vis.fmap);
    Fscalars_map fscalars = 
      mesh.add_property_map<face_descriptor, double>(
        "f:scalar", nan("")
      ).first;
    for(EMesh3::Face_index fi: mesh.faces()) {
      fscalars[fi] = fscalarmap_.first[fmap[fi]];
    }
  }
}

Rcpp::List clipping(EMesh3& tm, EMesh3& clipper, const bool clipVolume) {
  if(!CGAL::is_triangle_mesh(tm)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  if(clipVolume) {
    if(PMP::does_self_intersect(tm)) {
      Rcpp::stop("The mesh self-intersects.");
    }
  }
  if(!CGAL::is_triangle_mesh(clipper)) {
    Rcpp::stop("The clipping mesh is not triangle.");
  }
  if(!CGAL::is_closed(clipper)) {
    Rcpp::stop("The clipping mesh is not closed.");
  }
  if(PMP::does_self_intersect(clipper)) {
    Rcpp::stop("The clipping mesh self-intersects.");
  }
    
  MaybeFcolorMap fcolorMap_ = 
    copy_prop<face_descriptor, std::string>(tm, "f:color");
  MaybeFcolorMap fcolorMap2_ = 
    copy_prop<face_descriptor, std::string>(clipper, "f:color");
  const bool hasColors = fcolorMap_.second && fcolorMap2_.second;
  MaybeFscalarMap fscalarMap_ = 
    copy_prop<face_descriptor, double>(tm, "f:scalar");
  MaybeFscalarMap fscalarMap2_ = 
    copy_prop<face_descriptor, double>(clipper, "f:scalar");
  const bool hasScalars = fscalarMap_.second && fscalarMap2_.second;

  std::size_t nfaces = tm.number_of_faces();
  std::size_t nfaces_clipper = clipper.number_of_faces();
  
  std::size_t undetermined = 999999;

  ClipVisitor vis;
  Face_index_map fimap = 
    tm.add_property_map<face_descriptor, std::size_t>("f:i", undetermined).first;
  const bool doNotModify = !clipVolume;
  const bool clipping = PMP::clip(
    tm, clipper,
    PMP::parameters::clip_volume(clipVolume)
      .visitor(vis).face_index_map(fimap),
    PMP::parameters::clip_volume(clipVolume).do_not_modify(doNotModify)
  );
  if(!clipping) {
    Rcpp::stop("Clipping has failed.");
  }

  tm.collect_garbage();

  if(!clipVolume) {
    MapBetweenFaces fmap = *(vis.fmap_tm);
    if(hasColors || hasScalars) {
      std::map<face_descriptor, std::string> fcolorMap;
      Fcolors_map newfcolor;
      std::map<face_descriptor, double> fscalarMap;
      Fscalars_map newfscalar;
      if(hasColors) {
        fcolorMap = fcolorMap_.first;
        newfcolor = 
          tm.add_property_map<face_descriptor, std::string>(
            "f:color", ""
          ).first;
      }
      if(hasScalars) {
        fscalarMap = fscalarMap_.first;
        newfscalar = 
          tm.add_property_map<face_descriptor, double>(
            "f:scalar", nan("")
          ).first;
      }
      for(EMesh3::Face_index fi : tm.faces()) {
        face_descriptor fd = CGAL::SM_Face_index(fimap[fi]);
        face_descriptor fdnew = size_t(fd) < nfaces ? fd : fmap[fd];
        if(hasColors) {
          newfcolor[fi] = fcolorMap[fdnew];
        }
        if(hasScalars) {
          newfscalar[fi] = fscalarMap[fdnew];
        }
      }
    }
    tm.remove_property_map(fimap);
    return Rcpp::List::create();
  }

  MapBetweenFaces fmap_tm = *(vis.fmap_tm);
  MapBetweenFaces fmap_clipper = *(vis.fmap_clipper);
  MapBetweenFaces ftargets = *(vis.ftargets);
  MapBetweenFaces zeros;
  std::map<face_descriptor, std::string> fcolorMap;
  std::map<face_descriptor, std::string> fcolorMap2;
  std::map<face_descriptor, std::string> fcolorMap_clipper;
  std::map<face_descriptor, double> fscalarMap;
  std::map<face_descriptor, double> fscalarMap2;
  std::map<face_descriptor, double> fscalarMap_clipper;
  if(hasColors) {
    fcolorMap = fcolorMap_.first;
    fcolorMap2 = fcolorMap2_.first;
  }
  if(hasScalars) {
    fscalarMap = fscalarMap_.first;
    fscalarMap2 = fscalarMap2_.first;
  }
  {
    int fdi = 0;
    for(auto it = ftargets.rbegin(); it != ftargets.rend(); ++it) {
      face_descriptor fd = it->second;
      if(std::size_t(fd) >= nfaces_clipper) {
        fd = fmap_clipper[fd];
      }
      if(hasColors) {
        fcolorMap_clipper[it->first] = fcolorMap2[fd];
      }
      if(hasScalars) {
        fscalarMap_clipper[it->first] = fscalarMap2[fd];
      }
      while(fimap[CGAL::SM_Face_index(fdi)] != undetermined) {
        fdi++;
      }
      zeros[CGAL::SM_Face_index(fdi++)] = it->first;
    }
  }

  Face_index_map whichPart = 
    tm.add_property_map<face_descriptor, std::size_t>("f:which").first;
  Fcolors_map newfcolor;
  Fscalars_map newfscalar;
  if(hasColors) {
    newfcolor = 
      tm.add_property_map<face_descriptor, std::string>("f:color", "").first;
  }
  if(hasScalars) {
    newfscalar = 
      tm.add_property_map<face_descriptor, double>("f:scalar", nan("")).first;
  }

  for(EMesh3::Face_index fi : tm.faces()) {
    std::size_t ifi = fimap[fi];
    face_descriptor fd = CGAL::SM_Face_index(ifi);
    if(fi != fd && ifi != undetermined) {
      whichPart[fi] = 0;
      if(hasColors) {
        newfcolor[fi] = fcolorMap[fmap_tm[fd]];
      }
      if(hasScalars) {
        newfscalar[fi] = fscalarMap[fmap_tm[fd]];
      }
    } else if(ifi != undetermined) {
      whichPart[fi] = 0;
      if(hasColors) {
        newfcolor[fi] = 
          ifi < nfaces ? fcolorMap[fd] : fcolorMap[fmap_tm[fd]];
      }
      if(hasScalars) {
        newfscalar[fi] = 
          ifi < nfaces ? fscalarMap[fd] : fscalarMap[fmap_tm[fd]];
      }
    } else if(ifi == undetermined) {
      whichPart[fi] = 1;
      if(hasColors) {
        newfcolor[fi] = fcolorMap_clipper[zeros[fi]];
      }
      if(hasScalars) {
        newfscalar[fi] = fscalarMap_clipper[zeros[fi]];
      }
    } 
  }

  tm.remove_property_map(fimap);

  EMesh3 tmesh;
  {
    Filtered_graph ffg(tm, 0, whichPart);
    MapBetweenFaceDescriptors f2fmap_;
    boost::associative_property_map<MapBetweenFaceDescriptors> f2fmap(f2fmap_);
    CGAL::copy_face_graph(
      ffg, tmesh, CGAL::parameters::face_to_face_map(f2fmap)
    );
    copy_property<ffg_face_descriptor, face_descriptor, std::string>(
      tm, tmesh, f2fmap_, "f:color"
    );
    copy_property<ffg_face_descriptor, face_descriptor, double>(
      tm, tmesh, f2fmap_, "f:scalar"
    );
  }

  EMesh3 tmesh2;
  {
    Filtered_graph ffg(tm, 1, whichPart);
    MapBetweenFaceDescriptors f2fmap_;
    boost::associative_property_map<MapBetweenFaceDescriptors> f2fmap(f2fmap_);
    CGAL::copy_face_graph(
      ffg, tmesh2, CGAL::parameters::face_to_face_map(f2fmap)
    );
    copy_property<ffg_face_descriptor, face_descriptor, std::string>(
      tm, tmesh2, f2fmap_, "f:color"
    );
    copy_property<ffg_face_descriptor, face_descriptor, double>(
      tm, tmesh2, f2fmap_, "f:scalar"
    );
  }

  Rcpp::List meshes = Rcpp::List::create(
    Rcpp::Named("mesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(tmesh), false),
    Rcpp::Named("mesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(tmesh2), false)
  );

  return meshes;
}

void clippingToPlane(EMesh3& tm, EPlane3 plane, const bool clipVolume) {
  if(!CGAL::is_triangle_mesh(tm)) {
    Rcpp::stop("The mesh is not triangle.");
  }
  if(PMP::does_self_intersect(tm)) {
    Rcpp::stop("The mesh self-intersects.");
  }

  std::size_t nfaces = tm.number_of_faces();

  MaybeFcolorMap fcolorMap_ = 
    copy_prop<face_descriptor, std::string>(tm, "f:color");
  const bool hasColors = fcolorMap_.second;
  MaybeFscalarMap fscalarMap_ = 
    copy_prop<face_descriptor, double>(tm, "f:scalar");
  const bool hasScalars = fscalarMap_.second;

  ClipVisitor vis;
  const bool clipping = PMP::clip(
    tm, plane,
    PMP::parameters::clip_volume(clipVolume).visitor(vis)
  );
  if(!clipping) {
    Rcpp::stop("Clipping has failed.");
  }

  Face_index_map fimap = 
    tm.add_property_map<face_descriptor, std::size_t>("f:i", 9999999).first;
  for(face_descriptor fd : tm.faces()) {
    fimap[fd] = std::size_t(fd);
  }

  tm.collect_garbage();

  if(!clipVolume){
    MapBetweenFaces fmap = *(vis.fmap_tm);

    if(hasColors || hasScalars) {
      std::map<face_descriptor, std::string> fcolorMap;
      Fcolors_map newfcolor;
      std::map<face_descriptor, double> fscalarMap;
      Fscalars_map newfscalar;
      if(hasColors) {
        fcolorMap = fcolorMap_.first;
        newfcolor = 
          tm.add_property_map<face_descriptor, std::string>(
            "f:color", ""
          ).first;
      }
      if(hasScalars) {
        fscalarMap = fscalarMap_.first;
        newfscalar = 
          tm.add_property_map<face_descriptor, double>(
            "f:scalar", nan("")
          ).first;
      }
      for(EMesh3::Face_index fi : tm.faces()) {
        std::size_t ffi = fimap[fi];
        face_descriptor fd = CGAL::SM_Face_index(ffi);
        fd = ffi < nfaces ? fd : fmap[fd];
        if(hasColors) {
          newfcolor[fi] = fcolorMap[fd];
        }
        if(hasScalars) {
          newfscalar[fi] = fscalarMap[fd];
        }
      }
    }
    tm.remove_property_map(fimap);
  }

  MapBetweenFaces fmap = *(vis.fmap_tm);
  MapBetweenFaces ftargets = *(vis.ftargets);

  if(hasColors || hasScalars) {
    std::map<face_descriptor, std::string> fcolorMap;
    Fcolors_map newfcolor;
    std::map<face_descriptor, double> fscalarMap;
    Fscalars_map newfscalar;
    if(hasColors) {
      fcolorMap = fcolorMap_.first;
      newfcolor = 
        tm.add_property_map<face_descriptor, std::string>(
          "f:color", ""
        ).first;
    }
    if(hasScalars) {
      fscalarMap = fscalarMap_.first;
      newfscalar = 
        tm.add_property_map<face_descriptor, double>(
          "f:scalar", nan("")
        ).first;
    }
    for(EMesh3::Face_index fi : tm.faces()) {
      std::size_t ffi = fimap[fi];
      face_descriptor fd = CGAL::SM_Face_index(ffi);
      if(auto search = ftargets.find(fd); search != ftargets.end()) {
        if(hasColors) {
          newfcolor[fi] = "black";
        }
        if(hasScalars) {
          newfscalar[fi] = nan("");
        }
        ftargets.erase(fd);
      } else {
        fd = ffi < nfaces ? fd : fmap[fd];
        if(hasColors) {
          newfcolor[fi] = fcolorMap[fd];
        }
        if(hasScalars) {
          newfscalar[fi] = fscalarMap[fd];
        }
      }
    }
  }
  tm.remove_property_map(fimap);

}