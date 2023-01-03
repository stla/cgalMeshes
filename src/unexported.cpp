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
  for(int i = 0; i < nf; i++) {
    Rcpp::IntegerVector intface = Rcpp::as<Rcpp::IntegerVector>(faces(i));
    const int sf = intface.size();
    std::vector<EMesh3::Vertex_index> face;
    face.reserve(sf);
    for(int k = 0; k < sf; k++) {
      face.emplace_back(CGAL::SM_Vertex_index(intface(k)));
    }
    mesh.add_face(face);
  }
  return mesh;
}


EMesh3 makeMesh(const Rcpp::NumericMatrix vertices,
                const Rcpp::List faces,
                bool soup,
                const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_,
                const Rcpp::Nullable<Rcpp::StringVector> &vcolors_,
                const Rcpp::Nullable<Rcpp::StringVector> &fcolors_) {
  Rcpp::Rcout << "soup: " << soup << "\n";
  if(soup) {
    return csoup2mesh<EMesh3, EPoint3>(
      matrix_to_points3<EPoint3>(vertices), 
      list_to_faces(faces), 
      true
    );
  }
  Rcpp::Rcout << "soup is false\n";
  EMesh3 mesh = vf2mesh(vertices, faces);
  if(normals_.isNotNull()) {
    Rcpp::NumericMatrix normals(normals_);
    if(mesh.number_of_vertices() != normals.ncol()) {
      Rcpp::stop("pb normals");
    }
    Rcpp::NumericVector def = {Rcpp::NumericVector::get_na(), Rcpp::NumericVector::get_na(), Rcpp::NumericVector::get_na()};
    Normals_map normalsmap = 
      mesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal", def).first;
    for(int j = 0; j < normals.ncol(); j++) {
      Rcpp::NumericVector normal = normals(Rcpp::_, j);
      normalsmap[CGAL::SM_Vertex_index(j)] = normal;
    }
  }
  if(vcolors_.isNotNull()) {
    Rcpp::StringVector vcolors(vcolors_);
    if(mesh.number_of_vertices() != vcolors.size()) {
      Rcpp::stop("pb vcolors");
    }
    Vcolors_map vcolorsmap = 
      mesh.add_property_map<vertex_descriptor, std::string>("v:color", "").first;
    for(int i = 0; i < vcolors.size(); i++) {
      vcolorsmap[CGAL::SM_Vertex_index(i)] = vcolors(i);
    }
  }
  if(fcolors_.isNotNull()) {
    Rcpp::Rcout << "fcolors is not null\n";
    Rcpp::StringVector fcolors(fcolors_);
    if(mesh.number_of_faces() != fcolors.size()) {
      Rcpp::stop("pb fcolors");
    }
    Fcolors_map fcolorsmap = 
      mesh.add_property_map<face_descriptor, std::string>("f:color", "").first;
    for(int i = 0; i < fcolors.size(); i++) {
      fcolorsmap[CGAL::SM_Face_index(i)] = fcolors(i);
    }
  }
  return mesh;
}

EMesh3 cloneMesh(EMesh3& mesh, const bool vnormals, const bool vcolors, const bool fcolors) {
  EMesh3 out;
  CGAL::copy_face_graph(mesh, out);
  if(fcolors) {
    std::pair<Fcolors_map, bool> fcolors_ = 
      mesh.property_map<face_descriptor, std::string>("f:color");
    if(fcolors_.second) {
      Fcolors_map fcolorsmap = 
        out.add_property_map<face_descriptor, std::string>("f:color", "").first;
      for(EMesh3::Face_index fi : out.faces()) {
        fcolorsmap[fi] = fcolors_.first[fi];
      }
    }
  }
  if(vcolors) {
    std::pair<Vcolors_map, bool> vcolors_ = 
      mesh.property_map<vertex_descriptor, std::string>("v:color");
    if(vcolors_.second) {
      Vcolors_map vcolorsmap = 
        out.add_property_map<vertex_descriptor, std::string>("v:color", "").first;
      for(EMesh3::Vertex_index vi : out.vertices()) {
        vcolorsmap[vi] = vcolors_.first[vi];
      }
    }
  }
  if(vnormals) {
    std::pair<Normals_map, bool> vnormals_ = 
      mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
    if(vnormals_.second) {
      Rcpp::NumericVector def = {Rcpp::NumericVector::get_na(), Rcpp::NumericVector::get_na(), Rcpp::NumericVector::get_na()};
      Normals_map vnormalsmap = 
        out.add_property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal", def).first;
      for(EMesh3::Vertex_index vi : out.vertices()) {
        vnormalsmap[vi] = vnormals_.first[vi];
      }
    }
  }
  return out;
}

void removeProperties(EMesh3& mesh, const bool vnormals, const bool vcolors, const bool fcolors) {
  if(fcolors) {
    std::pair<Fcolors_map, bool> fcolors_ = 
      mesh.property_map<face_descriptor, std::string>("f:color");
    if(fcolors_.second) {
      mesh.remove_property_map(fcolors_.first);
    }
  }
  if(vcolors) {
    std::pair<Vcolors_map, bool> vcolors_ = 
      mesh.property_map<vertex_descriptor, std::string>("v:color");
    if(vcolors_.second) {
      mesh.remove_property_map(vcolors_.first);
    }
  }
  if(vnormals) {
    std::pair<Normals_map, bool> vnormals_ = 
      mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
    if(vnormals_.second) {
      mesh.remove_property_map(vnormals_.first);
    }
  }
}

std::pair<std::map<face_descriptor, std::string>, bool> copy_fcolor(EMesh3& mesh) {
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

std::pair<std::map<vertex_descriptor, std::string>, bool> copy_vcolor(EMesh3& mesh) {
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

std::pair<std::map<vertex_descriptor, Rcpp::NumericVector>, bool> copy_vnormal(EMesh3& mesh) {
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
    
  MaybeFcolorMap fcolorMap_ = copy_fcolor(tm);
  MaybeFcolorMap fcolorMap2_ = copy_fcolor(clipper);
  std::size_t nfaces = tm.number_of_faces();
  std::size_t nfaces_clipper = clipper.number_of_faces();
  
  ClipVisitor vis;
  Face_index_map fimap = 
    tm.add_property_map<face_descriptor, std::size_t>("f:i", 0).first;
  const bool doNotModify = !clipVolume;
  const bool clipping = PMP::clip(
    tm, clipper,
    PMP::parameters::clip_volume(clipVolume).visitor(vis).face_index_map(fimap),
    PMP::parameters::clip_volume(clipVolume).do_not_modify(doNotModify)
  );
  if(!clipping) {
    Rcpp::stop("Clipping has failed.");
  }

  Rcpp::Rcout << "Clipping ok\n";

  tm.collect_garbage();

  Rcpp::LogicalVector fremoved(tm.number_of_faces());
  std::pair<EMesh3::Property_map<face_descriptor, bool>, bool> maybemap = tm.property_map<face_descriptor, bool>("f:removed");
  Rcpp::Rcout << "fremoved: " << maybemap.second << "\n";
  if(maybemap.second) {
    int fr = 0;
    for(EMesh3::Face_index fi : tm.faces()) {
      fremoved(fr++) = maybemap.first[fi];
    }
  }

  Rcpp::Rcout << "nfaces clipper: " << clipper.number_of_faces() << "\n";
  Rcpp::Rcout << "clipper has garbage: " << clipper.has_garbage() << "\n";

  if(!clipVolume){
    if(fcolorMap_.second) {
      std::map<face_descriptor, std::string> fcolorMap = fcolorMap_.first;
      MapBetweenFaces fmap = *(vis.fmap_tm);
      Fcolors_map newfcolor = 
        tm.add_property_map<face_descriptor, std::string>("f:color", "").first;
      for(EMesh3::Face_index fi : tm.faces()) {
        face_descriptor fd = CGAL::SM_Face_index(fimap[fi]);
        if(size_t(fd) < nfaces) {
          newfcolor[fi] = fcolorMap[fd];
        } else {
          newfcolor[fi] = fcolorMap[fmap[fd]];
        }
      }
    }
    return Rcpp::List::create(Rcpp::Named("o") = 1);
  }

  if(fcolorMap_.second) {
    std::map<face_descriptor, std::string> fcolorMap = fcolorMap_.first;
    std::map<face_descriptor, std::string> fcolorMap2 = fcolorMap2_.first;
    MapBetweenFaces fmap_tm = *(vis.fmap_tm);
    MapBetweenFaces fmap_clipper = *(vis.fmap_clipper);
    MapBetweenFaces ftargets = *(vis.ftargets);
    std::map<face_descriptor, std::string> fcolormap_clipper;
    int fdi = 0;
    for(auto it = ftargets.rbegin(); it != ftargets.rend(); ++it) {
      face_descriptor fd = it->second;
      if(std::size_t(fd) >= nfaces_clipper) {
        fd = fmap_clipper[fd];
      }
      fcolormap_clipper[it->first] = fcolorMap2[fd];
      while(fimap[CGAL::SM_Face_index(fdi)] != 0) {
        fdi++;
      }
      fimap[CGAL::SM_Face_index(fdi++)] = std::size_t(it->first);
    }
    Fcolors_map newfcolor = 
      tm.add_property_map<face_descriptor, std::string>("f:color", "").first;
    for(EMesh3::Face_index fi : tm.faces()) {
      std::size_t ifi = fimap[fi];
      face_descriptor fd = CGAL::SM_Face_index(ifi);
      if(auto search = fmap_clipper.find(fd); search != fmap_clipper.end()) {
        Rcpp::Rcout << "1\n";
        newfcolor[fi] = fcolorMap[fd];
      } else if(auto search = ftargets.find(fd); search != ftargets.end()) {
        Rcpp::Rcout << "2\n";
        newfcolor[fi] = fcolormap_clipper[fd];
      } else if(ifi < nfaces) {
        Rcpp::Rcout << "3\n";
        newfcolor[fi] = fcolorMap[fd];
      } else {
        Rcpp::Rcout << "4\n";
        newfcolor[fi] = fcolorMap[fmap_tm[fd]];
      }
    }
  }
  Rcpp::Rcout << "DONE\n";

  Rcpp::StringVector Action((*(vis.action)).size());
  for(int i = 0; i < (*(vis.action)).size(); i++) {
    Action(i) = (*(vis.action))[i];
  }
  Rcpp::IntegerMatrix Fmap1((*(vis.fmap_tm)).size(), 2);
  int iii = 0;
  for(const auto& [fnew, fsplit] : *(vis.fmap_tm)) {
    Fmap1(iii++, Rcpp::_) = Rcpp::IntegerVector({int(fnew), int(fsplit)});
  }
  Rcpp::IntegerMatrix Fmap2((*(vis.fmap_clipper)).size(), 2);
  int iii2 = 0;
  for(const auto& [fnew, fsplit] : *(vis.fmap_clipper)) {
    Fmap2(iii2++, Rcpp::_) = Rcpp::IntegerVector({int(fnew), int(fsplit)});
  }
  Rcpp::IntegerMatrix xxx((*(vis.ftargets)).size(), 2);
  int ii = 0;
  for(const auto& [fsrc, ftgt] : *(vis.ftargets)) {
    xxx(ii++, Rcpp::_) = Rcpp::IntegerVector({int(fsrc), int(ftgt)});
  }
  Face_index_map fimapnew = tm.property_map<face_descriptor, std::size_t>("f:i").first;
  Rcpp::IntegerVector Fimap(tm.number_of_faces());
  int fff = 0;
  for(EMesh3::Face_index fi : tm.faces()) {
    Fimap(fff++) = fimapnew[fi];
  }
  Rcpp::IntegerVector Nfaces((*(vis.nfaces)).size());
  for(int i = 0; i < (*(vis.nfaces)).size(); i++) {
    Nfaces(i) = (*(vis.nfaces))[i];
  }
  Rcpp::IntegerVector Nfaces2((*(vis.nfaces2)).size());
  for(int i = 0; i < (*(vis.nfaces2)).size(); i++) {
    Nfaces2(i) = (*(vis.nfaces2))[i];
  }
  return Rcpp::List::create(
    Rcpp::Named("xxx") = xxx, 
    Rcpp::Named("fimap") = Fimap,
//    Rcpp::Named("fimap2") = Fimap2,
    Rcpp::Named("fmap1") = Fmap1,
    Rcpp::Named("fmap2") = Fmap2,
    Rcpp::Named("nfaces") = Nfaces,
    Rcpp::Named("nfaces2") = Nfaces2,
//    Rcpp::Named("zz") = ZZ,
    Rcpp::Named("action") = Action,
    Rcpp::Named("fremoved") = fremoved,
    Rcpp::Named("clipper") = Rcpp::XPtr<EMesh3>(new EMesh3(clipper), false)
  );

  
  // Rcpp::IntegerMatrix out((*(vis.fmap)).size(), 2);
  // //Rcpp::IntegerVector out2((*(vis.i)).size());
  // int i = 0;
  // for(const auto& [fnew, fsplit] : *(vis.fmap)) {
  //   out(i++, Rcpp::_) = Rcpp::IntegerVector({int(fsplit), int(fnew)});
  // }
  // for(int j = 0; j < (*(vis.i)).size(); j++) {
  //   out2(j) = (*(vis.i))[j];
  // }
  // return Rcpp::List::create(Rcpp::Named("out") = out, Rcpp::Named("out2") = out2);
  Rcpp::Rcout << clipper.number_of_faces() << "\n";
  
  // Rcpp::IntegerVector findices(mesh.number_of_faces());
  // int ff = 0;
  // for(EMesh3::Face_index fi : mesh.faces()) {
  //   findices(ff++) = fmap[fi];
  // }
  
  //PMP::merge_duplicated_vertices_in_boundary_cycles(mesh);
  
  //mesh.collect_garbage();
  
  // Rcpp::IntegerVector fconnectivity(mesh.number_of_faces());
  // std::pair<Face_index_map, bool> xmaybemap = mesh.property_map<face_descriptor, std::size_t>("f:connectivity");
  // Rcpp::Rcout << "fconnectivity: " << xmaybemap.second << "\n";
  // if(xmaybemap.second) {
  //   int ffff = 0;
  //   for(EMesh3::Face_index fi : mesh.faces()) {
  //     fconnectivity(ffff++) = xmaybemap.first[fi];
  //   }
  // }
  
  // Rcpp::LogicalVector fremoved(mesh.number_of_faces());
  // std::pair<EMesh3::Property_map<face_descriptor, bool>, bool> xxmaybemap = mesh.property_map<face_descriptor, bool>("f:removed");
  // Rcpp::Rcout << "fremoved: " << xxmaybemap.second << "\n";
  // if(xxmaybemap.second) {
  //   int fffff = 0;
  //   for(EMesh3::Face_index fi : mesh.faces()) {
  //     fremoved(fffff++) = xxmaybemap.first[fi];
  //   }
  // }
  
  
  // std::pair<EMesh3::Property_map<face_descriptor, std::size_t>, bool> maybemap = mesh.property_map<face_descriptor, std::size_t>("f:i");
  // Rcpp::IntegerVector findices2(mesh.number_of_faces());
  // if(maybemap.second) {
  //   int fff = 0;
  //   for(EMesh3::Face_index fi : mesh.faces()) {
  //     findices2(fff++) = maybemap.first[fi];
  //   }
  // }
  
//     Rcpp::Rcout << "PROPERTIES" << "\n";
//     std::vector<std::string> props = mesh.properties<face_descriptor>();
//     for(int p = 0; p < props.size(); p++) {
//       Rcpp::Rcout << props[p] << "\n";
//     }
  
//     Rcpp::LogicalVector keep(clipper.number_of_faces());
//     int f = 0;
//     for(EMesh3::Face_index fi : clipper.faces()) {
//       keep(f++) = !clipper.is_removed(fi);
//     }
  
//     return Rcpp::List::create(
//       Rcpp::Named("out") = out, 
//       Rcpp::Named("keep") = keep,
//       Rcpp::Named("fconnectivity") = fconnectivity,
//       Rcpp::Named("fremoved") = fremoved
    
// //      Rcpp::Named("findices") = findices,
// //      Rcpp::Named("findices2") = findices2
//     );
  // if(fcolors.isNotNull()) {
  //   Rcpp::StringVector fcolors0(fcolors);
  //   const size_t nfaces = mcopy.number_of_faces();
  //   std::vector<Triangle> triangles;
  //   triangles.reserve(nfaces);
  //   for(EMesh3::Face_index fd : mcopy.faces()) {
  //     auto vd = vertices_around_face(mcopy.halfedge(fd), mcopy).begin();
  //     triangles.emplace_back(Triangle(
  //       mcopy.point(*(++vd)), mcopy.point(*(++vd)), mcopy.point(*vd)
  //     ));
  //   }
  // 
  //   const size_t nf = mesh.number_of_faces();
  //   Rcpp::StringVector fcolors1(nf);
  //   size_t j = 0;
  //   for(EMesh3::Face_index f : mesh.faces()) {
  //     auto vd = vertices_around_face(mesh.halfedge(f), mesh).begin();
  //     Triangle tr(mesh.point(*(++vd)), mesh.point(*(++vd)), mesh.point(*vd));
  //     EPoint3 c = CGAL::centroid(tr);
  //     size_t k;
  //     for(k = 0; k < nfaces; k++) {
  //       if(triangles[k].has_on(c)) {
  //         break;
  //       }
  //     }
  //     fcolors1(j++) = fcolors0(k);
  //   }
  //   fcolors = Rcpp::Nullable<Rcpp::StringVector>(fcolors1);
  // }
}
