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
    Face_index_map whichPart = 
      tm.add_property_map<face_descriptor, std::size_t>("f:which").first;
    Fcolors_map newfcolor = 
      tm.add_property_map<face_descriptor, std::string>("f:color", "").first;
    int nfaces_tmesh = 0;
    std::vector<std::string> fcolor_tmesh_vec;
    for(EMesh3::Face_index fi : tm.faces()) {
      std::size_t ifi = fimap[fi];
      face_descriptor fd = CGAL::SM_Face_index(ifi);
      if(auto search = fmap_clipper.find(fd); search != fmap_clipper.end()) {
        newfcolor[fi] = fcolorMap[fd];
        whichPart[fi] = 0;
        fcolor_tmesh_vec.push_back(fcolorMap[fd]);
        nfaces_tmesh++;
      } else if(auto search = ftargets.find(fd); search != ftargets.end()) {
        newfcolor[fi] = fcolormap_clipper[fd];
        whichPart[fi] = 1;
      } else if(ifi < nfaces) {
        newfcolor[fi] = fcolorMap[fd];
        whichPart[fi] = 0;
        fcolor_tmesh_vec.push_back(fcolorMap[fd]);
        nfaces_tmesh++;
      } else {
        newfcolor[fi] = fcolorMap[fmap_tm[fd]];
        whichPart[fi] = 0;
        fcolor_tmesh_vec.push_back(fcolorMap[fmap_tm[fd]]);
        nfaces_tmesh++;
      }
    }
    Rcpp::Rcout << "DONE\n";

    std::size_t sb = PMP::stitch_borders(tm);
    Rcpp::Rcout << "stitched borders: " << sb << "\n";

    std::vector<std::vector<std::size_t>> polygons;
    polygons.reserve(nfaces_tmesh);
    for(EMesh3::Face_index fi : tm.faces()) {
      if(whichPart[fi] == 0) {
        std::vector<std::size_t> face;
        for(EMesh3::Vertex_index vi :
            vertices_around_face(tm.halfedge(fi), tm)) {
          face.push_back(std::size_t(vi));
        }
        polygons.emplace_back(face);
      }
    }
    std::vector<EPoint3> points;
    points.reserve(tm.number_of_vertices());
    for(EMesh3::Vertex_index vi : tm.vertices()) {
      points.emplace_back(tm.point(vi));
    }

    SoupVisitor soupvis;
    bool success = PMP::orient_polygon_soup(points, polygons, CGAL::parameters::visitor(soupvis));
    Rcpp::Rcout << "orientation successful: " << success << "\n";

    EMesh3 tmesh;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, tmesh);
    Fcolors_map fcolor_tmesh = 
      tmesh.add_property_map<face_descriptor, std::string>("f:color", "").first;
    for(EMesh3::Face_index fi : tmesh.faces()) {
      fcolor_tmesh[fi] = fcolor_tmesh_vec[int(fi)];
    }

    // Face_index_map fm =
    //   tm.add_property_map<face_descriptor, std::size_t>("f:index").first;
    // Vertex_index_map vm =
    //   tm.add_property_map<vertex_descriptor, std::size_t>("v:index").first;
    // Halfedge_index_map hm =
    //   tm.add_property_map<halfedge_descriptor, std::size_t>("h:index").first;
    //Filtered_mesh ffg_tmesh(tm, 0, whichPart, PMP::parameters::face_index_map(fm).vertex_index_map(vm).halfedge_index_map(hm));
    //Filtered_graph ffg_clipper(tm, 1, whichPart);
    // EMesh3 tmesh = cloneMesh(tmbis, false, false, false);
    // for(int i = 0; i < tm.number_of_faces(); i++) {
    //   face_descriptor fd = CGAL::SM_Face_index(i);
    //   if(whichPart[fd] == 1) {
    //     tmesh.remove_face(fd);
    //   }
    // }
    // std::vector<EMesh3::Face_index> faces_to_remove;
    // std::vector<EMesh3::Vertex_index> vertices_to_remove;
    // for(EMesh3::Vertex_index vi : tmesh.vertices()) {
    //   EMesh3::Face_index fi = tmesh.face(tmesh.halfedge(vi));
    //   if(whichPart[fi] == 1) {
    //     if(!tmesh.is_removed(fi)) {
    //       faces_to_remove.push_back(fi);
    //     }
    //     //EMesh3::Vertex_index vi = tmesh.vertex(tmesh.edge(hi), 0);
    //     vertices_to_remove.push_back(vi);
    //   }
    // }
    // for(std::vector<EMesh3::Vertex_index>::iterator it = vertices_to_remove.begin(); it != vertices_to_remove.end(); ++it) {
    //   tmesh.remove_vertex(*it);
    // }
    // for(std::vector<EMesh3::Face_index>::iterator it = faces_to_remove.begin(); it != faces_to_remove.end(); ++it) {
    //   tmesh.remove_face(*it);
    // }
    Rcpp::Rcout << "number of faces: " << tmesh.number_of_faces() << "\n";
    Rcpp::Rcout << "number of vertices: " << tmesh.number_of_vertices() << "\n";
    Rcpp::Rcout << "number of removed faces: " << tmesh.number_of_removed_faces() << "\n";
    // Rcpp::LogicalVector fremoved(tmesh.number_of_faces());
    // std::pair<EMesh3::Property_map<face_descriptor, bool>, bool> maybe_fremoved = 
    //   tmesh.property_map<face_descriptor, bool>("f:removed");
    // Rcpp::Rcout << "fremoved present: " << maybe_fremoved.second << "\n";
    // if(maybe_fremoved.second) {
    //   int fr = 0;
    //   for(EMesh3::Face_index fi : tmesh.faces()) {
    //     fremoved(fr++) = maybe_fremoved.first[fi];
    //   }
    // }
    // if(Rcpp::is_false(Rcpp::any(fremoved))) {
    //   Rcpp::Rcout << "nothing removed\n";
    //   for(EMesh3::Face_index fi : tmesh.faces()) {
    //     if(whichPart[fi] == 1) {
    //       maybe_fremoved.first[fi] = true;
    //     }
    //   }
    // }
    // int fr = 0;
    // for(EMesh3::Face_index fi : tmesh.faces()) {
    //   fremoved(fr++) = tmesh.is_removed(fi);
    // }
    // if(Rcpp::is_false(Rcpp::any(fremoved))) {
    //   Rcpp::Rcout << "still nothing removed\n";
    // }
    // Rcpp::Rcout << "collect garbage\n";
    // tmesh.collect_garbage();
    //Rcpp::Rcout << "number of removed faces: " << tmesh.number_of_removed_faces() << "\n";


    EMesh3 tmesh2 = cloneMesh(tm, false, false, true);
    for(EMesh3::Face_index fi : tm.faces()) {
      if(whichPart[fi] == 0) {
        tmesh2.remove_face(fi);
      }
    }
    tmesh2.collect_garbage();

    Rcpp::LogicalVector Components(tm.number_of_faces());
    int fint = 0;
    for(EMesh3::Face_index fi : tm.faces()) {
      Components(fint++) = whichPart[fi] == 0;
    }

    std::size_t nvisolated = PMP::remove_isolated_vertices(tmesh);
    Rcpp::Rcout << "removed " << nvisolated << " isolated vertices\n";
    if(nvisolated > 0) {
      tmesh.collect_garbage();
    }
    std::size_t nvisolated2 = PMP::remove_isolated_vertices(tmesh2);
    Rcpp::Rcout << "removed " << nvisolated2 << " isolated vertices\n";

    std::vector<halfedge_descriptor> hds;
    std::back_insert_iterator<std::vector<halfedge_descriptor>> bii = std::back_inserter(hds);
    PMP::non_manifold_vertices(tmesh, bii);
    Rcpp::Rcout << "NMV: " << hds.size() << "\n";
    Rcpp::IntegerVector NMV(hds.size());
    for(int i = 0; i < hds.size(); i++) {
      NMV(i) = int(tmesh.target(hds[i]));
    }
    std::vector<halfedge_descriptor> hds2;
    PMP::non_manifold_vertices(tmesh2, std::back_insert_iterator<std::vector<halfedge_descriptor>>(hds2));
    Rcpp::Rcout << "NMV2: " << hds2.size() << "\n";
    Rcpp::IntegerVector NMV2(hds2.size());
    for(int i = 0; i < hds2.size(); i++) {
      NMV2(i) = int(tmesh2.target(hds2[i]));
    }


    return Rcpp::List::create(
      Rcpp::Named("mesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(tmesh), false),
      Rcpp::Named("mesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(tmesh2), false),
      Rcpp::Named("components") = Components,
      Rcpp::Named("nmv") = NMV,
      Rcpp::Named("nmv2") = NMV2
    );
  }



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
