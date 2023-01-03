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
  
  // EMesh3 mcopy;
  // if(fcolors.isNotNull()) {
  //   CGAL::copy_face_graph(mesh, mcopy);
  // }
  
  std::pair<std::map<face_descriptor, std::string>, bool> fcolorstdmap_ = copy_fcolor(tm);
  std::pair<std::map<face_descriptor, std::string>, bool> fcolorstdmap2_ = copy_fcolor(clipper);
  // std::pair<Fcolors_map, bool> fcolorsmap_ = 
  //   tm.property_map<face_descriptor, std::string>("f:color");
  // const bool has_fcolor = fcolorsmap_.second;
  // std::map<face_descriptor, std::string> fcolorstdmap;
  // if(has_fcolor) {
  //   for(EMesh3::Face_index fi : tm.faces()) {
  //     fcolorstdmap[fi] = fcolorsmap_.first[fi];
  //   }
  //   tm.remove_property_map(fcolorsmap_.first);
  // }    
  std::size_t nfaces = tm.number_of_faces();
  
  ClipVisitor vis;
  Face_index_map fimap = 
    tm.add_property_map<face_descriptor, std::size_t>("f:i", 0).first;
  Face_index_map fimap2 = 
    clipper.add_property_map<face_descriptor, std::size_t>("f:i").first;
  // for(EMesh3::Face_index fi : clipper.faces()) {
  //   fimap2[fi] = 10000;
  // }
  const bool doNotModify = !clipVolume;
  const bool clipping = PMP::clip(
    tm, clipper,
    PMP::parameters::clip_volume(clipVolume).visitor(vis).face_index_map(fimap),
    PMP::parameters::clip_volume(clipVolume).do_not_modify(doNotModify).face_index_map(fimap2)
  );
  if(!clipping) {
    Rcpp::stop("Clipping has failed.");
  }

  Rcpp::Rcout << "Clipping ok\n";

  Face_index_map fimapprime = tm.property_map<face_descriptor, std::size_t>("f:i").first;
  Rcpp::IntegerVector ZZ(tm.number_of_faces());
  int zz = 0;
  for(EMesh3::Face_index fi : tm.faces()) {
    ZZ(zz++) = fimapprime[fi];
  }  
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

  
  // normals = R_NilValue;
  // vcolors = R_NilValue;

  Rcpp::Rcout << "nfaces clipper: " << clipper.number_of_faces() << "\n";
  Rcpp::Rcout << "clipper has garbage: " << clipper.has_garbage() << "\n";

/*   if(fcolorstdmap_.second && !clipVolume) {
    std::map<face_descriptor, std::string> fcolorstdmap = fcolorstdmap_.first;
    std::map<face_descriptor, face_descriptor> fmap = *(vis.fmap);
    Rcpp::Rcout << "Test fcolorstdmap \n";
    Rcpp::Rcout << fcolorstdmap[CGAL::SM_Face_index(1)] << "\n";
    //mesh.remove_property_map(fcolorsmap_.first);
    // Rcpp::Rcout << "Test f:color map \n";
    // std::pair<Fcolors_map, bool> fcolormap_ = mesh.property_map<face_descriptor, std::string>("f:color");
    // Rcpp::Rcout << fcolormap_.second;
    Rcpp::Rcout << "\n";
    Face_index_map fimapnew = tm.property_map<face_descriptor, std::size_t>("f:i").first;
    Fcolors_map newfcolorsmap = 
      tm.add_property_map<face_descriptor, std::string>("f:color").first;
//    Rcpp::Rcout << "new fcolor map:\n";
//    Rcpp::Rcout << newfcolorsmap[CGAL::SM_Face_index(1)] << "\n";
    for(EMesh3::Face_index fi : tm.faces()) {
      face_descriptor fd = CGAL::SM_Face_index(fimapnew[fi]);
      if(size_t(fd) < nfaces) {
        newfcolorsmap[fi] = fcolorstdmap[fd];
      } else {
        face_descriptor ffd = fmap[fd];
        //   Rcpp::Rcout << ffd << "\n";
        newfcolorsmap[fi] = fcolorstdmap[ffd];
      }
    }
    return Rcpp::List::create(Rcpp::Named("o") = 1);
  }
 */  
  if(fcolorstdmap_.second) {
    std::map<face_descriptor, std::string> fcolorstdmap = fcolorstdmap_.first;
    std::map<face_descriptor, std::string> fcolorstdmap2 = fcolorstdmap2_.first;
    std::map<face_descriptor, face_descriptor> fmap1 = *(vis.fmap1);
    std::map<face_descriptor, face_descriptor> fmap2 = *(vis.fmap2);
    std::map<face_descriptor, bool> istm = *(vis.istm);
    std::map<face_descriptor, face_descriptor> XXX = *(vis.xxx);
    Rcpp::Rcout << "Test fcolorstdmap \n";
    Rcpp::Rcout << fcolorstdmap[CGAL::SM_Face_index(1)] << "\n";
    Rcpp::Rcout << fcolorstdmap2[CGAL::SM_Face_index(1)] << "\n";
    //mesh.remove_property_map(fcolorsmap_.first);
    // Rcpp::Rcout << "Test f:color map \n";
    // std::pair<Fcolors_map, bool> fcolormap_ = mesh.property_map<face_descriptor, std::string>("f:color");
    // Rcpp::Rcout << fcolormap_.second;
    Rcpp::Rcout << "\n";
    Face_index_map fimapnew = tm.property_map<face_descriptor, std::size_t>("f:i").first;
    Fcolors_map newfcolorsmap = 
      tm.add_property_map<face_descriptor, std::string>("f:color").first;
//    Rcpp::Rcout << "new fcolor map:\n";
//    Rcpp::Rcout << newfcolorsmap[CGAL::SM_Face_index(1)] << "\n";
    for(EMesh3::Face_index fi : tm.faces()) {
      std::size_t ifi = fimapnew[fi];
      std::size_t ifi2 = fimap[fi];
      if(ifi != ifi2) {
        Rcpp::Rcout << "ifi != ifi2\n";
      }
      face_descriptor fd = CGAL::SM_Face_index(ifi);
      if(auto search = fmap1.find(fd); search != fmap1.end()) {
        face_descriptor ffff = fmap1[fd];
        newfcolorsmap[fi] = fcolorstdmap[ffff];
      } else if(auto search = fmap2.find(fd); search != fmap2.end()) {
        face_descriptor ffff = fmap2[fd];
        newfcolorsmap[fi] = fcolorstdmap2[ffff];
      } else {
        if(ifi == 0) {
          newfcolorsmap[fi] = fcolorstdmap2[fi];
        } else {
          newfcolorsmap[fi] = fcolorstdmap[fd];
        }
        // if(size_t(fd) < nfaces) {
        //   newfcolorsmap[fi] = fcolorstdmap[fd];
        // } else {
        //   face_descriptor ffd = fmap[fd];
        //   Rcpp::Rcout << ffd << "\n";
        //   newfcolorsmap[fi] = istm[ffd] ? fcolorstdmap[ffd] : fcolorstdmap2[ffd];
        // }
      }
    }
  }
  Rcpp::StringVector Action((*(vis.action)).size());
  for(int i = 0; i < (*(vis.action)).size(); i++) {
    Action(i) = (*(vis.action))[i];
  }
  Rcpp::IntegerMatrix out((*(vis.istm)).size(), 2);
  int i = 0;
  for(const auto& [fsplit, bb] : *(vis.istm)) {
    out(i++, Rcpp::_) = Rcpp::IntegerVector({int(fsplit), int(bb)});
  }
  Rcpp::IntegerMatrix Fmap1((*(vis.fmap1)).size(), 2);
  int iii = 0;
  for(const auto& [fnew, fsplit] : *(vis.fmap1)) {
    Fmap1(iii++, Rcpp::_) = Rcpp::IntegerVector({int(fnew), int(fsplit)});
  }
  Rcpp::IntegerMatrix Fmap2((*(vis.fmap2)).size(), 2);
  int iii2 = 0;
  for(const auto& [fnew, fsplit] : *(vis.fmap2)) {
    Fmap2(iii2++, Rcpp::_) = Rcpp::IntegerVector({int(fnew), int(fsplit)});
  }
  Rcpp::IntegerMatrix xxx((*(vis.xxx)).size(), 2);
  int ii = 0;
  for(const auto& [fsrc, ftgt] : *(vis.xxx)) {
    xxx(ii++, Rcpp::_) = Rcpp::IntegerVector({int(fsrc), int(ftgt)});
  }
  Face_index_map fimapnew = tm.property_map<face_descriptor, std::size_t>("f:i").first;
  Rcpp::IntegerVector Fimap(tm.number_of_faces());
  int fff = 0;
  for(EMesh3::Face_index fi : tm.faces()) {
    Fimap(fff++) = fimapnew[fi];
  }
  Face_index_map fimap2new = clipper.property_map<face_descriptor, std::size_t>("f:i").first;
  Rcpp::IntegerVector Fimap2(clipper.number_of_faces());
  int ff2 = 0;
  for(EMesh3::Face_index fi : clipper.faces()) {
    Fimap2(ff2++) = fimap2new[fi];
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
    Rcpp::Named("out") = out, 
    Rcpp::Named("xxx") = xxx, 
    Rcpp::Named("fimap") = Fimap,
    Rcpp::Named("fimap2") = Fimap2,
    Rcpp::Named("fmap1") = Fmap1,
    Rcpp::Named("fmap2") = Fmap2,
    Rcpp::Named("nfaces") = Nfaces,
    Rcpp::Named("nfaces2") = Nfaces2,
    Rcpp::Named("zz") = ZZ,
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
