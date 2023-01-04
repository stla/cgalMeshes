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
  EMesh3& mesh, const bool vnormals, const bool vcolors, const bool fcolors
) {
  EMesh3 out;
  CGAL::copy_face_graph(mesh, out);
  if(fcolors) {
    std::pair<Fcolors_map, bool> fcolors_ = 
      mesh.property_map<face_descriptor, std::string>("f:color");
    if(fcolors_.second) {
      Fcolors_map fcolorsmap = 
        out.add_property_map<face_descriptor, std::string>(
          "f:color", ""
        ).first;
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
        out.add_property_map<vertex_descriptor, std::string>(
          "v:color", ""
        ).first;
      for(EMesh3::Vertex_index vi : out.vertices()) {
        vcolorsmap[vi] = vcolors_.first[vi];
      }
    }
  }
  if(vnormals) {
    std::pair<Normals_map, bool> vnormals_ = 
      mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
    if(vnormals_.second) {
      Rcpp::NumericVector def = defaultNormal();
      Normals_map vnormalsmap = 
        out.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
          "v:normal", def
        ).first;
      for(EMesh3::Vertex_index vi : out.vertices()) {
        vnormalsmap[vi] = vnormals_.first[vi];
      }
    }
  }
  return out;
}

void removeProperties(
  EMesh3& mesh, const bool vnormals, const bool vcolors, const bool fcolors
) {
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

MaybeFcolorMap copy_fcolor(EMesh3& mesh) {
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

void triangulateMesh(EMesh3& mesh) {
  MaybeFcolorMap fcolormap_ = copy_fcolor(mesh);
  const bool hasFcolors = fcolormap_.second;
  removeProperties(mesh, true, false, false);
  TriangulateVisitor vis;
  const bool success = 
    PMP::triangulate_faces(mesh, CGAL::parameters::visitor(vis));
  if(!success) {
    Rcpp::stop("Triangulation has failed.");
  }
  if(hasFcolors) {
    MapBetweenFaces fmap = *(vis.fmap);
    Fcolors_map fcolors = 
      mesh.add_property_map<face_descriptor, std::string>("f:color", "").first;
    for(EMesh3::Face_index fi: mesh.faces()) {
      fcolors[fi] = fcolormap_.first[fmap[fi]];
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
    
  MaybeFcolorMap fcolorMap_ = copy_fcolor(tm);
  MaybeFcolorMap fcolorMap2_ = copy_fcolor(clipper);
  const bool hasColors = fcolorMap_.second && fcolorMap2_.second;

  std::size_t nfaces = tm.number_of_faces();
  std::size_t nfaces_clipper = clipper.number_of_faces();
  
  ClipVisitor vis;
  Face_index_map fimap = 
    tm.add_property_map<face_descriptor, std::size_t>("f:i", 0).first;
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

  Rcpp::Rcout << "Clipping ok\n";

  tm.collect_garbage();

  std::size_t sb = PMP::stitch_borders(tm);
  if(sb > 0) {
    std::string msg;
    if(sb == 1) {
      msg = "Stitched one border.\n";
    } else {
      msg = "Stitched " + std::to_string(sb) + " borders.";
    }
    Message(msg);
  }

  Rcpp::Rcout << "nfaces clipper: " << clipper.number_of_faces() << "\n";
  Rcpp::Rcout << "clipper has garbage: " << clipper.has_garbage() << "\n";

  if(!clipVolume){
    if(hasColors) {
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
    tm.remove_property_map(fimap);
    return Rcpp::List::create();
  }

  MapBetweenFaces fmap_tm = *(vis.fmap_tm);
  MapBetweenFaces fmap_clipper = *(vis.fmap_clipper);
  MapBetweenFaces ftargets = *(vis.ftargets);
  std::map<face_descriptor, std::string> fcolorMap;
  std::map<face_descriptor, std::string> fcolorMap2;
  std::map<face_descriptor, std::string> fcolormap_clipper;
  if(hasColors) {
    fcolorMap = fcolorMap_.first;
    fcolorMap2 = fcolorMap2_.first;
  }
  int fdi = 0;
  for(auto it = ftargets.rbegin(); it != ftargets.rend(); ++it) {
    face_descriptor fd = it->second;
    if(std::size_t(fd) >= nfaces_clipper) {
      fd = fmap_clipper[fd];
    }
    if(hasColors) {
      fcolormap_clipper[it->first] = fcolorMap2[fd];
    }
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
  int nfaces_tmesh2 = 0;
  std::vector<std::string> fcolor_tmesh_vec;
  std::vector<std::string> fcolor_tmesh2_vec;
  for(EMesh3::Face_index fi : tm.faces()) {
    std::size_t ifi = fimap[fi];
    face_descriptor fd = CGAL::SM_Face_index(ifi);
    if(auto search = fmap_clipper.find(fd); search != fmap_clipper.end()) {
      whichPart[fi] = 0;
      if(hasColors) {
        std::string color = fcolorMap[fd]; 
        newfcolor[fi] = color;
        fcolor_tmesh_vec.push_back(color);
      }
      nfaces_tmesh++;
    } else if(auto search = ftargets.find(fd); search != ftargets.end()) {
      whichPart[fi] = 1;
      if(hasColors) {
        std::string color = fcolormap_clipper[fd]; 
        newfcolor[fi] = color;
        fcolor_tmesh2_vec.push_back(color);
      }
      nfaces_tmesh2++;
    } else if(ifi < nfaces) {
      whichPart[fi] = 0;
      if(hasColors) {
        std::string color = fcolorMap[fd]; 
        newfcolor[fi] = color;
        fcolor_tmesh_vec.push_back(color);
      }
      nfaces_tmesh++;
    } else {
      whichPart[fi] = 0;
      if(hasColors) {
        std::string color = fcolorMap[fmap_tm[fd]]; 
        newfcolor[fi] = color;
        fcolor_tmesh_vec.push_back(color);
      }
      nfaces_tmesh++;
    }
  }

  std::vector<EPoint3> points;
  points.reserve(tm.number_of_vertices());
  for(EMesh3::Vertex_index vi : tm.vertices()) {
    points.emplace_back(tm.point(vi));
  }

  EMesh3 tmesh;
  {
    std::vector<std::vector<std::size_t>> polygons;
    polygons.reserve(nfaces_tmesh);
    for(EMesh3::Face_index fi : tm.faces()) {
      if(whichPart[fi] == 0) {
        std::vector<std::size_t> face;
        face.reserve(3);
        for(EMesh3::Vertex_index vi :
            vertices_around_face(tm.halfedge(fi), tm)) {
          face.emplace_back(std::size_t(vi));
        }
        polygons.emplace_back(face);
      }
    }
    SoupVisitor soupvis;
    bool success = PMP::orient_polygon_soup(
      points, polygons, CGAL::parameters::visitor(soupvis)
    );
    if(!success) {
      Rcpp::warning("Polygons orientation failed.");
    }
    PMP::polygon_soup_to_polygon_mesh(points, polygons, tmesh);
    if(hasColors) {
      Fcolors_map fcolor_tmesh = 
        tmesh.add_property_map<face_descriptor, std::string>(
          "f:color", ""
        ).first;
      for(EMesh3::Face_index fi : tmesh.faces()) {
        fcolor_tmesh[fi] = fcolor_tmesh_vec[int(fi)];
      }
    }
  }

  EMesh3 tmesh2;
  {
    std::vector<std::vector<std::size_t>> polygons;
    polygons.reserve(nfaces_tmesh2);
    for(EMesh3::Face_index fi : tm.faces()) {
      if(whichPart[fi] == 1) {
        std::vector<std::size_t> face;
        face.reserve(3);
        for(EMesh3::Vertex_index vi :
            vertices_around_face(tm.halfedge(fi), tm)) {
          face.emplace_back(std::size_t(vi));
        }
        polygons.emplace_back(face);
      }
    }
    SoupVisitor soupvis;
    bool success = PMP::orient_polygon_soup(
      points, polygons, CGAL::parameters::visitor(soupvis)
    );
    if(!success) {
      Rcpp::warning("Polygons orientation failed.");
    }
    PMP::polygon_soup_to_polygon_mesh(points, polygons, tmesh2);
    if(hasColors) {
      Fcolors_map fcolor_tmesh2 = 
        tmesh2.add_property_map<face_descriptor, std::string>(
          "f:color", ""
        ).first;
      for(EMesh3::Face_index fi : tmesh2.faces()) {
        fcolor_tmesh2[fi] = fcolor_tmesh2_vec[int(fi)];
      }
    }
  }

  Rcpp::LogicalVector Components(tm.number_of_faces());
  int fint = 0;
  for(EMesh3::Face_index fi : tm.faces()) {
    Components(fint++) = whichPart[fi] == 0;
  }

  {
    std::size_t nvisolated = PMP::remove_isolated_vertices(tmesh);
    if(nvisolated > 0) {
      std::string msg;
      if(nvisolated == 1) {
        msg = "Removed one isolated vertex.";
      } else {
        msg = "Removed " + std::to_string(nvisolated) + " isolated vertices.";
      }
      Message(msg);
      tmesh.collect_garbage();
    }
  }
  {
    std::size_t nvisolated = PMP::remove_isolated_vertices(tmesh2);
    if(nvisolated > 0) {
      std::string msg;
      if(nvisolated == 1) {
        msg = "Removed one isolated vertex.";
      } else {
        msg = "Removed " + std::to_string(nvisolated) + " isolated vertices.";
      }
      Message(msg);
      tmesh2.collect_garbage();
    }
  }

  tm.remove_property_map(fimap);

  {
    std::vector<halfedge_descriptor> hds;
    std::back_insert_iterator<std::vector<halfedge_descriptor>> bii = std::back_inserter(hds);
    PMP::non_manifold_vertices(tmesh, bii);
    if(hds.size() > 0) {
      std::string msg;
      if(hds.size() == 1) {
        msg = "Found one non-manifold vertex.";
      } else {
        msg = "Found " + std::to_string(hds.size()) + " non-manifold vertices.";
      }
      Message(msg);
    }
  }
  {
    std::vector<halfedge_descriptor> hds;
    std::back_insert_iterator<std::vector<halfedge_descriptor>> bii = std::back_inserter(hds);
    PMP::non_manifold_vertices(tmesh2, bii);
    if(hds.size() > 0) {
      std::string msg;
      if(hds.size() == 1) {
        msg = "Found one non-manifold vertex.";
      } else {
        msg = "Found " + std::to_string(hds.size()) + " non-manifold vertices.";
      }
      Message(msg);
    }
  }

  Rcpp::List meshes = Rcpp::List::create(
    Rcpp::Named("mesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(tmesh), false),
    Rcpp::Named("mesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(tmesh2), false)
  );
  meshes.attr("components") = Components;

  return meshes;

  /////////////////////////////////////////////////////////////////////

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
    // Rcpp::Named("fremoved") = fremoved,
    Rcpp::Named("clipper") = Rcpp::XPtr<EMesh3>(new EMesh3(clipper), false)
  );
  
  //PMP::merge_duplicated_vertices_in_boundary_cycles(mesh);
}
