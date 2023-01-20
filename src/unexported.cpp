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
    face_descriptor fd = mesh.add_face(face);
    if(fd == EMesh3::null_face()) {
      Rcpp::stop("Cannot add face " + std::to_string(i+1) + ".");
    }
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


template <typename Keytype, typename Valuetype>
void removeProperty(EMesh3& mesh, std::string propname) {
  std::pair<EMesh3::Property_map<Keytype, Valuetype>, bool> pmap_ = 
    mesh.property_map<Keytype, Valuetype>(propname);
  if(pmap_.second) {
    mesh.remove_property_map(pmap_.first);
  }
}

template void removeProperty<face_descriptor, Color>(
  EMesh3&, std::string
);
template void removeProperty<vertex_descriptor, Color>(
  EMesh3&, std::string
);
template void removeProperty<vertex_descriptor, EVector3>(
  EMesh3&, std::string
);


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

template MaybeFcolorMap copy_prop<face_descriptor, std::string>(
  EMesh3&, std::string
);
template MaybeFscalarMap copy_prop<face_descriptor, double>(
  EMesh3&, std::string
);
template MaybeVcolorMap copy_prop<vertex_descriptor, std::string>(
  EMesh3&, std::string
);
template MaybeVscalarMap copy_prop<vertex_descriptor, double>(
  EMesh3&, std::string
);
template MaybeNormalMap copy_prop<vertex_descriptor, Rcpp::NumericVector>(
  EMesh3&, std::string
);
template std::pair<std::map<face_descriptor, Color>, bool> 
  copy_prop<face_descriptor, Color>(EMesh3&, std::string);
template std::pair<std::map<vertex_descriptor, Color>, bool> 
  copy_prop<vertex_descriptor, Color>(EMesh3&, std::string);
template std::pair<std::map<vertex_descriptor, EVector3>, bool> 
  copy_prop<vertex_descriptor, EVector3>(EMesh3&, std::string);


template <
  typename SourceDescriptor, typename TargetDescriptor, typename Valuetype
>
void copy_property(
  EMesh3& mesh, 
  EMesh3& fmesh, 
  std::map<SourceDescriptor, TargetDescriptor> dmap, 
  std::string propname
) {
  std::pair<EMesh3::Property_map<TargetDescriptor, Valuetype>, bool> pmap_ = 
    mesh.property_map<TargetDescriptor, Valuetype>(propname);
  bool has_prop = pmap_.second;
  if(has_prop) {
    EMesh3::Property_map<TargetDescriptor, Valuetype> pmap = 
      fmesh.add_property_map<TargetDescriptor, Valuetype>(
        propname
      ).first;
    for(const auto& [source_descriptor, target_descriptor] : dmap) {
      pmap[target_descriptor] = 
        pmap_.first[TargetDescriptor(int(source_descriptor))];
    }
  }
}

template void copy_property<
  ffg_vertex_descriptor, vertex_descriptor, Rcpp::NumericVector
>(EMesh3&, EMesh3&, MapBetweenVertexDescriptors, std::string);
template void copy_property<
  ffg_vertex_descriptor, vertex_descriptor, std::string
>(EMesh3&, EMesh3&, MapBetweenVertexDescriptors, std::string);
template void copy_property<
  ffg_vertex_descriptor, vertex_descriptor, double
>(EMesh3&, EMesh3&, MapBetweenVertexDescriptors, std::string);
template void copy_property<
  ffg_face_descriptor, face_descriptor, std::string
>(EMesh3&, EMesh3&, MapBetweenFaceDescriptors, std::string);
template void copy_property<
  ffg_face_descriptor, face_descriptor, double
>(EMesh3&, EMesh3&, MapBetweenFaceDescriptors, std::string);


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
