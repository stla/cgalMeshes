#ifndef _HEADER_
#include "cgalMesh.h"
#endif

class CGALmesh {
public:
  EMesh3 mesh;
  Rcpp::List mymesh;
  Rcpp::XPtr<EMesh3> xptr;
  Rcpp::Nullable<Rcpp::NumericMatrix> normals;
  Rcpp::Nullable<Rcpp::StringVector> vcolors;
  Rcpp::Nullable<Rcpp::StringVector> fcolors;
  CGALmesh(const Rcpp::NumericMatrix vertices,
           const Rcpp::List faces,
           bool soup,
           const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_,
           const Rcpp::Nullable<Rcpp::StringVector> &vcolors_,
           const Rcpp::Nullable<Rcpp::StringVector> &fcolors_)
    : mesh(
        makeMesh(
          vertices, faces, soup, normals_, vcolors_, fcolors_
        )
      ),
      mymesh(
        Rcpp::List::create(
          Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(&mesh, false),
          Rcpp::Named("normals") = normals_,
          Rcpp::Named("vcolors") = vcolors_,
          Rcpp::Named("fcolors") = fcolors_
        )
      ), 
      xptr(Rcpp::XPtr<EMesh3>(&mesh, false)),
      normals(normals_),
      vcolors(vcolors_),
      fcolors(fcolors_) {}
  
  CGALmesh(Rcpp::XPtr<EMesh3> xptr_)
    : mesh(*(xptr_.get())),
      mymesh(
        Rcpp::List::create(
          Rcpp::Named("xptr") = xptr_,
          Rcpp::Named("normals") = R_NilValue,
          Rcpp::Named("vcolors") = R_NilValue,
          Rcpp::Named("fcolors") = R_NilValue
        )
      ),
      xptr(xptr_),
      normals(R_NilValue),
      vcolors(R_NilValue),
      fcolors(R_NilValue) {}
  
  CGALmesh(const std::string filename, const bool binary)
    : mesh(readMeshFile(filename)), 
      mymesh(
        Rcpp::List::create(
          Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(&mesh, false),
          Rcpp::Named("normals") = R_NilValue,
          Rcpp::Named("vcolors") = R_NilValue,
          Rcpp::Named("fcolors") = R_NilValue
        )
      ),
      xptr(Rcpp::XPtr<EMesh3>(&mesh, false)),
      normals(R_NilValue),
      vcolors(R_NilValue),
      fcolors(R_NilValue) {}


  double area() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The mesh self-intersects.");
    }
    const EK::FT ar = PMP::area(mesh);
    return CGAL::to_double<EK::FT>(ar);
  }


  void assignFaceColors(Rcpp::StringVector colors) {
    if(colors.size() != 1 && colors.size() != mesh.number_of_faces()) {
      Rcpp::stop(
        "The number of colors does not match the number of faces."
      );
    }
    removeProperties(mesh, {"f:color"});
    Fcolors_map fcolor = 
      mesh.add_property_map<face_descriptor, std::string>(
        "f:color", ""
      ).first;
    if(colors.size() == 1) {
      for(EMesh3::Face_index fi : mesh.faces()) {
        fcolor[fi] = colors(0);
      }
    } else {
      int i = 0;
      for(EMesh3::Face_index fi : mesh.faces()) {
        fcolor[fi] = colors(i++);
      }
    }
  }


  void assignFaceScalars(Rcpp::NumericVector scalars) {
    if(scalars.size() != mesh.number_of_faces()) {
      Rcpp::stop(
        "The number of scalars does not match the number of faces."
      );
    }
    removeProperties(mesh, {"f:scalar"});
    Fscalars_map fscalar = 
      mesh.add_property_map<face_descriptor, double>(
        "f:scalar", nan("")
      ).first;
    int i = 0;
    for(EMesh3::Face_index fi : mesh.faces()) {
      fscalar[fi] = scalars(i++);
    }
  }


  void assignVertexColors(Rcpp::StringVector colors) {
    if(colors.size() != mesh.number_of_vertices()) {
      Rcpp::stop(
        "The number of colors does not match the number of vertices."
      );
    }
    removeProperties(mesh, {"v:color"});
    Vcolors_map vcolor = 
      mesh.add_property_map<vertex_descriptor, std::string>(
        "v:color", ""
      ).first;
    int i = 0;
    for(EMesh3::Vertex_index vi : mesh.vertices()) {
      vcolor[vi] = colors(i++);
    }
  }


  void assignVertexScalars(Rcpp::NumericVector scalars) {
    if(scalars.size() != mesh.number_of_vertices()) {
      Rcpp::stop(
        "The number of scalars does not match the number of vertices."
      );
    }
    removeProperties(mesh, {"v:scalar"});
    Vscalars_map vscalar = 
      mesh.add_property_map<vertex_descriptor, double>(
        "v:scalar", nan("")
      ).first;
    int i = 0;
    for(EMesh3::Vertex_index vi : mesh.vertices()) {
      vscalar[vi] = scalars(i++);
    }
  }


  Rcpp::NumericVector centroid() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    Mesh3 epickCopy;
    CGAL::copy_face_graph(mesh, epickCopy);
    const Point3 centroid = PMP::centroid(epickCopy);
    // const EPoint3 centroid = PMP::centroid(mesh);
    Rcpp::NumericVector out(3);
    out(0) = centroid.x();
    out(1) = centroid.y();
    out(2) = centroid.z();
    // out(0) = CGAL::to_double<EK::FT>(centroid.x());
    // out(1) = CGAL::to_double<EK::FT>(centroid.y());
    // out(2) = CGAL::to_double<EK::FT>(centroid.z());
    return out;
  }


  void testsplit(Rcpp::XPtr<EMesh3> clipperXPtr) {
    EMesh3 splitter = *(clipperXPtr.get());
    PMP::split(mesh, splitter);
    mesh.collect_garbage();
//     MyVisitor vis;
//     PMP::split(mesh, splitter, CGAL::parameters::visitor(vis));
//     mesh.collect_garbage();
//     Rcpp::IntegerMatrix out((*(vis.vmap)).size(), 2);
// //    Rcpp::IntegerVector out2((*(vis.i)).size());
//     int i = 0;
//     for(const auto& [fsplit, fnew] : *(vis.vmap)) {
//       out(i++, Rcpp::_) = Rcpp::IntegerVector({int(fsplit), int(fnew)});
//     }
//     // for(int j = 0; j < (*(vis.i)).size(); j++) {
//     //   out2(j) = (*(vis.i))[j];
//     // }
//     return out;
  }


  Rcpp::List clipMesh(Rcpp::XPtr<EMesh3> clipperXPtr, const bool clipVolume) {
    EMesh3 clipper = *(clipperXPtr.get());
    return clipping(mesh, clipper, clipVolume);
  }


  void clipToPlane(
    Rcpp::NumericVector planePoint, 
    Rcpp::NumericVector planeNormal, 
    const bool clipVolume
  ) {
    EPoint3 point(planePoint[0], planePoint[1], planePoint[2]);
    EVector3 normal(planeNormal[0], planeNormal[1], planeNormal[2]);
    EPlane3 plane(point, normal);
    clippingToPlane(mesh, plane, clipVolume);
  }


  Rcpp::XPtr<EMesh3> doubleclip(Rcpp::XPtr<EMesh3> clipperXPtr) { // to remove
    EMesh3 meshcopy = cloneMesh(mesh, {"f:color"});
    EMesh3 clipper = *(clipperXPtr.get());
    EMesh3 clippercopy = cloneMesh(clipper, {"f:color"});
     
    clipping(mesh, clipper, false);
    clipping(clippercopy, meshcopy, false);

    const int nfaces = mesh.number_of_faces();
    const int nvertices = mesh.number_of_vertices();

    std::pair<Fcolors_map, bool> fcolorsmap1_ = 
      mesh.property_map<face_descriptor, std::string>("f:color");
    std::pair<Fcolors_map, bool> fcolorsmap2_ = 
      clippercopy.property_map<face_descriptor, std::string>("f:color");

    EMesh3 out;
    for(EMesh3::Vertex_index vi : mesh.vertices()) {
      out.add_vertex(mesh.point(vi));
    }
    for(EMesh3::Vertex_index vi : clippercopy.vertices()) {
      out.add_vertex(clippercopy.point(vi));
    }
    for(EMesh3::Face_index fi : mesh.faces()) {
      auto it = vertices_around_face(mesh.halfedge(fi), mesh);
      std::vector<EMesh3::Vertex_index> face(it.begin(), it.end());
      out.add_face(face);
    }
    for(EMesh3::Face_index fi : clippercopy.faces()) {
      auto it = vertices_around_face(clippercopy.halfedge(fi), clippercopy);
      std::vector<EMesh3::Vertex_index> face;
      face.reserve(3);
      for(EMesh3::Vertex_index v: it) {
        face.emplace_back(CGAL::SM_Vertex_index(nvertices + int(v)));
      }
      out.add_face(face);
    }
    if(fcolorsmap1_.second && fcolorsmap2_.second) {
      Fcolors_map fcolorsmap = 
        out.add_property_map<face_descriptor, std::string>("f:color").first;
      for(EMesh3::Face_index fi : mesh.faces()) {
        fcolorsmap[fi] = fcolorsmap1_.first[fi];
      }
      for(EMesh3::Face_index fi : clippercopy.faces()) {
        fcolorsmap[CGAL::SM_Face_index(nfaces + int(fi))] = fcolorsmap2_.first[fi];
      }
    }

    halfedge_descriptor h;
    for(EMesh3::Face_index fi : out.faces()) {
      h = out.halfedge(fi);
      if(out.is_border(h)) {
        break;
      }
    }

    std::vector<face_descriptor> fs;
    PMP::triangulate_hole(out, h, std::back_insert_iterator<std::vector<face_descriptor>>(fs));
    for(int i = 0; i < fs.size(); i++) {
      Rcpp::Rcout << fs[i] << "\n";
    }    

    return Rcpp::XPtr<EMesh3>(new EMesh3(out), false);
  }

  
  Rcpp::XPtr<EMesh3> clone() {
    EMesh3 copy = cloneMesh(mesh, {"f:color", "v:color", "v:normal"});
    return Rcpp::XPtr<EMesh3>(new EMesh3(copy), false);
  }


  void computeNormals() {
    std::pair<CGALnormals_map, bool> vnormals_ = 
      mesh.property_map<vertex_descriptor, EVector3>("v:normals");
    if(vnormals_.second) {
      mesh.remove_property_map(vnormals_.first);
    }
    CGALnormals_map vnormals = 
      mesh.add_property_map<vertex_descriptor, EVector3>(
                          "v:normals", CGAL::NULL_VECTOR
                        ).first;
    PMP::compute_vertex_normals(mesh, vnormals);
    removeProperties(mesh, {"v:normal"});
    Normals_map vnormal_map = 
      mesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
                          "v:normal", defaultNormal()
                        ).first;
    for(EMesh3::Vertex_index vi : mesh.vertices()) {
      Rcpp::NumericVector rcppnormal(3);
      const EVector3 normal = vnormals[vi];
      rcppnormal(0) = CGAL::to_double<EK::FT>(normal.x());
      rcppnormal(1) = CGAL::to_double<EK::FT>(normal.y());
      rcppnormal(2) = CGAL::to_double<EK::FT>(normal.z());
      vnormal_map[vi] = rcppnormal;
    }
  }


  Rcpp::List connectedComponents(const bool triangulate) {

    const bool really_triangulate = 
      triangulate && !CGAL::is_triangle_mesh(mesh);
    EMesh3 tmesh = cloneMesh(
      mesh, {"f:color", "v:color", "f:scalar", "v:scalar", "v:normal"}
    );
    if(really_triangulate) {
      triangulateMesh(tmesh);
    }

    std::pair<Normals_map, bool> vnormals_ = 
      tmesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
    std::pair<Vcolors_map, bool> vcolors_ = 
      tmesh.property_map<vertex_descriptor, std::string>("v:color");
    std::pair<Fcolors_map, bool> fcolors_ = 
      tmesh.property_map<face_descriptor, std::string>("f:color");
    std::pair<Vscalars_map, bool> vscalars_ = 
      tmesh.property_map<vertex_descriptor, double>("v:scalar");
    std::pair<Fscalars_map, bool> fscalars_ = 
      tmesh.property_map<face_descriptor, double>("f:scalar");
    const bool hasNormals = vnormals_.second && !really_triangulate;
    const bool hasVcolors = vcolors_.second;
    const bool hasFcolors = fcolors_.second;
    const bool hasVscalars = vscalars_.second;
    const bool hasFscalars = fscalars_.second;
    const bool hasProps = 
      hasNormals || hasVcolors || hasFcolors || hasVscalars || hasFscalars;
    // EMesh3::Property_map<face_descriptor, std::vector<vertex_descriptor>> faceVertices = 
    //   tmesh.property_map<face_descriptor, std::vector<vertex_descriptor>>(
    //     "f:vertices"
    //   ).first;

    /*
    If there's no normals and no vertex colors, we use connected components 
    decomposition made with components of face indices.
    */
    Rcpp::Rcout << "starting\n";
    if(true){ //!hasNormals && !hasVcolors) {
      Face_index_map fccmap = 
        tmesh.add_property_map<face_descriptor, std::size_t>("f:CC", 0).first;
      const std::size_t ncc = PMP::connected_components(tmesh, fccmap);
      // std::vector<EMesh3> cc_meshes;
      // PMP::split_connected_components(tmesh, cc_meshes, CGAL::parameters::face_patch_map(fccmap));
      // Rcpp::Rcout << "splitting ok\n";
      // const std::size_t ncc = cc_meshes.size();
      if(ncc == 1) {
        Message("Only one component found.\n");
      } else {
        const std::string msg = "Found " + std::to_string(ncc) + " components.\n";
        Message(msg);
      }

      std::vector<MapBetweenFaces> fcomponents(ncc);
      std::vector<int> counters(ncc, 0);
      if(hasProps) {
        for(EMesh3::Face_index fi : tmesh.faces()) {
          std::size_t c = fccmap[fi];
          face_descriptor cf = CGAL::SM_Face_index(counters[c]++);
          fcomponents[c].insert(std::make_pair(cf, fi));
        }
      }

      Fcolors_map fcolors;
      if(hasFcolors) {
        fcolors = fcolors_.first;
      }
      Vcolors_map vcolors;
      if(hasVcolors) {
        vcolors = vcolors_.first;
      }
      Fscalars_map fscalars;
      if(hasFscalars) {
        fscalars = fscalars_.first;
      }
      Vscalars_map vscalars;
      if(hasVscalars) {
        vscalars = vscalars_.first;
      }
      Normals_map vnormals;
      if(hasNormals) {
        vnormals = vnormals_.first;
      }

      Rcpp::List xptrs(ncc);

      for(std::size_t c = 0; c < ncc; c++) {
//        EMesh3 cmesh = cc_meshes[c];
        Filtered_graph ffg(tmesh, c, fccmap);
        MapBetweenVertexDescriptors v2vmap_;
        boost::associative_property_map<MapBetweenVertexDescriptors> v2vmap(v2vmap_);
        MapBetweenFaceDescriptors f2fmap_;
        boost::associative_property_map<MapBetweenFaceDescriptors> f2fmap(f2fmap_);
        //boost::associative_property_map<std::map<boost::graph_traits<Filtered_graph>::face_descriptor, face_descriptor>> f2fmap(f2fmap_);
        // EMesh3::Property_map<boost::graph_traits<Filtered_graph>::vertex_descriptor, vertex_descriptor> v2vmap = 
        //   tmesh.add_property_map<boost::graph_traits<Filtered_graph>::vertex_descriptor, vertex_descriptor>("v:v").first;
        // EMesh3::Property_map<boost::graph_traits<Filtered_graph>::face_descriptor, face_descriptor> f2fmap = 
        //   tmesh.add_property_map<boost::graph_traits<Filtered_graph>::face_descriptor, face_descriptor>("f:f").first;
        EMesh3 cmesh;
        CGAL::copy_face_graph(
          ffg, cmesh, CGAL::parameters::vertex_to_vertex_map(v2vmap).face_to_face_map(f2fmap)
        );

        copy_property<ffg_vertex_descriptor, vertex_descriptor, Rcpp::NumericVector>(tmesh, cmesh, v2vmap_, "v:normal");
        copy_property<ffg_vertex_descriptor, vertex_descriptor, std::string>(tmesh, cmesh, v2vmap_, "v:color");
        copy_property<ffg_vertex_descriptor, vertex_descriptor, double>(tmesh, cmesh, v2vmap_, "v:scalar");
        copy_property<ffg_face_descriptor, face_descriptor, std::string>(tmesh, cmesh, f2fmap_, "f:color");
        copy_property<ffg_face_descriptor, face_descriptor, double>(tmesh, cmesh, f2fmap_, "f:scalar");

        // reverse the map between vertices
        // std::map<vertex_descriptor, vertex_descriptor> rv2vmap;
        // for(const auto& [vertex_ffg, vertex_cmesh] : v2vmap_) {
        //   rv2vmap.insert(std::make_pair(
        //     vertex_cmesh, CGAL::SM_Vertex_index(int(vertex_ffg))
        //   ));
        // }

        // MapBetweenFaces fcomponent = fcomponents[c];

        // if(hasFcolors) {
        //   Fcolors_map cfcolor = 
        //     cmesh.add_property_map<face_descriptor, std::string>(
        //       "f:color", ""
        //     ).first;
        //   for(EMesh3::Face_index cf : cmesh.faces()) {
        //     face_descriptor fd = fcomponent[cf];
        //     cfcolor[cf] = fcolors[fd];
        //   }          
        // }
        // if(hasFscalars) {
        //   Fscalars_map cfscalar = 
        //     cmesh.add_property_map<face_descriptor, double>(
        //       "f:scalar", nan("")
        //     ).first;
        //   for(EMesh3::Face_index cf : cmesh.faces()) {
        //     face_descriptor fd = fcomponent[cf];
        //     cfscalar[cf] = fscalars[fd];
        //   }          
        // }

        // Normals_map cvnormal;
        // Vcolors_map cvcolor;
        // Vscalars_map cvscalar;
        // if(hasNormals) {
        //   cvnormal = 
        //     cmesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
        //       "v:normal", defaultNormal()
        //     ).first;
        // }
        // if(hasVcolors) {
        //   cvcolor = 
        //     cmesh.add_property_map<vertex_descriptor, std::string>(
        //       "v:color", ""
        //     ).first;
        // }
        // if(hasVscalars) {
        //   cvscalar = 
        //     cmesh.add_property_map<vertex_descriptor, double>(
        //       "v:scalar", nan("")
        //     ).first;
        // }

        // if(hasNormals || hasVcolors || hasVscalars) {
        //   for(EMesh3::Vertex_index cv : cmesh.vertices()) {
        //     face_descriptor cf = cmesh.face(cmesh.halfedge(cv));
        //     face_descriptor pf = fcomponent[cf]; // face index of cv in the parent mesh
        //     // std::vector<vertex_descriptor> face = faceVertices[pf];
        //     // int vindex = 0;
        //     // for(EMesh3::Vertex_index v : vertices_around_face(cmesh.halfedge(cf), cmesh)) {
        //     //   if(v == cv) {
        //     //     break;
        //     //   } else {
        //     //     vindex++;
        //     //   }
        //     // }
        //     // auto pvs = vertices_around_face(tmesh.halfedge(pf), tmesh).begin();
        //     // for(EMesh3::Vertex_index v : vertices_around_face(cmesh.halfedge(cf), cmesh)) {
        //     //   if(v == cv) {
        //     //     break;
        //     //   } else {
        //     //     pvs++;
        //     //   }
        //     // }
        //     // vertex_descriptor VVV = *pvs;
        //     // int vi = CGAL::vertex_index_in_face(cv, cf, cmesh);
        //     // halfedge_descriptor ph = tmesh.halfedge(pf);
        //     // halfedge_descriptor pph = tmesh.halfedge(pf);
        //     // //auto vs = vertices_around_target(tmesh.halfedge(pf), tmesh).begin();
        //     // halfedge_descriptor ch = cmesh.halfedge(cf);
        //     // halfedge_descriptor cch = cmesh.halfedge(cf);
        //     // vertex_descriptor vmaster = tmesh.source(ph);
        //     // vertex_descriptor vvmaster;
        //     // bool b = cmesh.source(ch) == cv;
        //     // int i = 0;
        //     // while(!b) {
        //     //   cch = cmesh.next(cch);
        //     //   pph = tmesh.next(pph);
        //     //   b = cmesh.source(cch) == cv;
        //     //   i++;
        //     //   if(i > 3) {
        //     //     Rcpp::stop("xxxxxxxxx");
        //     //   }
        //     // }
        //     // vvmaster = tmesh.source(pph);
        //     // int vvi = vi;
        //     // while(vi > 0) {
        //     //   //vs++;
        //     //   ch = cmesh.next(ch);
        //     //   vi--;
        //     //   vvi--;
        //     //   if(cmesh.source(ch) == cv) {
        //     //     vmaster = tmesh.source(ph);
        //     //     vi = 0;
        //     //   }
        //     //   ph = tmesh.next(ph);
        //     // }
        //     // vertex_descriptor ccv = cmesh.source(ch);
        //     // vertex_descriptor pv = tmesh.source(ph); // *vs;
        //     // vertex_descriptor ppv = tmesh.target(ph);
        //     EMesh3::Vertex_index pv = rv2vmap[cv]; // vertex index of cv in the parent mesh
        //     if(hasNormals) {
        //       cvnormal[cv] = vnormals[pv];
        //     }
        //     if(hasVcolors) {
        //       cvcolor[cv] = vcolors[pv];
        //     }
        //     if(hasVscalars) {
        //       cvscalar[cv] = vscalars[pv];
        //     }
        //     // for(const auto& [key, value] : v2vmap_) {
        //     //   if(value == cv) {
        //     //     cvnormal[cv] = vnormals[CGAL::SM_Vertex_index(int(key))];
        //     //     break;
        //     //   }
        //     // }
        //     if(c == 1 && int(cv) == 1) {
        //       // Rcpp::Rcout << "f2fmap_:\n";
        //       // for(const auto& [key, value] : f2fmap_) {
        //       //   Rcpp::Rcout << key << " - " << value << "\n";
        //       // }
        //       // Rcpp::Rcout << "f2fmap:\n";
        //       // for(EMesh3::Face_index cf : cmesh.faces()) {
        //       //   int k = int(fcomponents[c][cf]);
        //       //   boost::graph_traits<Filtered_graph>::face_descriptor fd(k);
        //       //   Rcpp::Rcout << fd << " - " << f2fmap[fd] << "\n";
        //       // }
        //       Rcpp::Rcout << "v2vmap_:\n";
        //       for(const auto& [key, value] : v2vmap_) {
        //         Rcpp::Rcout << key << " - " << value << "\n";
        //       }
        //       // Rcpp::Rcout << "vmasterprev: " << tmesh.point(tmesh.source(tmesh.prev(ph))) << "\n";
        //       // Rcpp::Rcout << "vvi: " << vvi << "\n";
        //       // Rcpp::Rcout << "vi: " << vi << "\n";
        //       // Rcpp::Rcout << "i: " << i << "\n";
        //       // Rcpp::Rcout << "vvmaster: " << tmesh.point(vvmaster) << "\n";
        //       // Rcpp::Rcout << "vmaster: " << tmesh.point(vmaster) << "\n";
        //       // Rcpp::Rcout << "cv: " << cmesh.point(cv) << "\n";
        //       // Rcpp::Rcout << "ccv: " << cmesh.point(ccv) << "\n";
        //       // Rcpp::Rcout << "pv: " << tmesh.point(pv) << "\n";
        //       // Rcpp::Rcout << "ppv: " << tmesh.point(ppv) << "\n";
        //       // Rcpp::Rcout << "halfedges: " << mesh.halfedge(pf) << " - " << tmesh.halfedge(pf) << "\n";
        //       for(EMesh3::Vertex_index v : vertices_around_face(cmesh.halfedge(cf), cmesh)) {
        //         Rcpp::Rcout << tmesh.point(rv2vmap[v]) << "\n";
        //       }
        //       // for(EMesh3::Vertex_index v : vertices_around_face(mesh.halfedge(pf), mesh)) {
        //       //   Rcpp::Rcout << v << "\n";
        //       //   Rcpp::Rcout << mesh.point(v) << "\n";
        //       // }
        //       // for(int j = 0; j < 3; j++) {
        //       //   Rcpp::Rcout << tmesh.point(face[j]) << "\n";
        //       // }
        //       for(EMesh3::Vertex_index v : vertices_around_face(tmesh.halfedge(pf), tmesh)) {
        //         Rcpp::Rcout << v << "\n";
        //         Rcpp::Rcout << tmesh.point(v) << "\n";
        //       }
        //     }
        //   }
        // }
        xptrs(c) = Rcpp::XPtr<EMesh3>(new EMesh3(cmesh), false);
      }
      return xptrs;
    }


    /*
    Otherwise we use the connected components decomposition made with 
    components of vertex indices.
    */
    Vertex_index_map vccmap = 
      tmesh.add_property_map<vertex_descriptor, std::size_t>("v:CC").first;
    const std::size_t ncc = connected_components(tmesh, vccmap);
    std::vector<std::vector<EMesh3::Vertex_index>> components(ncc);
    for(EMesh3::Vertex_index v : tmesh.vertices()) {
      std::size_t c = vccmap[v];
      components[c].push_back(v);
    }
    if(ncc == 1) {
      Message("Only one component found.\n");
    } else {
      const std::string msg = "Found " + std::to_string(ncc) + " components.\n";
      Message(msg);
    }
    
    Rcpp::List mymeshes(ncc);
    
    std::vector<EMesh3> ccmeshes0(ncc);
    std::map<EMesh3::Face_index, int> mapfaces;
    int iii = 0;
    for(EMesh3::Face_index f : tmesh.faces()) {
      mapfaces[f] = iii++;
    }
    Rcpp::Rcout << "mapfaces done\n";
    //
    std::vector<std::vector<std::vector<EMesh3::Vertex_index>>> connec(ncc);
    std::vector<std::vector<EMesh3::Face_index>> todelete(ncc);
    std::vector<std::size_t> nfaces(ncc);
    for(std::size_t c = 0; c < ncc; c++) {
      std::vector<EMesh3::Vertex_index> component = components[c];
      std::size_t csize = component.size();
      std::map<EMesh3::Vertex_index, int> newinds;
      for(int j = 0; j < csize; j++) {
        newinds[component[j]] = j;
      }
      Rcpp::Rcout << "newinds done\n";
      for(const auto& [faceindex, intindex] : mapfaces) {
        auto it = vertices_around_face(tmesh.halfedge(faceindex), tmesh);
        if(vccmap[*(it.begin())] == c) {
          std::vector<EMesh3::Vertex_index> facenewinds;
          for(EMesh3::Vertex_index v : it) {
            //Rcpp::Rcout << "newind: " << newinds[v] << "\n";
            facenewinds.push_back(CGAL::SM_Vertex_index(newinds[v]));
          }
          //Rcpp::Rcout << "facenewinds done\n";
          connec[c].push_back(facenewinds);
          todelete[c].push_back(faceindex);
        }
      }
      nfaces[c] = todelete[c].size();
      for(int fi = 0; fi < nfaces[c]; fi++) {
        mapfaces.erase(todelete[c][fi]);
        //Rcpp::Rcout << "facedindex erased";
      }
      EMesh3 dd;
      for(int j = 0; j < csize; j++) {
        EPoint3 pt = tmesh.point(component[j]);
        dd.add_vertex(pt);
      }
      Rcpp::Rcout << "vertices added to dd\n";
      ccmeshes0[c] = dd;
    }
    std::vector<EMesh3> ccmeshes(ncc);
    for(std::size_t c = 0; c < ncc; c++) {
      EMesh3 dd = ccmeshes0[c];
      for(int j = 0; j < nfaces[c]; j++) {
        dd.add_face(connec[c][j]);
      }
      Rcpp::Rcout << "faces added to dd\n";
      ccmeshes[c] = dd;

      Rcpp::Rcout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
    }

    if(hasNormals) {
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Vertex_index> component = components[c];
        Normals_map vnormals = ccmeshes[c]
          .add_property_map<vertex_descriptor, Rcpp::NumericVector>(
            "v:normal"
          ).first;
        for(EMesh3::Vertex_index vi : ccmeshes[c].vertices()) {
          vertex_descriptor vmaster = 
            CGAL::SM_Vertex_index(component[int(vi)]);
          vnormals[vi] = vnormals_.first[vmaster];
        }
      }
    }
    if(hasVcolors) {
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Vertex_index> component = components[c];
        Vcolors_map vcolors = ccmeshes[c]
          .add_property_map<vertex_descriptor, std::string>("v:color").first;
        for(EMesh3::Vertex_index vi : ccmeshes[c].vertices()) {
          vertex_descriptor vmaster = 
            CGAL::SM_Vertex_index(component[int(vi)]);
          vcolors[vi] = vcolors_.first[vmaster];
        }
      }
    }
    if(hasFcolors) {
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Face_index> fcomponent = todelete[c];
        Fcolors_map fcolors = ccmeshes[c]
          .add_property_map<face_descriptor, std::string>("f:color").first;
        for(EMesh3::Face_index fi : ccmeshes[c].faces()) {
          face_descriptor fmaster = 
            CGAL::SM_Face_index(fcomponent[int(fi)]);
          fcolors[fi] = fcolors_.first[fmaster];
        }
      }
    }

    for(std::size_t c = 0; c < ncc; c++) {
      mymeshes(c) = Rcpp::XPtr<EMesh3>(new EMesh3(ccmeshes[c]), false);
    }
    
    return mymeshes;
  }


  Rcpp::List convexParts(const bool triangulate) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    NefPol nef(mesh);
    CGAL::convex_decomposition_3(nef);
    std::list<EMesh3> convex_parts;
    // the first volume is the outer volume, ignored in the decomposition
    NefPol::Volume_const_iterator ci = ++nef.volumes_begin();
    for( ; ci != nef.volumes_end(); ++ci) {
      if(ci->mark()) {
        EPolyhedron pol;
        nef.convert_inner_shell_to_polyhedron(ci->shells_begin(), pol);
        EMesh3 cmesh;
        CGAL::copy_face_graph(pol, cmesh);
        convex_parts.push_back(cmesh);
      }
    }
    const size_t ncp = convex_parts.size();
    std::string msg;
    if(ncp == 1) {
      msg = "Only one convex part found.";
    } else {
      msg = "Found " + std::to_string(ncp) + " convex parts.";
    }
    Message(msg);
    Rcpp::List out(ncp);
    size_t i = 0;
    for(EMesh3 cmesh : convex_parts) {
      if(triangulate && !CGAL::is_triangle_mesh(cmesh)) {
        if(!PMP::triangulate_faces(cmesh)) {
          Rcpp::stop("Triangulation has failed.");
        }
      }
      out(i++) = Rcpp::XPtr<EMesh3>(new EMesh3(cmesh), false);
    }
    return out;
  }
  
  
  Rcpp::NumericVector distance(Rcpp::NumericMatrix points) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    const size_t npoints = points.ncol();
    Rcpp::NumericVector distances(npoints);
    for(size_t i = 0; i < npoints; i++){
      Rcpp::NumericVector point_i = points(Rcpp::_, i);
      std::vector<EPoint3> pt = {EPoint3(point_i(0), point_i(1), point_i(2))};
      distances(i) = 
        PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(
          pt, mesh
        );
    }
    return distances;
  }  


  bool doesBoundVolume() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    return PMP::does_bound_a_volume(mesh);
  }


  bool doesSelfIntersect() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    return PMP::does_self_intersect(mesh);
  }


  Rcpp::XPtr<EMesh3> dual() {
    EMesh3 dualmesh = dualMesh(mesh);
    return Rcpp::XPtr<EMesh3>(new EMesh3(dualmesh), false);
  }


  Rcpp::DataFrame edges() {
    return getEdges<EK, EMesh3, EPoint3>(mesh);
  }


  Rcpp::IntegerVector facesAroundVertex(int v) {
    if(v >= mesh.number_of_vertices()) {
      Rcpp::stop("Too high vertex index.");
    }
    vertex_descriptor vd = CGAL::SM_Vertex_index(v);
    // Rcpp::Rcout << "degree: " << mesh.degree(vd) << "\n";
    // Rcpp::Rcout << "halfedge: " << mesh.halfedge(vd) << "\n";
    Rcpp::Rcout << "mesh.halfedge(vd): " << mesh.halfedge(vd) << "halfedge(vd, mesh): " << halfedge(vd, mesh) << "\n";
    Rcpp::IntegerVector faces;
    for(face_descriptor fd : CGAL::faces_around_target(mesh.halfedge(vd), mesh)) {
      //Rcpp::Rcout << "face is removed: " << mesh.is_removed(fd) << "\n";
      if(fd != EMesh3::null_face()) {
        faces.push_back(int(fd) + 1);
      }
      // Rcpp::Rcout << "face: " << fd << "\n";
      // Rcpp::Rcout << "isnull: " << (fd == EMesh3::null_face()) << "\n";
      // Rcpp::Rcout << "asint: " << int(fd) << "\n";
    }
    return faces;
    // for(halfedge_descriptor hd : CGAL::halfedges_around_target(vd, mesh)) {
    //   Rcpp::Rcout << "halfedge around target: " << hd << "\n";
    // }
    // for(halfedge_descriptor hd : CGAL::halfedges_around_source(vd, mesh)) {
    //   Rcpp::Rcout << "halfedge around source: " << hd << "\n";
    // }
  }


  void fair(Rcpp::IntegerVector indices) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    std::list<vertex_descriptor> selectedVertices;
    const int nindices = indices.size();
    const int nvertices = mesh.number_of_vertices();
    for(int i = 0; i < nindices; i++) {
      const int idx = indices(i);
      if(idx >= nvertices) {
        Rcpp::stop("Too large index.");
      }
      selectedVertices.push_back(*(mesh.vertices().begin() + idx));
    }
    removeProperties(mesh, {"v:normal"});
    const bool success = PMP::fair(mesh, selectedVertices);
    if(!success) {
      Rcpp::stop("Failed to fair the mesh.");
    }
  }


  Rcpp::List filterMesh(Rcpp::IntegerVector selectedFaces) {
    Face_index_map fimap = 
      mesh.add_property_map<face_descriptor, std::size_t>("f:i", 2).first;
    const int nfaces = mesh.number_of_faces();
    for(int i = 0; i < selectedFaces.size(); i++) {
      const int idx = selectedFaces(i);
      if(idx >= nfaces) {
        Rcpp::stop("Too large face index.");
      }
      fimap[CGAL::SM_Face_index(idx)] = 1;
    }
    EMesh3 fmesh1;
    {
      Filtered_graph ffg(mesh, 1, fimap);
      MapBetweenVertexDescriptors v2vmap_;
      boost::associative_property_map<MapBetweenVertexDescriptors> 
        v2vmap(v2vmap_);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> 
        f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, fmesh1, 
        CGAL::parameters::vertex_to_vertex_map(v2vmap).face_to_face_map(f2fmap)
      );
      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, Rcpp::NumericVector
      >(mesh, fmesh1, v2vmap_, "v:normal");
      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, std::string
      >(mesh, fmesh1, v2vmap_, "v:color");
      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, double
      >(mesh, fmesh1, v2vmap_, "v:scalar");
      copy_property<
        ffg_face_descriptor, face_descriptor, std::string
      >(mesh, fmesh1, f2fmap_, "f:color");
      copy_property<
        ffg_face_descriptor, face_descriptor, double
      >(mesh, fmesh1, f2fmap_, "f:scalar");
    }
    EMesh3 fmesh2;
    {
      Filtered_graph ffg(mesh, 2, fimap);
      MapBetweenVertexDescriptors v2vmap_;
      boost::associative_property_map<MapBetweenVertexDescriptors> 
        v2vmap(v2vmap_);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> 
        f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, fmesh2, 
        CGAL::parameters::vertex_to_vertex_map(v2vmap).face_to_face_map(f2fmap)
      );
      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, Rcpp::NumericVector
      >(mesh, fmesh2, v2vmap_, "v:normal");
      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, std::string
      >(mesh, fmesh2, v2vmap_, "v:color");
      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, double
      >(mesh, fmesh2, v2vmap_, "v:scalar");
      copy_property<
        ffg_face_descriptor, face_descriptor, std::string
      >(mesh, fmesh2, f2fmap_, "f:color");
      copy_property<
        ffg_face_descriptor, face_descriptor, double
      >(mesh, fmesh2, f2fmap_, "f:scalar");
    }
    mesh.remove_property_map(fimap);
    return Rcpp::List::create(
      Rcpp::Named("fmesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(fmesh1), false),
      Rcpp::Named("fmesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(fmesh2), false)
    );
  }


//   Rcpp::IntegerVector filterGraph() {
//     trivial_edge_predicate filter_edge;
//     positive_vertex_scalar filter_vertex(mesh.property_map<vertex_descriptor, double>("v:scalar").first);
//     boost::filtered_graph<EMesh3, trivial_edge_predicate, positive_vertex_scalar> fg(mesh, filter_edge, filter_vertex);
//     // EMesh3 submesh;
//     // CGAL::copy_face_graph(fg, submesh);
//     // return Rcpp::XPtr<EMesh3>(new EMesh3(submesh), false);
//     boost::property_map<boost::filtered_graph<EMesh3, trivial_edge_predicate, positive_vertex_scalar>, CGAL::vertex_index_t>::type pmap = 
//       get(CGAL::vertex_index, fg);
//     std::pair<boost::filtered_graph<EMesh3, trivial_edge_predicate, positive_vertex_scalar>::vertex_iterator, boost::filtered_graph<EMesh3, trivial_edge_predicate, positive_vertex_scalar>::vertex_iterator> vv = boost::vertices(fg);
// //    std::vector<bool> keep_vertex(mesh.number_of_vertices(), false);
//     for(auto vd : boost::make_iterator_range(vv)) {
//       Rcpp::Rcout << vd << " - " << pmap[vd] << "\n";
// //      keep_vertex[int(pmap[vd])] = true;
//       //remove_vertex(vd, mesh);
//       //Rcpp::Rcout << "face of removed vertex: " << fd << "\n";
//     }
//     boost::property_map<boost::filtered_graph<EMesh3, trivial_edge_predicate, positive_vertex_scalar>, CGAL::edge_index_t>::type epmap = 
//       get(CGAL::edge_index, fg);
//     std::pair<boost::filtered_graph<EMesh3, trivial_edge_predicate, positive_vertex_scalar>::edge_iterator, boost::filtered_graph<EMesh3, trivial_edge_predicate, positive_vertex_scalar>::edge_iterator> ee = boost::edges(fg);
//     Rcpp::IntegerVector Edges;
//     int nedges = 0;
//     for(auto ed : boost::make_iterator_range(ee)) {
//       Rcpp::Rcout << ed << " - " << epmap[ed] << "\n";
//       Rcpp::Rcout << "source: " << boost::source(ed, fg) << " - target: " << boost::target(ed, fg) << "\n";
//       Edges.push_back(int(boost::source(ed, fg))); Edges.push_back(int(boost::target(ed, fg)));
//       nedges++;
//     }
//     Edges.attr("dim") = Rcpp::Dimension(2, nedges);
//     return Edges;
//     // Rcpp::Rcout << "mesh has garbage: " << mesh.has_garbage() << "\n";
//     // mesh.collect_garbage();
//     // int nvertices = mesh.number_of_vertices();
//     // Face_index_map fimap = 
//     //   mesh.add_property_map<face_descriptor, std::size_t>("f:i", 0).first;
//     // EMesh3::Face_range faces = mesh.faces();
//     // for(face_descriptor fd : faces) {
//     //   bool include = true;
//     //   for(vertex_descriptor vd : vertices_around_face(mesh.halfedge(fd), mesh)) {
//     //     if(!keep_vertex[int(vd)]) {
//     //       include = false;
//     //       fimap[fd] = 1;
//     //       break;
//     //     }
//     //     // if(int(vd) >= nvertices) {
//     //     //   //remove_face(fd, mesh);
//     //     //   //Rcpp::Rcout << "face of removed face: " << f << "\n";
//     //     //   break;
//     //     // }
//     //   }
//     // }
//     // Rcpp::Rcout << "mesh has garbage: " << mesh.has_garbage() << "\n";
//     // mesh.collect_garbage();
//     // Filtered_graph ffg(mesh, 0, fimap);
//     // EMesh3 meshcopy;
//     // CGAL::copy_face_graph(ffg, meshcopy);
//     // EMesh3 submesh;
//     // boost::copy_graph(fg, submesh);
//     //return Rcpp::XPtr<EMesh3>(new EMesh3(meshcopy), false);
//     // Rcpp::Rcout << "vv.second:\n";
//     // for(vertex_descriptor v : vv.second) {
//     //   Rcpp::Rcout << pmap[v];
//     // }
//   }


  void fixManifoldness() {
    std::size_t nmv = PMP::duplicate_non_manifold_vertices(mesh);
    if(nmv > 0) {
      std::string msg;
      if(nmv = 1) {
        msg = "One non-manifold vertex duplicated.";
      } else {
        msg = "Duplicated " + std::to_string(nmv) + " non-manifold vertices.";
      }
      Message(msg);
    }
  }


  Rcpp::NumericVector geoDists(const int index) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    const int nvertices = mesh.number_of_vertices();
    if(index >= nvertices) {
      Rcpp::stop("Too large index.");
    }
    // property map for the distance values to the source set
    Vertex_distance_map vertex_distance = 
      mesh.add_property_map<vertex_descriptor, double>("v:distance", 0).first;
    vertex_descriptor source = *(std::next(mesh.vertices().begin(), index));
    CGAL::Heat_method_3::estimate_geodesic_distances(
      mesh, vertex_distance, source
    );
    Rcpp::NumericVector gdistances(nvertices);
    int i = 0;
    for(vertex_descriptor vd : mesh.vertices()) {
      gdistances(i++) = get(vertex_distance, vd);
    }
    mesh.remove_property_map(vertex_distance);
    return gdistances;    
  }


  Rcpp::List getBorders() {
    std::vector<halfedge_descriptor> border_cycles;
    PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
    int nborders = border_cycles.size();
    if(nborders == 0) {
      return Rcpp::List::create();
    }
    Rcpp::List Borders(nborders);
    Rcpp::CharacterVector colnames = {"edge", "v1", "v2"};
    for(int i = 0; i < nborders; i++) {
      Rcpp::IntegerVector border_i;
      halfedge_descriptor h = border_cycles[i];
      int nedges = 0;
      for(halfedge_descriptor hd : halfedges_around_face(h, mesh)) {
        edge_descriptor ed = mesh.edge(hd);
        vertex_descriptor v1 = source(ed, mesh);
        vertex_descriptor v2 = target(ed, mesh);
        border_i.push_back(int(ed) + 1);
        border_i.push_back(int(v1) + 1);
        border_i.push_back(int(v2) + 1);
        nedges++;
      }
      border_i.attr("dim") = Rcpp::Dimension(3, nedges);
      Rcpp::IntegerMatrix Border_i = 
        Rcpp::transpose(Rcpp::as<Rcpp::IntegerMatrix>(border_i));
      Rcpp::colnames(Border_i) = colnames;
      Borders(i) = Border_i;
    }
    return Borders;
  }


  Rcpp::IntegerMatrix getFacesMatrix() {
    const size_t nfaces = mesh.number_of_faces();
    if(CGAL::is_triangle_mesh(mesh)) {
      Rcpp::IntegerMatrix Faces(3, nfaces);
      {
        int i = 0;
        for(EMesh3::Face_index fi : mesh.faces()) {
          auto vs = vertices_around_face(mesh.halfedge(fi), mesh).begin();
          Rcpp::IntegerVector col_i = 
            {int(*(vs++)) + 1, int(*(vs++)) + 1, int(*vs) + 1};
          Faces(Rcpp::_, i++) = col_i;
        }
      }
      return Rcpp::transpose(Faces);
    } else if(CGAL::is_quad_mesh(mesh)) {
      Rcpp::IntegerMatrix Faces(4, nfaces);
      {
        int i = 0;
        for(EMesh3::Face_index fi : mesh.faces()) {
          auto vs = vertices_around_face(mesh.halfedge(fi), mesh).begin();
          Rcpp::IntegerVector col_i = 
            {int(*(vs++)) + 1, int(*(vs++)) + 1, int(*(vs++)) + 1, int(*vs) + 1};
          Faces(Rcpp::_, i++) = col_i;
        }
      }
      return Rcpp::transpose(Faces);
    } else {
      Rcpp::stop("This function can be used with triangle or quad meshes only.");
    }
  }


  Rcpp::List getFacesList() {
    return getFaces<EMesh3>(mesh);
  }


  Rcpp::Nullable<Rcpp::StringVector> getFcolors() {
    std::pair<Fcolors_map, bool> fcolorsmap_ = 
      mesh.property_map<face_descriptor, std::string>("f:color");
    if(!fcolorsmap_.second) {
      return R_NilValue;
    }
    Rcpp::StringVector Fcolors(mesh.number_of_faces());
    int i = 0;
    for(EMesh3::Face_index fi : mesh.faces()) {
      Fcolors(i++) = fcolorsmap_.first[fi];
    }
    return Rcpp::Nullable<Rcpp::StringVector>(Fcolors);
  }


  Rcpp::Nullable<Rcpp::NumericVector> getFscalars() {
    std::pair<Fscalars_map, bool> fscalarsmap_ = 
      mesh.property_map<face_descriptor, double>("f:scalar");
    if(!fscalarsmap_.second) {
      return R_NilValue;
    }
    Rcpp::NumericVector Fscalars(mesh.number_of_faces());
    int i = 0;
    for(EMesh3::Face_index fi : mesh.faces()) {
      Fscalars(i++) = fscalarsmap_.first[fi];
    }
    return Rcpp::Nullable<Rcpp::NumericVector>(Fscalars);
  }


  Rcpp::IntegerMatrix getHalfedges() {
    Rcpp::IntegerMatrix Halfedges(mesh.number_of_halfedges(), 6);
    int i = 0;
    for(halfedge_descriptor hd : mesh.halfedges()) {
      Rcpp::IntegerVector row = {int(hd), int(mesh.source(hd)), int(mesh.target(hd)), int(mesh.opposite(hd)), int(mesh.face(hd)), int(mesh.is_border(hd))};
      Halfedges(i++, Rcpp::_) = row;
    }
    return Halfedges;
  }  


  Rcpp::Nullable<Rcpp::NumericMatrix> getNormals() {
    std::pair<Normals_map, bool> normalsmap_ = 
      mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
    if(normalsmap_.second) {
      Normals_map normalsmap = normalsmap_.first;
      Rcpp::NumericMatrix Normals(3, mesh.number_of_vertices());
      for(int i = 0; i < mesh.number_of_vertices(); i++) {
        Normals(Rcpp::_, i) = normalsmap[CGAL::SM_Vertex_index(i)];
      }
      return Rcpp::Nullable<Rcpp::NumericMatrix>(Rcpp::transpose(Normals));
    }
    return R_NilValue;
  }


  Rcpp::Nullable<Rcpp::StringVector> getVcolors() {
    std::pair<Vcolors_map, bool> vcolorsmap_ = 
      mesh.property_map<vertex_descriptor, std::string>("v:color");
    if(!vcolorsmap_.second) {
      return R_NilValue;
    }
    Rcpp::StringVector Vcolors(mesh.number_of_vertices());
    int i = 0;
    for(EMesh3::Vertex_index vi : mesh.vertices()) {
      Vcolors(i++) = vcolorsmap_.first[vi];
    }
    return Rcpp::Nullable<Rcpp::StringVector>(Vcolors);
  }


  Rcpp::Nullable<Rcpp::NumericVector> getVscalars() {
    std::pair<Vscalars_map, bool> vscalarsmap_ = 
      mesh.property_map<vertex_descriptor, double>("v:scalar");
    if(!vscalarsmap_.second) {
      return R_NilValue;
    }
    Rcpp::NumericVector Vscalars(mesh.number_of_vertices());
    int i = 0;
    for(EMesh3::Vertex_index vi : mesh.vertices()) {
      Vscalars(i++) = vscalarsmap_.first[vi];
    }
    return Rcpp::Nullable<Rcpp::NumericVector>(Vscalars);
  }


  Rcpp::NumericMatrix getVertices() {
    return getVertices_EK(mesh);
  }  


  Rcpp::List getRmesh() {
    std::pair<Normals_map, bool> normalsmap_ = 
      mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
    const bool there_is_normals = normalsmap_.second;
    Rcpp::List rmesh;
    if(CGAL::is_triangle_mesh(mesh)) {
      rmesh = RSurfEKMesh2(mesh, false, 3);
    } else if(CGAL::is_quad_mesh(mesh)) {
      rmesh = RSurfEKMesh2(mesh, false, 4);
    } else {
      rmesh = RSurfEKMesh(mesh, false);
    }
    if(there_is_normals) {
      Normals_map normalsmap = normalsmap_.first;
      Rcpp::NumericMatrix Normals(3, mesh.number_of_vertices());
      for(int i = 0; i < mesh.number_of_vertices(); i++) {
        Normals(Rcpp::_, i) = normalsmap[CGAL::SM_Vertex_index(i)];
      }
      rmesh["normals"] = Normals;
    }
    std::pair<Vcolors_map, bool> vcolorsmap_ = 
      mesh.property_map<vertex_descriptor, std::string>("v:color");
    std::pair<Fcolors_map, bool> fcolorsmap_ = 
      mesh.property_map<face_descriptor, std::string>("f:color");
    if(vcolorsmap_.second || fcolorsmap_.second) {
      if(vcolorsmap_.second) {
        Rcpp::StringVector Colors(mesh.number_of_vertices());
        int i = 0;
        for(EMesh3::Vertex_index v : mesh.vertices()) {
          Colors(i++) = vcolorsmap_.first[v];
        }
        rmesh["colors"] = Colors;
      } else {
        Rcpp::StringVector Colors(mesh.number_of_faces());
        int i = 0;
        for(EMesh3::Face_index f : mesh.faces()) {
          Colors(i++) = fcolorsmap_.first[f];
        }
        rmesh["colors"] = Colors;
      }
    }
    return rmesh;
  }


  Rcpp::XPtr<EMesh3> intersection(Rcpp::XPtr<EMesh3> mesh2XPtr) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The reference mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The reference mesh self-intersects.");
    }
    EMesh3 mesh2 = *(mesh2XPtr.get());
    if(!CGAL::is_triangle_mesh(mesh2)) {
      Rcpp::stop("The second mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh2)) {
      Rcpp::stop("The second mesh self-intersects.");
    }
    EMesh3 imesh;
    const bool success = PMP::corefine_and_compute_intersection(
      mesh, mesh2, imesh
    );
    if(!success) {
      Rcpp::stop("Intersection computation has failed.");
    }
    return Rcpp::XPtr<EMesh3>(new EMesh3(imesh), false);
  }


  bool isClosed() {
    return CGAL::is_closed(mesh);
  }


  bool isOutwardOriented() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    return PMP::is_outward_oriented(mesh);
  }


  bool isQuad() {
    return CGAL::is_quad_mesh(mesh);
  }


  bool isTriangle() {
    return CGAL::is_triangle_mesh(mesh);
  }


  bool isValid() {
    return mesh.is_valid();
  }


  bool isValidFaceGraph() {
    return CGAL::is_valid_face_graph(mesh, true);
  }


  bool isValidHalfedgeGraph() {
    return CGAL::is_valid_halfedge_graph(mesh, true);
  }


  bool isValidPolygonMesh() {
    return CGAL::is_valid_polygon_mesh(mesh, true);
  }


  void orientToBoundVolume() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    PMP::orient_to_bound_a_volume(mesh);
    // faut-il updater normals?
  }


  void print() {
    Rcpp::Rcout << "Mesh with " << mesh.number_of_vertices() 
                << " vertices and " << mesh.number_of_faces() << " faces.\n";
  }


  void removeSelfIntersections() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    PMP::experimental::remove_self_intersections(mesh);
    mesh.collect_garbage();
  }


  void reverseFaceOrientations() {
    PMP::reverse_face_orientations(mesh);
    // update normals
  }


  Rcpp::XPtr<EMesh3> subtract(Rcpp::XPtr<EMesh3> mesh2XPtr) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The reference mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The reference mesh self-intersects.");
    }
    EMesh3 mesh2 = *(mesh2XPtr.get());
    if(!CGAL::is_triangle_mesh(mesh2)) {
      Rcpp::stop("The second mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh2)) {
      Rcpp::stop("The second mesh self-intersects.");
    }
    EMesh3 imesh;
    const bool success = PMP::corefine_and_compute_difference(
      mesh, mesh2, imesh
    );
    if(!success) {
      Rcpp::stop("Difference computation has failed.");
    }
    return Rcpp::XPtr<EMesh3>(new EMesh3(imesh), false);
  }


  void triangulate() {
    triangulateMesh(mesh);
  }


  Rcpp::List Union(Rcpp::XPtr<EMesh3> mesh2XPtr) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The reference mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The reference mesh self-intersects.");
    }
    EMesh3 mesh2 = *(mesh2XPtr.get());
    if(!CGAL::is_triangle_mesh(mesh2)) {
      Rcpp::stop("The second mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh2)) {
      Rcpp::stop("The second mesh self-intersects.");
    }

    EMesh3 mcopy;
    CGAL::copy_face_graph(mesh, mcopy);

    int nfaces1 = mesh.number_of_faces();
    int nfaces2 = mesh2.number_of_faces();

    std::pair<Fcolors_map, bool> fcolor1_ = 
      mesh.property_map<face_descriptor, std::string>("f:color");
    std::pair<Fcolors_map, bool> fcolor2_ = 
      mesh2.property_map<face_descriptor, std::string>("f:color");
    Fcolors_map fcolor1; 
    Fcolors_map fcolor2;
    const bool hasColors = fcolor1_.second && fcolor2_.second;
    if(hasColors) {
      fcolor1 = fcolor1_.first;
      fcolor2 = fcolor2_.first;
    }

    std::pair<Fscalars_map, bool> fscalar1_ = 
      mesh.property_map<face_descriptor, double>("f:scalar");
    std::pair<Fscalars_map, bool> fscalar2_ = 
      mesh2.property_map<face_descriptor, double>("f:scalar");
    Fscalars_map fscalar1; 
    Fscalars_map fscalar2;
    const bool hasScalars = fscalar1_.second && fscalar2_.second;
    if(hasScalars) {
      fscalar1 = fscalar1_.first;
      fscalar2 = fscalar2_.first;
    }

    // MaybeFcolorMap fcolor2_ = copy_fcolor(mesh2);
    // std::map<face_descriptor, std::string> fcolor2 = fcolor2_.first;
    // for(const auto& [fd, col] : fcolor2) {
    //   Rcpp::Rcout << fd << "->" << col << "\n";
    // }


    UnionVisitor vis;
    // Face_index_map fimap1 = 
    //   mesh.add_property_map<face_descriptor, std::size_t>("f:i").first;
    // Face_index_map fimap2 = 
    //   mesh2.add_property_map<face_descriptor, std::size_t>("f:i").first;


    EMesh3 umesh;
    const bool success = PMP::corefine_and_compute_union(
      mesh, mesh2, umesh,
      PMP::parameters::visitor(vis)
    );
    if(!success) {
      Rcpp::stop("Union computation has failed.");
    }

    Rcpp::Rcout << "nfaces mesh2: " << mesh2.number_of_faces() << "\n";

    Rcpp::Rcout << "nfaces mesh: " << mesh.number_of_faces() << "\n";
    Rcpp::Rcout << "mesh has garbage: " << mesh.has_garbage() << "\n";

    // if(true) {
    //   //Rcpp::StringVector fcolors0(fcolors);
    //   const int nfaces1 = mcopy.number_of_faces();
    //   std::vector<Triangle> triangles1;
    //   triangles1.reserve(nfaces1);
    //   for(EMesh3::Face_index fd : mcopy.faces()) {
    //     auto vd = vertices_around_face(mcopy.halfedge(fd), mcopy).begin();
    //     triangles1.emplace_back(Triangle(
    //       mcopy.point(*(++vd)), mcopy.point(*(++vd)), mcopy.point(*vd)
    //     ));
    //   }
    //   const int nfaces2 = mesh2.number_of_faces();
    //   std::vector<Triangle> triangles2;
    //   triangles2.reserve(nfaces2);
    //   for(EMesh3::Face_index fd : mesh2.faces()) {
    //     auto vd = vertices_around_face(mesh2.halfedge(fd), mesh2).begin();
    //     triangles2.emplace_back(Triangle(
    //       mesh2.point(*(++vd)), mesh2.point(*(++vd)), mesh2.point(*vd)
    //     ));
    //   }
    
    //   const size_t nf = imesh.number_of_faces();
    //   //Rcpp::StringVector fcolors1(nf);
    //   //size_t j = 0;
    //   for(EMesh3::Face_index f : imesh.faces()) {
    //     auto vd = vertices_around_face(imesh.halfedge(f), imesh).begin();
    //     Triangle tr(imesh.point(*(++vd)), imesh.point(*(++vd)), imesh.point(*vd));
    //     EPoint3 c = CGAL::centroid(tr);
    //     int k1, k2;
    //     bool found1 = false, found2 = false;
    //     for(k1 = 0; k1 < nfaces1; k1++) {
    //       if(triangles1[k1].has_on(c)) {
    //         found1 = true;
    //         break;
    //       }
    //     }
    //     if(!found1) {
    //       k1 = -1;
    //     }
    //     for(k2 = 0; k2 < nfaces2; k2++) {
    //       if(triangles2[k2].has_on(c)) {
    //         found2 = true;
    //         break;
    //       }
    //     }
    //     if(!found2) {
    //       k2 = -1;
    //     }
    //     Rcpp::Rcout << f << ": " << k1 << ", " << k2 << "\n";

    //     //fcolors1(j++) = fcolors0(k);
    //   }
    //   //fcolors = Rcpp::Nullable<Rcpp::StringVector>(fcolors1);
    // }

    Rcpp::Rcout << "nfaces: " << umesh.number_of_faces() << "\n";

    Rcpp::Rcout << "fprev: " << *(vis.fprev) << "\n";

    Rcpp::Rcout << "PROPERTIES" << "\n";
    std::vector<std::string> props = umesh.properties<face_descriptor>();
    for(int p = 0; p < props.size(); p++) {
      Rcpp::Rcout << props[p] << "\n";
    }

    // Face_index_map fimap = mesh.property_map<face_descriptor, std::size_t>("f:i").first;
    // Rcpp::IntegerVector Fimap(mesh.number_of_faces());
    // int i = 0;
    // for(EMesh3::Face_index fi : mesh.faces()) {
    //   Fimap(i++) = fimap[fi];
    // }


    // std::vector<bool> remove_on_mesh1(mesh.number_of_faces(), true);
    // std::vector<bool> remove_on_mesh2(mesh2.number_of_faces(), true);
    // std::vector<halfedge_descriptor> halfedges_to_remove;
    std::vector<std::string> fcolor_mesh1;
    std::vector<std::string> fcolor_mesh2;
    MapBetweenFaces fmap_mesh1 = *(vis.fmap_mesh1);
    MapBetweenFaces fmap_mesh2 = *(vis.fmap_mesh2);
    MapBetweenFaces fmap_union = *(vis.fmap_union);
    const bool disjoint = fmap_mesh1.size() == 0 && fmap_mesh2.size() == 0;
    Rcpp::Rcout << "disjoint: " << disjoint << "\n";


    if(disjoint) {
      std::map<vertex_descriptor, vertex_descriptor> vmap_union = *(vis.vmap_union);
      // vertex normals
      std::pair<Normals_map, bool> vnormal1_ = 
        mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
      std::pair<Normals_map, bool> vnormal2_ = 
        mesh2.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
      if(vnormal1_.second && vnormal2_.second) {
        Normals_map vnormal1 = vnormal1_.first;
        Normals_map vnormal2 = vnormal2_.first;
        Normals_map vnormal = umesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
          "v:normal", defaultNormal()
        ).first;
        for(int i = 0; i < mesh.number_of_vertices(); i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vnormal[vi] = vnormal1[vd];
        }
        for(int i = mesh.number_of_vertices(); i < umesh.number_of_vertices(); i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vnormal[vi] = vnormal2[vd];
        }
      }
      // vertex colors
      std::pair<Vcolors_map, bool> vcolor1_ = 
        mesh.property_map<vertex_descriptor, std::string>("v:color");
      std::pair<Vcolors_map, bool> vcolor2_ = 
        mesh2.property_map<vertex_descriptor, std::string>("v:color");
      if(vcolor1_.second && vcolor2_.second) {
        Vcolors_map vcolor1 = vcolor1_.first;
        Vcolors_map vcolor2 = vcolor2_.first;
        Vcolors_map vcolor = umesh.add_property_map<vertex_descriptor, std::string>(
          "v:color", ""
        ).first;
        for(int i = 0; i < mesh.number_of_vertices(); i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vcolor[vi] = vcolor1[vd];
        }
        for(int i = mesh.number_of_vertices(); i < umesh.number_of_vertices(); i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vcolor[vi] = vcolor2[vd];
        }
      }
      // vertex scalars
      std::pair<Vscalars_map, bool> vscalar1_ = 
        mesh.property_map<vertex_descriptor, double>("v:scalar");
      std::pair<Vscalars_map, bool> vscalar2_ = 
        mesh2.property_map<vertex_descriptor, double>("v:scalar");
      if(vscalar1_.second && vscalar2_.second) {
        Vscalars_map vscalar1 = vscalar1_.first;
        Vscalars_map vscalar2 = vscalar2_.first;
        Vscalars_map vscalar = umesh.add_property_map<vertex_descriptor, double>(
          "v:scalar", nan("")
        ).first;
        for(int i = 0; i < mesh.number_of_vertices(); i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vscalar[vi] = vscalar1[vd];
        }
        for(int i = mesh.number_of_vertices(); i < umesh.number_of_vertices(); i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vscalar[vi] = vscalar2[vd];
        }
      }
    }

    if(disjoint && (!hasColors && !hasScalars)) {
      return Rcpp::List::create(
        Rcpp::Named("umesh") = Rcpp::XPtr<EMesh3>(new EMesh3(umesh), false)
      );
    }

    Rcpp::Rcout << "size fmap_union: " << fmap_union.size() << "\n";
    int nfaces_umesh1 = *(vis.nfaces_umesh1);
    Rcpp::Rcout << "nfaces_umesh1: " << nfaces_umesh1 << "\n";
    if(disjoint) {
      nfaces_umesh1 = nfaces1;
    }
    Fcolors_map fcolor;
    Fscalars_map fscalar;
    if(hasColors) {
      fcolor = umesh.add_property_map<face_descriptor, std::string>(
        "f:color", ""
      ).first;
    }
    if(hasScalars) {
      fscalar = umesh.add_property_map<face_descriptor, double>(
        "f:scalar", nan("")
      ).first;
    }

    // bool is_mesh1 = true;
    Face_index_map fwhich = umesh.add_property_map<face_descriptor, std::size_t>(
      "f:which", 0
    ).first;
    // std::string color;
    for(int i = 0; i < nfaces_umesh1; i++) {
      face_descriptor fi = CGAL::SM_Face_index(i);
      face_descriptor fd = fmap_union[fi];
      face_descriptor fd1 = int(fd) < nfaces1 ? fd : fmap_mesh1[fd];
      if(hasColors) {
        fcolor[fi] = fcolor1[fd1];
      }
      if(hasScalars) {
        fscalar[fi] = fscalar1[fd1];
      }
      fwhich[fi] = 1;
    }
    for(int i = nfaces_umesh1; i < umesh.number_of_faces(); i++) {
      face_descriptor fi = CGAL::SM_Face_index(i);
      face_descriptor fd = fmap_union[fi];
      face_descriptor fd2 = int(fd) < nfaces2 ? fd : fmap_mesh2[fd];
      if(hasColors) {
        fcolor[fi] = fcolor2[fd2];
      }
      if(hasScalars) {
        fscalar[fi] = fscalar2[fd2];
      }
      fwhich[fi] = 2;
    }

    if(disjoint) {
      return Rcpp::List::create(
        Rcpp::Named("umesh") = Rcpp::XPtr<EMesh3>(new EMesh3(umesh), false)
      );
    }
/*     for(EMesh3::Face_index fi : umesh.faces()) {
      face_descriptor fd = fmap_union[fi];
      if(is_mesh1) {
        is_mesh1 = int(fi) <= int(fd);
        face_descriptor fd1 = int(fd) < nfaces1 ? fd : fmap_mesh1[fd];
        if(hasColors) {
          if(is_mesh1) {
            color = fcolor1[fd1];
            fcolor_mesh1.push_back(color);
          } else {
            color = fcolor2[fd];
            fcolor_mesh2.push_back(color);
          }
          fcolor[fi] = color;
        }
        if(hasScalars) {
          fscalar[fi] = is_mesh1 ? fscalar1[fd1] : fscalar2[fd];
        }
        if(is_mesh1) {
          fwhich[fi] = 1;
          //remove_on_mesh1[int(fd)] = false;
        } else {
          fwhich[fi] = 2;
          //remove_on_mesh2[int(fd)] = false;
        }
        // fcolor[fi] = int(fd) < nfaces1 ? 
        //   fcolor1[fi] : 
        //   is_mesh1 ? fcolor1[fmap_mesh1[fd]] : fcolor2[fd];
        // if(int(fd) < nfaces1) {
        //   fcolor[fi] = fcolor1[fi];
        // } else {
        //   fcolor[fi] = is_mesh1 ? fcolor1[fmap_mesh1[fd]] : fcolor2[fd];  
        // }
      } else {
        face_descriptor fd2 = int(fd) < nfaces2 ? fd : fmap_mesh2[fd];
        if(hasColors) {
          color = fcolor2[fd2];
          fcolor[fi] = color; 
          fcolor_mesh2.push_back(color);
        }
        if(hasScalars) {
          fscalar[fi] = fscalar2[fd2];
        }
        //remove_on_mesh2[int(fd)] = false;
        fwhich[fi] = 2;
      }
    }
 */
    // nfaces1 = mesh.number_of_faces();
    // nfaces2 = mesh2.number_of_faces();
    // for(int i = nfaces1-1; i >= 0; i--) {
    //   if(remove_on_mesh1[i]) {
    //     CGAL::Euler::remove_face(mesh.halfedge(CGAL::SM_Face_index(i)), mesh);
    //     //mesh.remove_face(CGAL::SM_Face_index(i));
    //     // halfedge_descriptor h = mesh.halfedge(SM_Face_index(i));
    //     // halfedges_to_remove.push_back(h);
    //     // halfedges_to_remove.push_back(mesh.next(h));
    //     // halfedges_to_remove.push_back(mesh.next(mesh.next(h)));
    //   }
    // }
    // // for(int i = 0; i < halfedges_to_remove.size(); i++) {
    // //   mesh.remove_halfedge(halfedges_to_remove[i]);
    // // }
    // for(int i = 0; i < nfaces2; i++) {
    //   if(remove_on_mesh2[i]) {
    //     mesh2.remove_face(CGAL::SM_Face_index(i));
    //   }
    // }

    // Rcpp::Rcout << "n removed faces: " << mesh.number_of_removed_faces() << "\n";

    // mesh.collect_garbage();
    // mesh2.collect_garbage();

    Rcpp::Rcout << "starting filtering" << "\n";

    EMesh3 umesh1;
    {
      Filtered_graph ffg(umesh, 1, fwhich);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, umesh1, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(umesh, umesh1, f2fmap_, "f:color");
      copy_property<ffg_face_descriptor, face_descriptor, double>(umesh, umesh1, f2fmap_, "f:scalar");

      // // reverse the map between faces
      // // std::map<face_descriptor, face_descriptor> rf2fmap;
      // // for(const auto& [face_ffg, face_fmesh] : f2fmap_) {
      // //   rf2fmap.insert(std::make_pair(
      // //     face_fmesh, CGAL::SM_Vertex_index(int(face_ffg))
      // //   ));
      // // }
      // //
      // if(hasColors) {
      //   Rcpp::Rcout << "adding fcolorMap to umesh1" << "\n";
      //   Fcolors_map fcolorMap_fmesh = 
      //     umesh1.add_property_map<face_descriptor, std::string>(
      //       "f:color", ""
      //     ).first;
      //   Rcpp::Rcout << "starting fcolorMap for umesh1" << "\n";
      //   for(const auto& [face_ffg, face_fmesh] : f2fmap_) {
      //     fcolorMap_fmesh[face_fmesh] = fcolor[CGAL::SM_Face_index(int(face_ffg))];
      //     // rf2fmap.insert(std::make_pair(
      //     //   face_fmesh, CGAL::SM_Face_index(int(face_ffg))
      //     // ));
      //   }
      //   // for(face_descriptor fd : umesh1.faces()) {
      //   //   fcolorMap_fmesh[fd] = fcolor[rf2fmap[fd]];
      //   // }
      // }
      // if(hasScalars) {
      //   Fscalars_map fscalarMap_fmesh = 
      //     umesh1.add_property_map<face_descriptor, double>(
      //       "f:scalar", nan("")
      //     ).first;
      //   for(const auto& [face_ffg, face_fmesh] : f2fmap_) {
      //     fscalarMap_fmesh[face_fmesh] = fscalar[CGAL::SM_Face_index(int(face_ffg))];
      //   }
      // }
    }
    EMesh3 umesh2;
    {
      Filtered_graph ffg(umesh, 2, fwhich);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, umesh2, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(umesh, umesh2, f2fmap_, "f:color");
      copy_property<ffg_face_descriptor, face_descriptor, double>(umesh, umesh2, f2fmap_, "f:scalar");

      // // reverse the map between faces
      // // std::map<face_descriptor, face_descriptor> rf2fmap;
      // // for(const auto& [face_ffg, face_fmesh] : f2fmap_) {
      // //   rf2fmap.insert(std::make_pair(
      // //     face_fmesh, CGAL::SM_Face_index(int(face_ffg))
      // //   ));
      // // }
      // //
      // if(hasColors) {
      //   Fcolors_map fcolorMap_fmesh = 
      //     umesh2.add_property_map<face_descriptor, std::string>(
      //       "f:color", ""
      //     ).first;
      //   Rcpp::Rcout << "starting fcolorMap for umesh2" << "\n";
      //   for(const auto& [face_ffg, face_fmesh] : f2fmap_) {
      //     fcolorMap_fmesh[face_fmesh] = fcolor[CGAL::SM_Face_index(int(face_ffg))];
      //     // rf2fmap.insert(std::make_pair(
      //     //   face_fmesh, CGAL::SM_Face_index(int(face_ffg))
      //     // ));
      //   }
      //   // for(face_descriptor fd : umesh2.faces()) {
      //   //   fcolorMap_fmesh[fd] = fcolor[rf2fmap[fd]];
      //   // }
      // }
      // if(hasScalars) {
      //   Fscalars_map fscalarMap_fmesh = 
      //     umesh2.add_property_map<face_descriptor, double>(
      //       "f:scalar", nan("")
      //     ).first;
      //   for(const auto& [face_ffg, face_fmesh] : f2fmap_) {
      //     fscalarMap_fmesh[face_fmesh] = fscalar[CGAL::SM_Face_index(int(face_ffg))];
      //   }
      // }
    }
    umesh.remove_property_map(fwhich);


    // if(hasColors) {
    //   mesh.remove_property_map(fcolor1);
    //   mesh2.remove_property_map(fcolor2);
    //   Fcolors_map fcolorMap_mesh1 = 
    //     mesh.add_property_map<face_descriptor, std::string>(
    //       "f:color", ""
    //     ).first;
    //   Fcolors_map fcolorMap_mesh2 = 
    //     mesh2.add_property_map<face_descriptor, std::string>(
    //       "f:color", ""
    //     ).first;
    //   for(face_descriptor fd : mesh.faces()) {
    //     fcolorMap_mesh1[fd] = int(fd) < fcolor_mesh1.size() ? fcolor_mesh1[int(fd)] : "";
    //   }
    //   for(face_descriptor fd : mesh2.faces()) {
    //     fcolorMap_mesh2[fd] = int(fd) < fcolor_mesh2.size() ? fcolor_mesh2[int(fd)] : "";
    //   }
    // }

    {
      std::size_t nvisolated = PMP::remove_isolated_vertices(umesh1);
      if(nvisolated > 0) {
        std::string msg;
        if(nvisolated == 1) {
          msg = "Removed one isolated vertex.";
        } else {
          msg = "Removed " + std::to_string(nvisolated) + " isolated vertices.";
        }
        Message(msg);
        umesh1.collect_garbage();
      }
    }
    {
      std::size_t nvisolated = PMP::remove_isolated_vertices(umesh2);
      if(nvisolated > 0) {
        std::string msg;
        if(nvisolated == 1) {
          msg = "Removed one isolated vertex.";
        } else {
          msg = "Removed " + std::to_string(nvisolated) + " isolated vertices.";
        }
        Message(msg);
        umesh2.collect_garbage();
      }
    }




    // Rcpp::IntegerVector Fimap1(mesh.number_of_faces());
    // int i1 = 0;
    // for(EMesh3::Face_index fi : mesh.faces()) {
    //   Fimap1(i1++) = fimap1[fi];
    // }
    // Rcpp::IntegerVector Fimap2(mesh2.number_of_faces());
    // int i2 = 0;
    // for(EMesh3::Face_index fi : mesh2.faces()) {
    //   Fimap2(i2++) = fimap2[fi];
    // }
    // Rcpp::IntegerMatrix Fvis1(fmap_mesh1.size(), 2);
    // {
    //   int ii1 = 0;
    //   for(const auto& [fnew, fsplit] : fmap_mesh1) {
    //     Fvis1(ii1++, Rcpp::_) = Rcpp::IntegerVector({int(fsplit), int(fnew)});
    //   }
    // }
    // Rcpp::IntegerMatrix Fvis2(fmap_mesh2.size(), 2);
    // {
    //   int ii1 = 0;
    //   for(const auto& [fnew, fsplit] : fmap_mesh2) {
    //     Fvis2(ii1++, Rcpp::_) = Rcpp::IntegerVector({int(fsplit), int(fnew)});
    //   }
    // }
    // Rcpp::IntegerMatrix Fvis3(fmap_union.size(), 2);
    // {
    //   int ii1 = 0;
    //   for(const auto& [fnew, fsplit] : fmap_union) {
    //     Fvis3(ii1++, Rcpp::_) = Rcpp::IntegerVector({int(fsplit), int(fnew)});
    //   }
    // }
    // Rcpp::IntegerVector Nfaces(nfaces.size());
    // {
    //   for(int i = 0; i < nfaces.size(); i++) {
    //     Nfaces(i) = nfaces[i];
    //   }
    // }
    // Rcpp::IntegerVector Nfaces2(nfaces2.size());
    // {
    //   for(int i = 0; i < nfaces2.size(); i++) {
    //     Nfaces2(i) = nfaces2[i];
    //   }
    // }

//     Rcpp::List myimesh = Rcpp::List::create(
//       Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(umesh), false),
//       Rcpp::Named("normals") = R_NilValue,
//       Rcpp::Named("vcolors") = R_NilValue,
//       Rcpp::Named("fcolors") = R_NilValue,
//       //Rcpp::Named("fimap") = Fimap,
//       // Rcpp::Named("fimap1") = Fimap1,
//       // Rcpp::Named("fimap2") = Fimap2,
// //      Rcpp::Named("nfaces") = Nfaces,
// //      Rcpp::Named("nfaces2") = Nfaces2,
//       Rcpp::Named("fvis1") = Fvis1,
//       Rcpp::Named("fvis2") = Fvis2,
//       Rcpp::Named("fvis3") = Fvis3
//     );
    return Rcpp::List::create(
      Rcpp::Named("umesh") = Rcpp::XPtr<EMesh3>(new EMesh3(umesh), false),
      Rcpp::Named("mesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(umesh1), false),
      Rcpp::Named("mesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(umesh2), false)
    );
  }


  double volume() {
    if(!CGAL::is_closed(mesh)) {
      Rcpp::stop("The mesh is not closed.");
    }
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The mesh self-intersects.");
    }
    const EK::FT vol = PMP::volume(mesh);
    return CGAL::to_double<EK::FT>(vol);
  }


  void writeFile(
    Rcpp::String filename, const int precision, const bool binary
  ) {
    writeMeshFile(filename, precision, binary, mesh);
  }
  
};