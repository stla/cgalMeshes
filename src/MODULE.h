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
           const bool clean,
           const Rcpp::Nullable<Rcpp::NumericMatrix> &normals_,
           const Rcpp::Nullable<Rcpp::StringVector> &vcolors_,
           const Rcpp::Nullable<Rcpp::StringVector> &fcolors_)
    : mesh(
        csoup2mesh<EMesh3, EPoint3>(
            matrix_to_points3<EPoint3>(vertices), 
            list_to_faces(faces), 
            clean
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
  CGALmesh(Rcpp::List mymesh_)
    : mesh(*(Rcpp::as<Rcpp::XPtr<EMesh3>>(mymesh_["xptr"]).get())),
      mymesh(mymesh_), 
      xptr(mymesh_["xptr"]),
      normals(mymesh_["normals"]),
      vcolors(mymesh_["vcolors"]),
      fcolors(mymesh_["fcolors"]) {}
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
  
  void clipMesh(Rcpp::XPtr<EMesh3> clipperXPtr, const bool clipVolume) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(clipVolume) {
      if(PMP::does_self_intersect(mesh)) {
        Rcpp::stop("The mesh self-intersects.");
      }
    }
    EMesh3 clipper = *(clipperXPtr.get());
    if(!CGAL::is_triangle_mesh(clipper)) {
      Rcpp::stop("The clipping mesh is not triangle.");
    }
    if(!CGAL::is_closed(clipper)) {
      Rcpp::stop("The clipping mesh is not closed.");
    }
    if(PMP::does_self_intersect(clipper)) {
      Rcpp::stop("The clipping mesh self-intersects.");
    }
    
    normals = R_NilValue;
    vcolors = R_NilValue;

    EMesh3 meshcopy;
    if(fcolors.isNotNull()) {
      CGAL::copy_face_graph(mesh, meshcopy);
    }
    
    const bool doNotModify = !clipVolume;

    const bool clipping = PMP::clip(
      mesh, clipper, PMP::parameters::clip_volume(clipVolume),
      PMP::parameters::clip_volume(clipVolume).do_not_modify(doNotModify)
    );
    if(!clipping) {
      Rcpp::stop("Clipping has failed.");
    }

    mesh.collect_garbage();
    
    if(fcolors.isNotNull()) {
      Rcpp::StringVector fcolors0(fcolors);
      const size_t nfaces = meshcopy.number_of_faces();
      std::vector<Triangle> triangles;
      triangles.reserve(nfaces);
      for(EMesh3::Face_index fd : meshcopy.faces()) {
        auto it = vertices_around_face(meshcopy.halfedge(fd), meshcopy);
        auto vd = it.begin();
        triangles.emplace_back(Triangle(meshcopy.point(*(++vd)), meshcopy.point(*(++vd)), meshcopy.point(*(vd))));
      }
      
      const size_t nf = mesh.number_of_faces();
      Rcpp::StringVector fcolors1(nf);
      size_t j = 0;
      for(EMesh3::Face_index f : mesh.faces()) {
        auto it = vertices_around_face(mesh.halfedge(f), mesh);
        auto vd = it.begin();
        Triangle tr(mesh.point(*(++vd)), mesh.point(*(++vd)), mesh.point(*(vd)));
        EPoint3 c = CGAL::centroid(tr);
        size_t k;
        for(k = 0; k < nfaces; k++) {
          Triangle trk = triangles[k];
          if(trk.has_on(c)) {
            break;
          }
        }
        fcolors1(j++) = fcolors0(k);
      }
      Rcpp::Nullable<Rcpp::StringVector> nullable_fcolors(fcolors1);
      fcolors = nullable_fcolors;
    }
    
  }
  
  Rcpp::List clone() {
    EMesh3 copy;
    CGAL::copy_face_graph(mesh, copy);
    // MyMesh mycopy = MYMESH((
    //   xxx.mesh = copy,
    //   xxx.normals = normals,
    //   xxx.vcolors = vcolors,
    //   xxx.fcolors = fcolors
    // ));
    Rcpp::List mycopy = Rcpp::List::create(
      Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(&copy, false),
      Rcpp::Named("normals") = normals,
      Rcpp::Named("vcolors") = vcolors,
      Rcpp::Named("fcolors") = fcolors
    );
    return mycopy;
  }
  
  Rcpp::List connectedComponents(const bool triangulate) {
    // Face_index_map fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
    // const std::size_t ncc = PMP::connected_components(mesh, fccmap);
    Vertex_index_map vccmap = mesh.add_property_map<vertex_descriptor, std::size_t>("v:CC").first;
    Vertex_index_map vindexmap = mesh.add_property_map<vertex_descriptor, std::size_t>("vindex:I").first;
    const std::size_t ncc = connected_components(mesh, vccmap);
    std::vector<std::vector<EMesh3::Vertex_index>> components(ncc);
    std::vector<std::size_t> Indices;
    Indices.reserve(mesh.number_of_vertices());
    std::vector<int> Sizes(ncc, 0);
    Rcpp::NumericMatrix vert = getVertices_EK(mesh);
    int kk = 0;
    for(EMesh3::Vertex_index v : mesh.vertices()) {
      std::size_t c = vccmap[v];
      components[c].push_back(v);
      Indices.emplace_back(c);
      Sizes[c]++;
      int j = int(v);
      if(j < 4) {
        Rcpp::Rcout << "oooooooooo  " << j << "  ooooooooo\n";
        EPoint3 pt = mesh.point(v);
        Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.x()) << "\n";
        Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.y()) << "\n";
        Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.z()) << "\n";
        Rcpp::Rcout << "------------";
        Rcpp::Rcout << vert(0, j) << "\n";
        Rcpp::Rcout << vert(1, j) << "\n";
        Rcpp::Rcout << vert(2, j) << "\n";
        Rcpp::Rcout << "------------";
      }
    }
    // std::vector<int> counters(ncc, 0);
    // for(EMesh3::Vertex_index v : mesh.vertices()) {
    //   std::size_t c = vccmap[v];
    //   vindexmap[v] = components[c][counters[c]++];
    // }
    Face_index_map fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
    std::vector<std::vector<int>> fcomponents(ncc);
    std::vector<std::vector<std::vector<int>>> vcomponents(ncc);
    std::vector<int> fsizes(ncc, 0);
    for(EMesh3::Face_index f : mesh.faces()) {
      auto it = vertices_around_face(mesh.halfedge(f), mesh);
      auto v = it.begin();
      std::size_t c = vccmap[*v];
      fccmap[f] = c;
      fcomponents[c].push_back(int(f));
      fsizes[c]++;
      std::vector<int> theface(3);
      theface[0] = *(v);
      theface[1] = *(++v);
      theface[2] = *(++v);
      vcomponents[c].push_back(theface);
    }
    if(ncc == 1) {
      Message("Only one component found.\n");
    } else {
      const std::string msg = "Found " + std::to_string(ncc) + " components.\n";
      Message(msg);
    }
    const bool really_triangulate = 
      triangulate && !CGAL::is_triangle_mesh(mesh);
    
    Rcpp::List mymeshes(ncc);
    // Face_index_map fmap = mesh.add_property_map<face_descriptor, std::size_t>("f:i").first;
    // Halfedge_index_map hmap = mesh.add_property_map<halfedge_descriptor, std::size_t>("h:i").first;
    // for(EMesh3::Face_index f : mesh.faces()) {
    //   fmap[f] = std::size_t(f);
    // }
    // for(EMesh3::Halfedge_index h : mesh.halfedges()) {
    //   hmap[h] = std::size_t(h);
    // }

    // EMesh3 DD;
    // DD.reserve(mesh.number_of_vertices(), mesh.number_of_edges(), mesh.number_of_faces());
    // for(EMesh3::Vertex_index vi : mesh.vertices()) {
    //   const EPoint3 vertex = mesh.point(vi);
    //   DD.add_vertex(vertex);
    // }
    // Rcpp::Rcout << "\n DD NUmBER VERTICES: " << DD.number_of_vertices() << "\n";
    
    std::vector<EMesh3> ddlist(ncc);
    std::map<EMesh3::Face_index, int> mapfaces;
    int iii = 0;
    for(EMesh3::Face_index f : mesh.faces()) {
      mapfaces[f] = iii++;
    }
    Rcpp::Rcout << "mapfaces done\n";
    //
    std::vector<std::vector<std::vector<EMesh3::Vertex_index>>> connec(ncc);
    for(std::size_t c = 0; c < ncc; c++) {
      std::vector<EMesh3::Vertex_index> component = components[c];
      std::map<EMesh3::Vertex_index,int> newinds;
      for(int j = 0; j < component.size(); j++) {
        newinds[component[j]] = j;
      }
      Rcpp::Rcout << "newinds done\n";
      std::vector<EMesh3::Face_index> todelete;
      for(const auto& [faceindex, intindex] : mapfaces) {
        auto it = vertices_around_face(mesh.halfedge(faceindex), mesh);
        if(vccmap[*(it.begin())] == c) {
          std::vector<EMesh3::Vertex_index> facenewinds;
          for(EMesh3::Vertex_index v : it) {
            //Rcpp::Rcout << "newind: " << newinds[v] << "\n";
            facenewinds.push_back(CGAL::SM_Vertex_index(newinds[v]));
          }
          //Rcpp::Rcout << "facenewinds done\n";
          connec[c].push_back(facenewinds);
          todelete.push_back(faceindex);
        }
      }
      for(int fi = 0; fi < todelete.size(); fi++) {
        mapfaces.erase(todelete[fi]);
        //Rcpp::Rcout << "facedindex erased";
      }
      EMesh3 dd;
      for(int j = 0; j < component.size(); j++) {
        EPoint3 pt = mesh.point(component[j]);
        dd.add_vertex(pt);
      }
      Rcpp::Rcout << "vertices added to dd\n";
      ddlist[c] = dd;
    }
    std::vector<EMesh3> ddlist2(ncc);
    for(std::size_t c = 0; c < ncc; c++) {
      EMesh3 dd = ddlist[c];
      for(int j = 0; j < connec[c].size(); j++) {
        dd.add_face(connec[c][j]);
      }
      Rcpp::Rcout << "faces added to dd\n";
      ddlist2[c] = dd;
      //Rcpp::Rcout << "component " + std::to_string(c) << "\n";
      // Filtered_graph ffg(mesh, c, fccmap);//, CGAL::parameters::face_index_map(fmap).halfedge_index_map(hmap).vertex_index_map(vindexmap));
      // assert(ffg.is_selection_valid());
      // EMesh3 cc;
      // CGAL::copy_face_graph(ffg, cc);
      // Rcpp::Rcout << "\nGRAPH COPIED\n";
      // const size_t nvertices = cc.number_of_vertices();
      // const size_t nedges    = cc.number_of_edges();
      // const size_t nfaces    = cc.number_of_faces();
      // EMesh3 dd;
      // dd.reserve(mesh.number_of_vertices(), mesh.number_of_edges(), mesh.number_of_faces());
      // for(EMesh3::Vertex_index vi : mesh.vertices()) {
      //   const EPoint3 vertex = mesh.point(vi);
      //   dd.add_vertex(vertex);
      // }
      // CGAL::copy_face_graph(DD, dd);
      // Rcpp::Rcout << "\n dd NUmBER VERTICES: " << dd.number_of_vertices() << "\n";
      // for(int ii = 0; ii < fsizes[c]; ii++) {
      //   std::vector<int> ddface = vcomponents[c][ii];
      //   Rcpp::Rcout << "FFFFFFFFFFFFFFFFF  " << ddface.size() << "  FFFFFFFFFFFFFFFFFFFFFFFF\n";
      //   for(int uuu = 0; uuu < 3; uuu++) {
      //     Rcpp::Rcout << "GGGGGGGGGGGG  " << ddface[uuu] << "  GGGGGGGGGGGGGG\n";
      //   }
      //   dd.add_face(
      //     CGAL::SM_Vertex_index(ddface[0]), 
      //     CGAL::SM_Vertex_index(ddface[1]), 
      //     CGAL::SM_Vertex_index(ddface[2])
      //   );
      // }
      // std::size_t nrem = PMP::remove_isolated_vertices(dd);
      // dd.reserve(nvertices, nedges, nfaces);
      // Rcpp::Rcout << "\ntttttttttttttt  " << nfaces << "  tttttttttttt\n";
      // for(int ii = 0; ii < nfaces; ii++) {
      //   Rcpp::Rcout << "FFFFFFFFFFFFFFFFF  " << ii << "  FFFFFFFFFFFFFFFFFFFFFFFF\n";
      //   std::vector<int> ddface = vcomponents[c][ii];
      //   Rcpp::Rcout << "FFFFFFFFFFFFFFFFF  " << ddface.size() << "  FFFFFFFFFFFFFFFFFFFFFFFF\n";
      //   for(int uuu = 0; uuu < 3; uuu++) {
      //     Rcpp::Rcout << "GGGGGGGGGGGG  " << ddface[uuu] << "  GGGGGGGGGGGGGG\n";
      //   }
      //   const EPoint3 vertex1 = mesh.point(CGAL::SM_Vertex_index(ddface[0]));
      //   const EPoint3 vertex2 = mesh.point(CGAL::SM_Vertex_index(ddface[1]));
      //   const EPoint3 vertex3 = mesh.point(CGAL::SM_Vertex_index(ddface[2]));
      //   dd.add_vertex(EPoint3(vertex1.x(), vertex1.y(), vertex1.z()));
      //   dd.add_vertex(EPoint3(vertex2.x(), vertex2.y(), vertex2.z()));
      //   dd.add_vertex(EPoint3(vertex3.x(), vertex3.y(), vertex3.z()));
      // }
      // ddlist[c] = dd;
      //////
      Rcpp::Rcout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
      if(c == 0) {
        for(EMesh3::Vertex_index v : dd.vertices()) {
          int j = int(v);
          if(j < 12) {
            Rcpp::Rcout << "XXXXXXXXXX  " << j << "  XXXXXXXXXX\n";
            EPoint3 pt = dd.point(v);
            Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.x()) << "\n";
            Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.y()) << "\n";
            Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.z()) << "\n";
            Rcpp::Rcout << "------------";
            Rcpp::Rcout << vert(0, int(components[c][j])) << "\n";
            Rcpp::Rcout << vert(1, int(components[c][j])) << "\n";
            Rcpp::Rcout << vert(2, int(components[c][j])) << "\n";
            Rcpp::Rcout << "------------";
          }
        }
        //////////
        // for(EMesh3::Vertex_index v : dd.vertices()) {
        //   int j = int(v);
        //   if(j < 12) {
        //     Rcpp::Rcout << "YYYYYYYYYYY  " << j << "  YYYYYYYYYYY\n";
        //     EPoint3 pt = dd.point(v);
        //     Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.x()) << "\n";
        //     Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.y()) << "\n";
        //     Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.z()) << "\n";
        //     Rcpp::Rcout << "------------";
        //     Rcpp::Rcout << vert(0, components[c][j]) << "\n";
        //     Rcpp::Rcout << vert(1, components[c][j]) << "\n";
        //     Rcpp::Rcout << vert(2, components[c][j]) << "\n";
        //     Rcpp::Rcout << "------------";
        //   }
        // }
      }
      // EMesh3 cc;
      // CGAL::copy_face_graph(mesh, cc);
      // std::vector<std::size_t> keep = {i};
      // Vertex_index_map vmap = mesh.add_property_map<vertex_descriptor, std::size_t>("v:index").first;
      // std::size_t k = 0;
      // for(EMesh3::Vertex_index v : mesh.vertices()) {
      //   vmap[v] = k++;
      // }
      // PMP::keep_connected_components(cc, keep, fccmap, PMP::parameters::vertex_index_map(vmap));
      // Rcpp::Rcout << "component " + std::to_string(i) << "\n";
      // for(EMesh3::Vertex_index v : cc.vertices()) {
      //   Rcpp::Rcout << vmap[v] << "\n";
      // }
      if(really_triangulate) {
        const bool success = PMP::triangulate_faces(dd);
        if(!success) {
          const std::string msg = "Triangulation has failed (component " +
            std::to_string(c + 1) + ").";
          Rcpp::stop(msg);
        }
      }
      
      //      MyMesh* mycc = new MyMesh; 
      // Rcpp::NumericMatrix cnormals(3, Sizes[c]);
      // Rcpp::NumericMatrix dnormals(3, Sizes[c]);
      // dd.collect_garbage();
      // Rcpp::Rcout << "Sizes[c]  " << Sizes[c] << "  Sizes[c]\n";
      // Rcpp::Rcout << "Sizes[c]  " << dd.number_of_vertices() << "  Sizes[c]\n";
      // Rcpp::Rcout << "Sizes[c]  " << dd.number_of_faces() << "  Sizes[c]\n";
      // Rcpp::Rcout << "Sizes[c]  " << fsizes[c] << "  Sizes[c]\n";
      // Rcpp::Rcout << "Sizes[c]  " << fcomponents[c].size() << "  Sizes[c]\n";
      // if(normals.isNotNull()) {
      //   Rcpp::NumericMatrix normals0(normals);
      //   
      //   for(EMesh3::Face_index fd : dd.faces()) {
      //     std::vector<int> ccface;
      //     for(EMesh3::Vertex_index vd : vertices_around_face(dd.halfedge(fd), dd)) {
      //       ccface.push_back(int(vd));
      //       Rcpp::Rcout << "VD  " << vd << "  VD\n";
      //     }
      //     int ff = fcomponents[c][int(fd)];
      //     Rcpp::Rcout << "DDDDDDDDDDD  " << ff << "  DDDDDDDDDDDD\n";
      //     std::vector<int> theface = vcomponents[c][ff];
      //     Rcpp::Rcout << "EEEE  " << theface[0] << " " << theface[1] << " " << theface[2] << "  DDDDDDDDDDDD\n";
      //     dnormals(Rcpp::_, ccface[0]) = normals0(Rcpp::_, theface[0]);
      //     dnormals(Rcpp::_, ccface[1]) = normals0(Rcpp::_, theface[1]);
      //     dnormals(Rcpp::_, ccface[2]) = normals0(Rcpp::_, theface[2]);
      //   }
        
        // int i = 0;
        // for(EMesh3::Vertex_index v : mesh.vertices()) {
        //   if(vccmap[v] == c){
        //     Rcpp::NumericVector col = normals0(Rcpp::_, int(v));
        //     cnormals(Rcpp::_, i++) = col;
        //   };
        // }
        // for(EMesh3::Face_index fd : cc.faces()) {
        //   std::vector<int> ccface;
        //   for(EMesh3::Vertex_index vd : vertices_around_face(cc.halfedge(fd), cc)) {
        //     ccface.push_back(int(vd));
        //   }
        //   int ff = fcomponents[c][int(fd)];
        //   Rcpp::Rcout << "DDDDDDDDDDD  " << ff << "  DDDDDDDDDDDD\n";
        //   std::vector<int> theface = vcomponents[c][ff];
        //   Rcpp::Rcout << "EEEE  " << theface[0] << " " << theface[1] << " " << theface[2] << "  DDDDDDDDDDDD\n";
        //   cnormals(Rcpp::_, ccface[0]) = normals0(Rcpp::_, theface[0]);
        //   cnormals(Rcpp::_, ccface[1]) = normals0(Rcpp::_, theface[1]);
        //   cnormals(Rcpp::_, ccface[2]) = normals0(Rcpp::_, theface[2]);
        // }
        // for(int i = 0; i < Sizes[c]; i++) {
        //   Rcpp::NumericVector col = normals0(Rcpp::_, components[c][i]);
        //   cnormals(Rcpp::_, i) = col;
        // }
        
        // for(int i = 0; i < Sizes[c]; i++) {
        //   Rcpp::NumericVector col = normals0(Rcpp::_, components[c][i]);
        //   cnormals(Rcpp::_, i) = col;
        // }
//        Rcpp::Rcout << Sizes[c];
//        Rcpp::Rcout << cnormals.ncol();
//        mycc->normals = Rcpp::Nullable<Rcpp::NumericMatrix>(cnormals); // !!!! achtung si tu triangules !!!
//      }
//      Rcpp::Nullable<Rcpp::NumericMatrix> test = mycc->normals;
      // Rcpp::Rcout << "\n" << test.isNotNull();
      // Rcpp::RObject robj = Rcpp::wrap(mycc->normals);
      // int obj_type = robj.sexp_type();
      // Rcpp::Rcout << "\n---" << obj_type << "\n---";
      //delete mycc;
    }
    std::vector<Rcpp::Nullable<Rcpp::NumericMatrix>> NORMALS(ncc);
    std::vector<Rcpp::Nullable<Rcpp::StringVector>> VCOLORS(ncc);
    std::vector<Rcpp::Nullable<Rcpp::StringVector>> FCOLORS(ncc);
    if(normals.isNotNull()) {
      Rcpp::NumericMatrix normals0(normals);
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Vertex_index> component = components[c];
        size_t csize = component.size();
        Rcpp::NumericMatrix cnormals(3, csize);
        for(size_t j = 0; j < csize; j++) {
          cnormals(Rcpp::_, j) = normals0(Rcpp::_, int(component[j]));
        }
        NORMALS[c] = Rcpp::Nullable<Rcpp::NumericMatrix>(cnormals);
      }
    } else {
      for(std::size_t c = 0; c < ncc; c++) {
        NORMALS[c] = R_NilValue;
      }
    }
    if(vcolors.isNotNull()) {
      Rcpp::StringVector vcolors0(vcolors);
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Vertex_index> component = components[c];
        size_t csize = component.size();
        Rcpp::StringVector cvcolors(csize);
        for(size_t j = 0; j < csize; j++) {
          cvcolors(j) = vcolors0(int(component[j]));
        }
        VCOLORS[c] = Rcpp::Nullable<Rcpp::StringVector>(cvcolors);
      }
    } else {
      for(std::size_t c = 0; c < ncc; c++) {
        VCOLORS[c] = R_NilValue;
      }
    }
    if(fcolors.isNotNull()) {
      Rcpp::StringVector fcolors0(fcolors);
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Vertex_index> component = components[c];
        size_t csize = component.size();
        Rcpp::StringVector cfcolors(csize);
        for(size_t j = 0; j < csize; j++) {
          cfcolors(j) = fcolors0(int(component[j]));
        }
        FCOLORS[c] = Rcpp::Nullable<Rcpp::StringVector>(cfcolors);
      }
    } else {
      for(std::size_t c = 0; c < ncc; c++) {
        FCOLORS[c] = R_NilValue;
      }
    }
    for(std::size_t c = 0; c < ncc; c++) {
      Rcpp::List MM = Rcpp::List::create(
        Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(ddlist2[c]), false),
        Rcpp::Named("normals") = NORMALS[c],
                                        Rcpp::Named("vcolors") = VCOLORS[c],
                                                                        Rcpp::Named("fcolors") = FCOLORS[c]
      );
      mymeshes(c) = MM;
    }
    
    // for(int yyy = 0; yyy < ncc; yyy++) {
    //   EMesh3 dd = ddlist[yyy];
    //   for(int ii = 0; ii < fsizes[yyy]; ii++) {
    //     std::vector<int> ddface = vcomponents[yyy][ii];
    //     dd.add_face(
    //       CGAL::SM_Vertex_index(ddface[0]), 
    //       CGAL::SM_Vertex_index(ddface[1]), 
    //       CGAL::SM_Vertex_index(ddface[2])
    //     );
    //   }
    // }
    
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
      // MyMesh mycmesh = MYMESH((
      //   xxx.mesh = cmesh
      // ));
      Rcpp::List mycmesh = Rcpp::List::create(
        Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(cmesh), false),
        Rcpp::Named("normals") = R_NilValue,
        Rcpp::Named("vcolors") = R_NilValue,
        Rcpp::Named("fcolors") = R_NilValue
      );
      out(i) = mycmesh;
      i++;
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
      distances(i) = PMP::max_distance_to_triangle_mesh<CGAL::Sequential_tag>(
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
    const bool success = PMP::fair(mesh, selectedVertices);
    if(!success) {
      Rcpp::stop("Failed to fair the mesh.");
    }
    normals = R_NilValue;
    vcolors = R_NilValue;
    fcolors = R_NilValue;
  }
  
  Rcpp::NumericVector geoDists(const int index) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    int nvertices = mesh.number_of_vertices();
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
      gdistances(i) = get(vertex_distance, vd);
      i++;
    }
    return gdistances;    
  }
  
  Rcpp::Nullable<Rcpp::NumericMatrix> getNormals() {
    if(normals.isNotNull()) {
      Rcpp::NumericMatrix normals0(normals);
      Rcpp::Rcout << normals0.ncol() << " " << normals0.nrow();
      return normals;
    } else {
      Rcpp::Rcout << "null";
      return R_NilValue;
    }
  }
  
  Rcpp::List getRmesh() {
    const bool there_is_normals = normals.isNotNull();
    Rcpp::List rmesh;
    if(CGAL::is_triangle_mesh(mesh)) {
      rmesh = RSurfEKMesh2(mesh, !there_is_normals, 3);
    } else if(CGAL::is_quad_mesh(mesh)) {
      rmesh = RSurfEKMesh2(mesh, !there_is_normals, 4);
    } else {
      rmesh = RSurfEKMesh(mesh, !there_is_normals);
    }
    if(there_is_normals) {
      rmesh["normals"] = normals;
    }
    if(vcolors.isNotNull() || fcolors.isNotNull()) {
      const size_t nf = mesh.number_of_faces();
      if(vcolors.isNotNull()) {
        Rcpp::StringVector colors(vcolors);
        if(colors.size() == nf) {
          rmesh["colors"] = vcolors;
        }
      } else {
        rmesh["colors"] = fcolors;
      }
    }
    return rmesh;
  }
  
  Rcpp::List intersection(Rcpp::XPtr<EMesh3> mesh2XPtr) {
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
    // MyMesh myimesh = MYMESH((
    //   xxx.mesh = imesh
    // ));
    Rcpp::List myimesh = Rcpp::List::create(
      Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(imesh), false),
      Rcpp::Named("normals") = R_NilValue,
      Rcpp::Named("vcolors") = R_NilValue,
      Rcpp::Named("fcolors") = R_NilValue
    );
    return myimesh;
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
  
  bool isTriangle() {
    return CGAL::is_triangle_mesh(mesh);
  }

  bool isValid() {
    return CGAL::is_valid_polygon_mesh(mesh);
  }
  
  void orientToBoundVolume() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    PMP::orient_to_bound_a_volume(mesh);
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
    normals = R_NilValue;
    vcolors = R_NilValue;
    fcolors = R_NilValue;
    mesh.collect_garbage();
  }
  
  void reverseFaceOrientations() {
    PMP::reverse_face_orientations(mesh);  
  }
  
  Rcpp::List subtract(Rcpp::XPtr<EMesh3> mesh2XPtr) {
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
    // MyMesh myimesh = MYMESH((
    //   xxx.mesh = imesh
    // ));
    Rcpp::List myimesh = Rcpp::List::create(
      Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(imesh), false),
      Rcpp::Named("normals") = R_NilValue,
      Rcpp::Named("vcolors") = R_NilValue,
      Rcpp::Named("fcolors") = R_NilValue
    );
    return myimesh;
  }
  
  void triangulate() {
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
    normals = R_NilValue;
    fcolors = R_NilValue;
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
    EMesh3 imesh;
    const bool success = PMP::corefine_and_compute_union(
      mesh, mesh2, imesh
    );
    if(!success) {
      Rcpp::stop("Union computation has failed.");
    }
    // MyMesh myimesh = MYMESH((
    //   xxx.mesh = imesh
    // ));
    Rcpp::List myimesh = Rcpp::List::create(
      Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(imesh), false),
      Rcpp::Named("normals") = R_NilValue,
      Rcpp::Named("vcolors") = R_NilValue,
      Rcpp::Named("fcolors") = R_NilValue
    );
    return myimesh;
  }
  
  Rcpp::NumericMatrix vertices() {
    return getVertices_EK(mesh);
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