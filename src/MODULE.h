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
//     EMesh3 splitter = *(clipperXPtr.get());
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
  
  void clipMesh(Rcpp::XPtr<EMesh3> clipperXPtr, const bool clipVolume) {
    EMesh3 clipper = *(clipperXPtr.get());
    clipping(mesh, clipper, clipVolume);
  }

  Rcpp::XPtr<EMesh3> doubleclip(Rcpp::XPtr<EMesh3> clipperXPtr) {
    EMesh3 tm1;
    CGAL::copy_face_graph(mesh, tm1);
    EMesh3 clipper = *(clipperXPtr.get());
    EMesh3 clipper1;
    CGAL::copy_face_graph(clipper, clipper1);

    Fcolors_map fcolorsmap1 = 
      mesh.property_map<face_descriptor, std::string>("f:color").first;
    Fcolors_map fcolorsmap11 = 
      tm1.add_property_map<face_descriptor, std::string>("f:color").first;
    for(EMesh3::Face_index fi : mesh.faces()) {
      fcolorsmap11[fi] = fcolorsmap1[fi];
    }
    Fcolors_map fcolorsmap2 = 
      clipper.property_map<face_descriptor, std::string>("f:color").first;
    Fcolors_map fcolorsmap22 = 
      clipper1.add_property_map<face_descriptor, std::string>("f:color").first;
    for(EMesh3::Face_index fi : clipper.faces()) {
      fcolorsmap22[fi] = fcolorsmap2[fi];
    }
     
    clipping(mesh, clipper, false);
    clipping(clipper1, tm1, false);
    const int nfaces = mesh.number_of_faces();
    std::map<face_descriptor, std::string> fcolorstdmap;
    Fcolors_map fcolorsmap111 = 
      mesh.property_map<face_descriptor, std::string>("f:color").first;
    for(EMesh3::Face_index fi : mesh.faces()) {
      fcolorstdmap[fi] = fcolorsmap111[fi];
    }
    Fcolors_map fcolorsmap222 = 
      clipper1.property_map<face_descriptor, std::string>("f:color").first;
    for(EMesh3::Face_index fi : clipper1.faces()) {
      fcolorstdmap[CGAL::SM_Face_index(nfaces + int(fi))] = fcolorsmap222[fi];
    }

    const int nvertices = mesh.number_of_vertices();



    EMesh3 out;
    for(EMesh3::Vertex_index vi : mesh.vertices()) {
      out.add_vertex(mesh.point(vi));
    }
    for(EMesh3::Vertex_index vi : clipper1.vertices()) {
      out.add_vertex(clipper1.point(vi));
    }
    for(EMesh3::Face_index fi : mesh.faces()) {
      auto it = vertices_around_face(mesh.halfedge(fi), mesh);
      std::vector<EMesh3::Vertex_index> face(it.begin(), it.end());
      out.add_face(face);
    }
    for(EMesh3::Face_index fi : clipper1.faces()) {
      auto it = vertices_around_face(clipper1.halfedge(fi), clipper1);
      std::vector<EMesh3::Vertex_index> face;
      for(EMesh3::Vertex_index v: it) {
        face.push_back(CGAL::SM_Vertex_index(nvertices + int(v)));
      }
      out.add_face(face);
    }
    Fcolors_map fcolorsmap = 
      out.add_property_map<face_descriptor, std::string>("f:color").first;
    for(EMesh3::Face_index fi : out.faces()) {
      fcolorsmap[fi] = fcolorstdmap[fi];
    }

    return Rcpp::XPtr<EMesh3>(new EMesh3(out), false);
  }

  
  Rcpp::List clone() {
    EMesh3 copy;
    CGAL::copy_face_graph(mesh, copy);
    Rcpp::List mycopy = Rcpp::List::create(
      Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(copy), false),
      Rcpp::Named("normals") = normals,
      Rcpp::Named("vcolors") = vcolors,
      Rcpp::Named("fcolors") = fcolors
    );
    return mycopy;
  }
  
  Rcpp::List connectedComponents(const bool triangulate) {
    
/*     Vertex_index_map vmap = mesh.add_property_map<vertex_descriptor, std::size_t>("v:i").first;
    Face_index_map fmap = mesh.add_property_map<face_descriptor, std::size_t>("f:i").first;
    Rcpp::Rcout << "PROPERTY MAPS CREATED" << "\n";
    
    std::vector<EMesh3> cc_meshes;
    PMP::split_connected_components(mesh, cc_meshes);//, PMP::parameters::vertex_index_map(vmap).face_index_map(fmap));

    Rcpp::Rcout << "FPROPERTIES" << "\n";
    std::vector<std::string> props = mesh.properties<face_descriptor>();
    for(int p = 0; p < props.size(); p++) {
      Rcpp::Rcout << props[p] << "\n";
    }
    Rcpp::Rcout << "VPROPERTIES" << "\n";
    std::vector<std::string> vprops = mesh.properties<vertex_descriptor>();
    for(int p = 0; p < vprops.size(); p++) {
      Rcpp::Rcout << vprops[p] << "\n";
    }

    Rcpp::IntegerVector findex(mesh.number_of_faces());
    std::pair<Face_index_map, bool> xmaybemap = mesh.property_map<face_descriptor, std::size_t>("f:i");
    Rcpp::Rcout << "fi: " << xmaybemap.second << "\n";
    if(xmaybemap.second) {
      int ffff = 0;
      for(EMesh3::Face_index fi : mesh.faces()) {
        findex(ffff++) = xmaybemap.first[fi];
        Rcpp::Rcout << xmaybemap.first[fi] << "\n";
      }
    }
    Rcpp::IntegerVector vindex(mesh.number_of_vertices());
    std::pair<Vertex_index_map, bool> xxmaybemap = mesh.property_map<vertex_descriptor, std::size_t>("v:i");
    Rcpp::Rcout << "vi: " << xxmaybemap.second << "\n";
    if(xxmaybemap.second) {
      int ffff = 0;
      for(EMesh3::Vertex_index vi : mesh.vertices()) {
        vindex(ffff++) = xxmaybemap.first[vi];
        Rcpp::Rcout << xxmaybemap.first[vi] << "\n";
      }
    }
    
    Rcpp::IntegerVector fconnectivity(mesh.number_of_faces());
    std::pair<Face_index_map, bool> xxxmaybemap = mesh.property_map<face_descriptor, std::size_t>("f:connectivity");
    Rcpp::Rcout << "fconnectivity: " << xxxmaybemap.second << "\n";
    if(xxxmaybemap.second) {
      int ffff = 0;
      for(EMesh3::Face_index fi : mesh.faces()) {
        fconnectivity(ffff++) = xxxmaybemap.first[fi];
        Rcpp::Rcout << xxxmaybemap.first[fi] << "\n";
      }
    }

    Rcpp::LogicalVector fremoved(mesh.number_of_faces());
    std::pair<EMesh3::Property_map<face_descriptor, bool>, bool> xxxxmaybemap = mesh.property_map<face_descriptor, bool>("f:removed");
    Rcpp::Rcout << "fremoved: " << xxxxmaybemap.second << "\n";
    if(xxxxmaybemap.second) {
      int fffff = 0;
      for(EMesh3::Face_index fi : mesh.faces()) {
        fremoved(fffff++) = xxxxmaybemap.first[fi];
        Rcpp::Rcout << xxxxmaybemap.first[fi] << "\n";
      }
    }
    

    Rcpp::Rcout << "FPROPERTIES CONCOM" << "\n";
    std::vector<std::string> cprops = cc_meshes[0].properties<face_descriptor>();
    for(int p = 0; p < cprops.size(); p++) {
      Rcpp::Rcout << cprops[p] << "\n";
    }
    Rcpp::Rcout << "VPROPERTIES CONCOM" << "\n";
    std::vector<std::string> vcprops = cc_meshes[0].properties<vertex_descriptor>();
    for(int p = 0; p < vcprops.size(); p++) {
      Rcpp::Rcout << vcprops[p] << "\n";
    }
 */    
    // std::pair<EMesh3::Property_map<face_descriptor, std::size_t>, bool> maybemap = mesh.property_map<face_descriptor, std::size_t>("f:i");
    // Rcpp::IntegerVector findices2(mesh.number_of_faces());
    // if(maybemap.second) {
    //   int fff = 0;
    //   for(EMesh3::Face_index fi : mesh.faces()) {
    //     findices2(fff++) = maybemap.first[fi];
    //   }
    // }
    
    //return Rcpp::List::create(Rcpp::Named("fuck") = 0);
    ///////////////////////////////////////////////////
    
/*     const size_t ncc = cc_meshes.size();
    if(ncc == 1) {
      Message("Only one component found.\n");
    } else {
      const std::string msg = "Found " + std::to_string(ncc) + " components.\n";
      Message(msg);
    }
    const bool really_triangulate = 
      triangulate && !CGAL::is_triangle_mesh(mesh);
    Rcpp::List xptrs(ncc);
    int i = 0;
    for(auto cc = cc_meshes.begin(); cc != cc_meshes.end(); ++cc) {
      if(really_triangulate) {
        const bool success = PMP::triangulate_faces(*cc);
        if(!success) {
          const std::string msg = "Triangulation has failed (component " +
            std::to_string(i + 1) + ").";
          Rcpp::stop(msg);
        }
      }
      xptrs(i) = Rcpp::XPtr<EMesh3>(new EMesh3(*cc), false);
      i++;
    }
    return xptrs;
 */    
    
    
//     // Face_index_map fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
//     // const std::size_t ncc = PMP::connected_components(mesh, fccmap);

    std::pair<Fcolors_map, bool> test_ = mesh.property_map<face_descriptor, std::string>("f:color");
    Rcpp::Rcout << "has fcolor: " << test_.second << "\n";

    Vertex_index_map vccmap = mesh.add_property_map<vertex_descriptor, std::size_t>("v:CC").first;
    // Vertex_index_map vindexmap = mesh.add_property_map<vertex_descriptor, std::size_t>("vindex:I").first;
    const std::size_t ncc = connected_components(mesh, vccmap);
    std::vector<std::vector<EMesh3::Vertex_index>> components(ncc);
    // std::vector<std::size_t> Indices;
    // Indices.reserve(mesh.number_of_vertices());
    // std::vector<int> Sizes(ncc, 0);
    //Rcpp::NumericMatrix vert = getVertices_EK(mesh);
    // int kk = 0;
    for(EMesh3::Vertex_index v : mesh.vertices()) {
      std::size_t c = vccmap[v];
      components[c].push_back(v);
      // Indices.emplace_back(c);
      // Sizes[c]++;
      // int j = int(v);
      // if(j < 4) {
      //   Rcpp::Rcout << "oooooooooo  " << j << "  ooooooooo\n";
      //   EPoint3 pt = mesh.point(v);
      //   Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.x()) << "\n";
      //   Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.y()) << "\n";
      //   Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.z()) << "\n";
      //   Rcpp::Rcout << "------------";
      //   Rcpp::Rcout << vert(0, j) << "\n";
      //   Rcpp::Rcout << vert(1, j) << "\n";
      //   Rcpp::Rcout << vert(2, j) << "\n";
      //   Rcpp::Rcout << "------------";
      // }
    }
    // std::vector<int> counters(ncc, 0);
    // for(EMesh3::Vertex_index v : mesh.vertices()) {
    //   std::size_t c = vccmap[v];
    //   vindexmap[v] = components[c][counters[c]++];
    // }
    // Face_index_map fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
    // std::vector<std::vector<int>> fcomponents(ncc);
    // std::vector<std::vector<std::vector<int>>> vcomponents(ncc);
    // std::vector<int> fsizes(ncc, 0);
    // for(EMesh3::Face_index f : mesh.faces()) {
    //   auto it = vertices_around_face(mesh.halfedge(f), mesh);
    //   auto v = it.begin();
    //   std::size_t c = vccmap[*v];
    //   fccmap[f] = c;
    //   fcomponents[c].push_back(int(f));
    //   fsizes[c]++;
    //   std::vector<int> theface(3);
    //   theface[0] = *(v);
    //   theface[1] = *(++v);
    //   theface[2] = *(++v);
    //   vcomponents[c].push_back(theface);
    // }
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
    
    std::vector<EMesh3> ccmeshes0(ncc);
    std::map<EMesh3::Face_index, int> mapfaces;
    int iii = 0;
    for(EMesh3::Face_index f : mesh.faces()) {
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
        auto it = vertices_around_face(mesh.halfedge(faceindex), mesh);
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
        EPoint3 pt = mesh.point(component[j]);
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
      // if(c == 0) {
      //   for(EMesh3::Vertex_index v : dd.vertices()) {
      //     int j = int(v);
      //     if(j < 12) {
      //       Rcpp::Rcout << "XXXXXXXXXX  " << j << "  XXXXXXXXXX\n";
      //       EPoint3 pt = dd.point(v);
      //       Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.x()) << "\n";
      //       Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.y()) << "\n";
      //       Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.z()) << "\n";
      //       Rcpp::Rcout << "------------";
      //       Rcpp::Rcout << vert(0, int(components[c][j])) << "\n";
      //       Rcpp::Rcout << vert(1, int(components[c][j])) << "\n";
      //       Rcpp::Rcout << vert(2, int(components[c][j])) << "\n";
      //       Rcpp::Rcout << "------------";
      //     }
      //   }
      //   //////////
      //   // for(EMesh3::Vertex_index v : dd.vertices()) {
      //   //   int j = int(v);
      //   //   if(j < 12) {
      //   //     Rcpp::Rcout << "YYYYYYYYYYY  " << j << "  YYYYYYYYYYY\n";
      //   //     EPoint3 pt = dd.point(v);
      //   //     Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.x()) << "\n";
      //   //     Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.y()) << "\n";
      //   //     Rcpp::Rcout << CGAL::to_double<EK::FT>(pt.z()) << "\n";
      //   //     Rcpp::Rcout << "------------";
      //   //     Rcpp::Rcout << vert(0, components[c][j]) << "\n";
      //   //     Rcpp::Rcout << vert(1, components[c][j]) << "\n";
      //   //     Rcpp::Rcout << vert(2, components[c][j]) << "\n";
      //   //     Rcpp::Rcout << "------------";
      //   //   }
      //   // }
      // }
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



    std::pair<Normals_map, bool> vnormals_ = mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
    if(vnormals_.second) {
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Vertex_index> component = components[c];
        Normals_map vnormals = ccmeshes[c].add_property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal").first;
        for(EMesh3::Vertex_index vi : ccmeshes[c].vertices()) {
          vnormals[vi] = vnormals_.first[CGAL::SM_Vertex_index(component[int(vi)])];
        }
      }
    }
    std::pair<Vcolors_map, bool> vcolors_ = mesh.property_map<vertex_descriptor, std::string>("v:color");
    if(vcolors_.second) {
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Vertex_index> component = components[c];
        Vcolors_map vcolors = ccmeshes[c].add_property_map<vertex_descriptor, std::string>("v:color").first;
        for(EMesh3::Vertex_index vi : ccmeshes[c].vertices()) {
          vcolors[vi] = vcolors_.first[CGAL::SM_Vertex_index(component[int(vi)])];
        }
      }
    }
    std::pair<Fcolors_map, bool> fcolors_ = mesh.property_map<face_descriptor, std::string>("f:color");
    Rcpp::Rcout << "has fcolor: " << fcolors_.second << "\n";
    if(fcolors_.second) {
      for(std::size_t c = 0; c < ncc; c++) {
        std::vector<EMesh3::Face_index> fcomponent = todelete[c];
        Fcolors_map fcolors = ccmeshes[c].add_property_map<face_descriptor, std::string>("f:color").first;
        for(EMesh3::Face_index fi : ccmeshes[c].faces()) {
          fcolors[fi] = fcolors_.first[CGAL::SM_Face_index(fcomponent[int(fi)])];
        }
      }
    }

    for(std::size_t c = 0; c < ncc; c++) {
      mymeshes(c) = Rcpp::XPtr<EMesh3>(new EMesh3(ccmeshes[c]), false);
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
        size_t csize = todelete[c].size();
        Rcpp::StringVector cfcolors(csize);
        for(size_t j = 0; j < csize; j++) {
          cfcolors(j) = fcolors0(int(todelete[c][j]));
        }
        FCOLORS[c] = Rcpp::Nullable<Rcpp::StringVector>(cfcolors);
      }
    } else {
      for(std::size_t c = 0; c < ncc; c++) {
        FCOLORS[c] = R_NilValue;
      }
    }
    // for(std::size_t c = 0; c < ncc; c++) {
    //   Rcpp::List MM = Rcpp::List::create(
    //     Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(ddlist2[c]), false),
    //     Rcpp::Named("normals") = NORMALS[c],
    //     Rcpp::Named("vcolors") = VCOLORS[c],
    //     Rcpp::Named("fcolors") = FCOLORS[c]
    //   );
    //   mymeshes(c) = MM;
    // }
    
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
    return gdistances;    
  }
  
  Rcpp::Nullable<Rcpp::NumericMatrix> getNormals() {
    return normals;
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

  Rcpp::List getRmesh() {
    std::pair<Normals_map, bool> normalsmap_ = 
      mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
    const bool there_is_normals = normalsmap_.second;
    Rcpp::List rmesh;
    if(CGAL::is_triangle_mesh(mesh)) {
      rmesh = RSurfEKMesh2(mesh, !there_is_normals, 3);
    } else if(CGAL::is_quad_mesh(mesh)) {
      rmesh = RSurfEKMesh2(mesh, !there_is_normals, 4);
    } else {
      rmesh = RSurfEKMesh(mesh, !there_is_normals);
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
    normals = R_NilValue;
    vcolors = R_NilValue;
    fcolors = R_NilValue;
    mesh.collect_garbage();
  }
  
  void reverseFaceOrientations() {
    PMP::reverse_face_orientations(mesh);
    // update normals
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

    EMesh3 mcopy;
    CGAL::copy_face_graph(mesh, mcopy);


    UnionVisitor vis1;
    Face_index_map fimap1 = 
      mesh.add_property_map<face_descriptor, std::size_t>("f:i").first;
    Face_index_map fimap2 = 
      mesh2.add_property_map<face_descriptor, std::size_t>("f:i").first;


    EMesh3 imesh;
    const bool success = PMP::corefine_and_compute_union(
      mesh, mesh2, imesh,
      PMP::parameters::face_index_map(fimap1).visitor(vis1),
      PMP::parameters::face_index_map(fimap2)
    );
    if(!success) {
      Rcpp::stop("Union computation has failed.");
    }

    if(true) {
      //Rcpp::StringVector fcolors0(fcolors);
      const int nfaces1 = mcopy.number_of_faces();
      std::vector<Triangle> triangles1;
      triangles1.reserve(nfaces1);
      for(EMesh3::Face_index fd : mcopy.faces()) {
        auto vd = vertices_around_face(mcopy.halfedge(fd), mcopy).begin();
        triangles1.emplace_back(Triangle(
          mcopy.point(*(++vd)), mcopy.point(*(++vd)), mcopy.point(*vd)
        ));
      }
      const int nfaces2 = mesh2.number_of_faces();
      std::vector<Triangle> triangles2;
      triangles2.reserve(nfaces2);
      for(EMesh3::Face_index fd : mesh2.faces()) {
        auto vd = vertices_around_face(mesh2.halfedge(fd), mesh2).begin();
        triangles2.emplace_back(Triangle(
          mesh2.point(*(++vd)), mesh2.point(*(++vd)), mesh2.point(*vd)
        ));
      }
    
      const size_t nf = imesh.number_of_faces();
      //Rcpp::StringVector fcolors1(nf);
      //size_t j = 0;
      for(EMesh3::Face_index f : imesh.faces()) {
        auto vd = vertices_around_face(imesh.halfedge(f), imesh).begin();
        Triangle tr(imesh.point(*(++vd)), imesh.point(*(++vd)), imesh.point(*vd));
        EPoint3 c = CGAL::centroid(tr);
        int k1, k2;
        bool found1 = false, found2 = false;
        for(k1 = 0; k1 < nfaces1; k1++) {
          if(triangles1[k1].has_on(c)) {
            found1 = true;
            break;
          }
        }
        if(!found1) {
          k1 = -1;
        }
        for(k2 = 0; k2 < nfaces2; k2++) {
          if(triangles2[k2].has_on(c)) {
            found2 = true;
            break;
          }
        }
        if(!found2) {
          k2 = -1;
        }
        Rcpp::Rcout << f << ": " << k1 << ", " << k2 << "\n";

        //fcolors1(j++) = fcolors0(k);
      }
      //fcolors = Rcpp::Nullable<Rcpp::StringVector>(fcolors1);
    }

//    Rcpp::Rcout << "nfaces: " << imesh.number_of_faces() << "\n";

    Rcpp::Rcout << "PROPERTIES" << "\n";
    std::vector<std::string> props = imesh.properties<face_descriptor>();
    for(int p = 0; p < props.size(); p++) {
      Rcpp::Rcout << props[p] << "\n";
    }

    // Face_index_map fimap = mesh.property_map<face_descriptor, std::size_t>("f:i").first;
    // Rcpp::IntegerVector Fimap(mesh.number_of_faces());
    // int i = 0;
    // for(EMesh3::Face_index fi : mesh.faces()) {
    //   Fimap(i++) = fimap[fi];
    // }

    Rcpp::IntegerVector Fimap1(mesh.number_of_faces());
    int i1 = 0;
    for(EMesh3::Face_index fi : mesh.faces()) {
      Fimap1(i1++) = fimap1[fi];
    }
    Rcpp::IntegerVector Fimap2(mesh2.number_of_faces());
    int i2 = 0;
    for(EMesh3::Face_index fi : mesh2.faces()) {
      Fimap2(i2++) = fimap2[fi];
    }
    Rcpp::IntegerMatrix Fvis1((*(vis1.fmap)).size(), 2);
    int ii1 = 0;
    for(const auto& [fnew, fsplit] : *(vis1.fmap)) {
      Fvis1(ii1++, Rcpp::_) = Rcpp::IntegerVector({int(fsplit), int(fnew)});
    }

    Rcpp::List myimesh = Rcpp::List::create(
      //Rcpp::Named("xptr") = Rcpp::XPtr<EMesh3>(new EMesh3(imesh), false),
      Rcpp::Named("normals") = R_NilValue,
      Rcpp::Named("vcolors") = R_NilValue,
      Rcpp::Named("fcolors") = R_NilValue,
      //Rcpp::Named("fimap") = Fimap,
      Rcpp::Named("fimap1") = Fimap1,
      Rcpp::Named("fimap2") = Fimap2,
      Rcpp::Named("fvis1") = Fvis1
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