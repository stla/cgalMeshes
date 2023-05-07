#ifndef _HEADER_
#include "cgalMesh.h"
#endif

class CGALmesh {
public:
  EMesh3 mesh;
  Rcpp::XPtr<EMesh3> xptr;

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
      xptr(Rcpp::XPtr<EMesh3>(&mesh, false))
      {}
  
  CGALmesh(Rcpp::XPtr<EMesh3> xptr_)
    : mesh(*(xptr_.get())),
      xptr(Rcpp::XPtr<EMesh3>(&mesh, false))
      {}
  
  CGALmesh(const std::string filename, bool binary, bool soup)
    : mesh(readMeshFile(filename, binary, soup)), 
      xptr(Rcpp::XPtr<EMesh3>(&mesh, false))
      {}
  

  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  void assignNormals(Rcpp::NumericMatrix normals) {
    if(normals.ncol() != mesh.number_of_vertices()) {
      Rcpp::stop(
        "The number of normals does not match the number of vertices."
      );
    }
    removeProperties(mesh, {"v:normal"});
    Normals_map vnormal = 
      mesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
        "v:normal", defaultNormal()
      ).first;
    int i = 0;
    for(EMesh3::Vertex_index vi : mesh.vertices()) {
      Rcpp::NumericVector normal = normals(Rcpp::_, i++);
      vnormal[vi] = normal;
    }
  }


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  Rcpp::List boundingBox() {
    Bbox3 bbox = PMP::bbox(mesh);
    Rcpp::NumericVector lcorner = {bbox.xmin(), bbox.ymin(), bbox.zmin()};
    Rcpp::NumericVector ucorner = {bbox.xmax(), bbox.ymax(), bbox.zmax()};
    return Rcpp::List::create(
      Rcpp::Named("lcorner") = lcorner,
      Rcpp::Named("ucorner") = ucorner
    );
  }


  // ----------------------------------------------------------------------- //
  void CatmullClark(unsigned int iterations) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    removeProperties(
      mesh, {"v:normal", "v:scalar", "v:color", "f:scalar", "f:color"}
    );
    CGAL::Subdivision_method_3::CatmullClark_subdivision(
      mesh, CGAL::parameters::number_of_iterations(iterations)
    );
  }


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  Rcpp::List clipMesh(Rcpp::XPtr<EMesh3> clipperXPtr, const bool clipVolume) {
    EMesh3 clipper = *(clipperXPtr.get());
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(clipVolume) {
      if(PMP::does_self_intersect(mesh)) {
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
      copy_prop<face_descriptor, std::string>(mesh, "f:color");
    MaybeFcolorMap fcolorMap2_ = 
      copy_prop<face_descriptor, std::string>(clipper, "f:color");
    const bool hasColors = fcolorMap_.second && fcolorMap2_.second;
    MaybeFscalarMap fscalarMap_ = 
      copy_prop<face_descriptor, double>(mesh, "f:scalar");
    MaybeFscalarMap fscalarMap2_ = 
      copy_prop<face_descriptor, double>(clipper, "f:scalar");
    const bool hasScalars = fscalarMap_.second && fscalarMap2_.second;

    std::size_t nfaces = mesh.number_of_faces();
    std::size_t nfaces_clipper = clipper.number_of_faces();

    ClipVisitor vis; 
    
    MapBetweenFaces mesh_map;
    MapBetweenFaces clipper_map;
    for(std::size_t i = 0; i < nfaces; i++) {
      face_descriptor fi(i);
      mesh_map[fi] = fi;
    }
    for(std::size_t i = 0; i < nfaces_clipper; i++) {
      face_descriptor fi(i);
      clipper_map[fi] = fi;
    }
    vis.FACEMAPS[&mesh]    = &mesh_map;
    vis.FACEMAPS[&clipper] = &clipper_map;
    
    std::size_t undetermined = 999999999;
    Face_index_map fimap = 
      mesh.add_property_map<face_descriptor, std::size_t>(
        "f:i", undetermined
      ).first;
    
    const bool doNotModify = !clipVolume;
    const bool clipping = PMP::clip(
      mesh, clipper,
      PMP::parameters::clip_volume(clipVolume)
        .visitor(vis).face_index_map(fimap),
      PMP::parameters::clip_volume(clipVolume).do_not_modify(doNotModify)
    );
    if(!clipping) {
      Rcpp::stop("Clipping has failed.");
    }

    mesh.collect_garbage();

    /* --------------- clipVolume is false --------------- */
    if(!clipVolume) {
      if(hasColors || hasScalars) {
        std::map<face_descriptor, std::string> fcolorMap;
        Fcolors_map newfcolor;
        std::map<face_descriptor, double> fscalarMap;
        Fscalars_map newfscalar;
        if(hasColors) {
          fcolorMap = fcolorMap_.first;
          newfcolor = 
            mesh.add_property_map<face_descriptor, std::string>(
              "f:color", ""
            ).first;
        }
        if(hasScalars) {
          fscalarMap = fscalarMap_.first;
          newfscalar = 
            mesh.add_property_map<face_descriptor, double>(
              "f:scalar", nan("")
            ).first;
        }
        for(EMesh3::Face_index fi : mesh.faces()) {
          face_descriptor fd = CGAL::SM_Face_index(fimap[fi]);
          face_descriptor fdnew = mesh_map[fd];
          if(hasColors) {
            newfcolor[fi] = fcolorMap[fdnew];
          }
          if(hasScalars) {
            newfscalar[fi] = fscalarMap[fdnew];
          }
        }
      }
      mesh.remove_property_map(fimap);
      return Rcpp::List::create();
    }

    /* --------------- clipVolume is true --------------- */
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
      for(auto it = ftargets.begin(); it != ftargets.end(); ++it) {
        face_descriptor fd = clipper_map[it->second];
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
      mesh.add_property_map<face_descriptor, std::size_t>("f:which").first;

    Fcolors_map newfcolor;
    Fscalars_map newfscalar;
    if(hasColors) {
      newfcolor = mesh.add_property_map<face_descriptor, std::string>(
        "f:color", ""
      ).first;
    }
    if(hasScalars) {
      newfscalar = mesh.add_property_map<face_descriptor, double>(
        "f:scalar", nan("")
      ).first;
    }

    for(EMesh3::Face_index fi : mesh.faces()) {
      std::size_t ifi = fimap[fi];
      face_descriptor fd = mesh_map[CGAL::SM_Face_index(ifi)];
      if(ifi != undetermined) {
        whichPart[fi] = 0;
        if(hasColors) {
          newfcolor[fi] = fcolorMap[fd]; 
        }
        if(hasScalars) {
          newfscalar[fi] = fscalarMap[fd];
        }
      } else {
        whichPart[fi] = 1;
        if(hasColors) {
          newfcolor[fi] = fcolorMap_clipper[zeros[fi]];
        }
        if(hasScalars) {
          newfscalar[fi] = fscalarMap_clipper[zeros[fi]];
        }
      } 
    }

    mesh.remove_property_map(fimap);

    EMesh3 tmesh;
    {
      Filtered_graph ffg(mesh, 0, whichPart);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, tmesh, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(
        mesh, tmesh, f2fmap_, "f:color"
      );
      copy_property<ffg_face_descriptor, face_descriptor, double>(
        mesh, tmesh, f2fmap_, "f:scalar"
      );
    }

    EMesh3 tmesh2;
    {
      Filtered_graph ffg(mesh, 1, whichPart);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, tmesh2, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(
        mesh, tmesh2, f2fmap_, "f:color"
      );
      copy_property<ffg_face_descriptor, face_descriptor, double>(
        mesh, tmesh2, f2fmap_, "f:scalar"
      );
    }

    Rcpp::List meshes = Rcpp::List::create(
      Rcpp::Named("mesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(tmesh), false),
      Rcpp::Named("mesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(tmesh2), false)
    );

    return meshes;
  }


  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  Rcpp::List clipToPlane(
    Rcpp::NumericVector planePoint, 
    Rcpp::NumericVector planeNormal, 
    const bool clipVolume
  ) {
    EPoint3 point(planePoint[0], planePoint[1], planePoint[2]);
    EVector3 normal(planeNormal[0], planeNormal[1], planeNormal[2]);
    EPlane3 plane(point, normal);
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(clipVolume && PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The mesh self-intersects.");
    }

    std::size_t nfaces = mesh.number_of_faces();

    MaybeFcolorMap fcolorMap_ = 
      copy_prop<face_descriptor, std::string>(mesh, "f:color");
    const bool hasColors = fcolorMap_.second;
    MaybeFscalarMap fscalarMap_ = 
      copy_prop<face_descriptor, double>(mesh, "f:scalar");
    const bool hasScalars = fscalarMap_.second;

    ClipPlaneVisitor vis;
    MapBetweenFaces mesh_map;
    for(std::size_t i = 0; i < nfaces; i++) {
      face_descriptor fi(i);
      mesh_map[fi] = fi;
    }
    vis.FACEMAP = &mesh_map;

    Face_index_map fdmap = 
      mesh.add_property_map<face_descriptor, std::size_t>(
          "f:dummy", 0
      ).first;
    
    const bool clipping = PMP::clip(
      mesh, plane,
      PMP::parameters::clip_volume(clipVolume)
                      .visitor(vis)
                      .allow_self_intersections(!clipVolume)
    );
    if(!clipping) {
      Rcpp::stop("Clipping has failed.");
    }

    mesh.remove_property_map(fdmap);
    
    Face_index_map fimap = 
      mesh.add_property_map<face_descriptor, std::size_t>(
        "f:i", 999999999
      ).first;
    for(face_descriptor fd : mesh.faces()) {
      fimap[fd] = std::size_t(fd);
    }

    mesh.collect_garbage();

    /* --------------- clipVolume is false --------------- */
    if(!clipVolume){

      if(hasColors || hasScalars) {
        std::map<face_descriptor, std::string> fcolorMap;
        Fcolors_map newfcolor;
        std::map<face_descriptor, double> fscalarMap;
        Fscalars_map newfscalar;
        if(hasColors) {
          fcolorMap = fcolorMap_.first;
          newfcolor = 
            mesh.add_property_map<face_descriptor, std::string>(
              "f:color", ""
            ).first;
        }
        if(hasScalars) {
          fscalarMap = fscalarMap_.first;
          newfscalar = 
            mesh.add_property_map<face_descriptor, double>(
              "f:scalar", nan("")
            ).first;
        }

        for(EMesh3::Face_index fi : mesh.faces()) {
          std::size_t ffi = fimap[fi];
          face_descriptor fd = mesh_map[CGAL::SM_Face_index(ffi)];
          if(hasColors) {
            newfcolor[fi] = fcolorMap[fd];
          }
          if(hasScalars) {
            newfscalar[fi] = fscalarMap[fd];
          }
        }
      }

      mesh.remove_property_map(fimap);
      return Rcpp::List::create();
    }

    
    /* --------------- clipVolume is true --------------- */
    MapBetweenFaces ftargets = *(vis.ftargets);

    std::map<face_descriptor, std::string> fcolorMap;
    Fcolors_map newfcolor;
    std::map<face_descriptor, double> fscalarMap;
    Fscalars_map newfscalar;
    if(hasColors) {
      fcolorMap = fcolorMap_.first;
      newfcolor = 
        mesh.add_property_map<face_descriptor, std::string>(
          "f:color", ""
        ).first;
    }
    if(hasScalars) {
      fscalarMap = fscalarMap_.first;
      newfscalar = 
        mesh.add_property_map<face_descriptor, double>(
          "f:scalar", nan("")
        ).first;
    }

    Face_index_map fwhich = 
      mesh.add_property_map<face_descriptor, std::size_t>("f:which", 1).first;

    for(EMesh3::Face_index fi : mesh.faces()) {
      std::size_t ffi = fimap[fi];
      face_descriptor fd = CGAL::SM_Face_index(ffi);
      if(auto search = ftargets.find(fd); search != ftargets.end()) {
        fwhich[fi] = 2;
        ftargets.erase(fd);
      } else {
        fd = mesh_map[fd];
        if(hasColors) {
          newfcolor[fi] = fcolorMap[fd];
        }
        if(hasScalars) {
          newfscalar[fi] = fscalarMap[fd];
        }
      }
    }

    EMesh3 mesh1;
    {
      Filtered_graph ffg(mesh, 1, fwhich);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> 
        f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, mesh1, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(
        mesh, mesh1, f2fmap_, "f:color"
      );
      copy_property<ffg_face_descriptor, face_descriptor, double>(
        mesh, mesh1, f2fmap_, "f:scalar"
      );
    }

    EMesh3 mesh2;
    {
      Filtered_graph ffg(mesh, 2, fwhich);
      CGAL::copy_face_graph(ffg, mesh2);
    }

    mesh.remove_property_map(fimap);
    mesh.remove_property_map(fwhich);

    return Rcpp::List::create(
      Rcpp::Named("mesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(mesh1), false),
      Rcpp::Named("mesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(mesh2), false)
    );
  }


  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  Rcpp::List clipToIsoCuboid(
    Rcpp::NumericVector lcorner, 
    Rcpp::NumericVector ucorner, 
    const bool clipVolume
  ) {
    EPoint3 lpoint(lcorner[0], lcorner[1], lcorner[2]);
    EPoint3 upoint(ucorner[0], ucorner[1], ucorner[2]);
    IsoCuboid3 isocuboid(lpoint, upoint, 0);
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(clipVolume && PMP::does_self_intersect(mesh)) {
      Rcpp::stop("The mesh self-intersects.");
    }

    std::size_t nfaces = mesh.number_of_faces();

    MaybeFcolorMap fcolorMap_ = 
      copy_prop<face_descriptor, std::string>(mesh, "f:color");
    const bool hasColors = fcolorMap_.second;
    MaybeFscalarMap fscalarMap_ = 
      copy_prop<face_descriptor, double>(mesh, "f:scalar");
    const bool hasScalars = fscalarMap_.second;

    ClipVisitor vis;
    MapBetweenFaces mesh_map;
    for(std::size_t i = 0; i < nfaces; i++) {
      face_descriptor fi(i);
      mesh_map[fi] = fi;
    }
    vis.FACEMAPS[&mesh] = &mesh_map;

    const bool clipping = PMP::clip(
      mesh, isocuboid,
      PMP::parameters::clip_volume(clipVolume)
                      .visitor(vis)
                      .allow_self_intersections(!clipVolume)
    );
    if(!clipping) {
      Rcpp::stop("Clipping has failed.");
    }

    Face_index_map fimap = 
      mesh.add_property_map<face_descriptor, std::size_t>(
        "f:i", 999999999
      ).first;
    for(face_descriptor fd : mesh.faces()) {
      fimap[fd] = std::size_t(fd);
    }

    mesh.collect_garbage();

    /* --------------- clipVolume is false --------------- */
    if(!clipVolume){

      if(hasColors || hasScalars) {
        std::map<face_descriptor, std::string> fcolorMap;
        Fcolors_map newfcolor;
        std::map<face_descriptor, double> fscalarMap;
        Fscalars_map newfscalar;
        if(hasColors) {
          fcolorMap = fcolorMap_.first;
          newfcolor = 
            mesh.add_property_map<face_descriptor, std::string>(
              "f:color", ""
            ).first;
        }
        if(hasScalars) {
          fscalarMap = fscalarMap_.first;
          newfscalar = 
            mesh.add_property_map<face_descriptor, double>(
              "f:scalar", nan("")
            ).first;
        }

        for(EMesh3::Face_index fi : mesh.faces()) {
          std::size_t ffi = fimap[fi];
          face_descriptor fd = mesh_map[CGAL::SM_Face_index(ffi)];
          if(hasColors) {
            newfcolor[fi] = fcolorMap[fd];
          }
          if(hasScalars) {
            newfscalar[fi] = fscalarMap[fd];
          }
        }
      }

      mesh.remove_property_map(fimap);
      return Rcpp::List::create();
    }

    /* --------------- clipVolume is true --------------- */
    MapBetweenFaces ftargets = *(vis.ftargets);

    std::map<face_descriptor, std::string> fcolorMap;
    Fcolors_map newfcolor;
    std::map<face_descriptor, double> fscalarMap;
    Fscalars_map newfscalar;
    if(hasColors) {
      fcolorMap = fcolorMap_.first;
      newfcolor = 
        mesh.add_property_map<face_descriptor, std::string>(
          "f:color", ""
        ).first;
    }
    if(hasScalars) {
      fscalarMap = fscalarMap_.first;
      newfscalar = 
        mesh.add_property_map<face_descriptor, double>(
          "f:scalar", nan("")
        ).first;
    }

    Face_index_map fwhich = 
      mesh.add_property_map<face_descriptor, std::size_t>("f:which", 1).first;

    for(EMesh3::Face_index fi : mesh.faces()) {
      std::size_t ffi = fimap[fi];
      face_descriptor fd = CGAL::SM_Face_index(ffi);
      if(auto search = ftargets.find(fd); search != ftargets.end()) {
        fwhich[fi] = 2;
        ftargets.erase(fd);
      } else {
        fd = mesh_map[fd];
        if(hasColors) {
          newfcolor[fi] = fcolorMap[fd];
        }
        if(hasScalars) {
          newfscalar[fi] = fscalarMap[fd];
        }
      }
    }

    EMesh3 mesh1;
    {
      Filtered_graph ffg(mesh, 1, fwhich);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, mesh1, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(
        mesh, mesh1, f2fmap_, "f:color"
      );
      copy_property<ffg_face_descriptor, face_descriptor, double>(
        mesh, mesh1, f2fmap_, "f:scalar"
      );
    }

    EMesh3 mesh2;
    {
      Filtered_graph ffg(mesh, 2, fwhich);
      CGAL::copy_face_graph(ffg, mesh2);
    }

    mesh.remove_property_map(fimap);
    mesh.remove_property_map(fwhich);

    return Rcpp::List::create(
      Rcpp::Named("mesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(mesh1), false),
      Rcpp::Named("mesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(mesh2), false)
    );
  }

  
  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  Rcpp::XPtr<EMesh3> clone() {
    EMesh3 copy = cloneMesh(
      mesh, {"f:color", "v:color", "f:scalar", "v:scalar", "v:normal"}
    );
    return Rcpp::XPtr<EMesh3>(new EMesh3(copy), false);
  }


  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  void collectGarbage() {
    Rcpp::Rcout << "Mesh has garbage: " << mesh.has_garbage() << ".\n";
    mesh.collect_garbage();
  }


  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  Rcpp::List connectedComponents(const bool triangulate) {

    EMesh3 tmesh = cloneMesh(
      mesh, {"f:color", "v:color", "f:scalar", "v:scalar", "v:normal"}
    );

    const bool really_triangulate = 
      triangulate && !CGAL::is_triangle_mesh(mesh);
    if(really_triangulate) {
      triangulateMesh(tmesh);
    }

    Face_index_map fccmap = 
      tmesh.add_property_map<face_descriptor, std::size_t>("f:CC", 0).first;

    const std::size_t ncc = PMP::connected_components(tmesh, fccmap);

    if(ncc == 1) {
      Message("Only one component found.\n");
    } else {
      const std::string msg = "Found " + std::to_string(ncc) + " components.\n";
      Message(msg);
    }

    Rcpp::List xptrs(ncc);

    for(std::size_t c = 0; c < ncc; c++) {

      Filtered_graph ffg(tmesh, c, fccmap);
      MapBetweenVertexDescriptors v2vmap_;
      boost::associative_property_map<MapBetweenVertexDescriptors> 
        v2vmap(v2vmap_);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> 
        f2fmap(f2fmap_);

      EMesh3 cmesh;
      CGAL::copy_face_graph(
        ffg, cmesh, 
        CGAL::parameters::vertex_to_vertex_map(v2vmap).face_to_face_map(f2fmap)
      );

      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, Rcpp::NumericVector
      >(tmesh, cmesh, v2vmap_, "v:normal");
      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, std::string
      >(tmesh, cmesh, v2vmap_, "v:color");
      copy_property<
        ffg_vertex_descriptor, vertex_descriptor, double
      >(tmesh, cmesh, v2vmap_, "v:scalar");
      copy_property<
        ffg_face_descriptor, face_descriptor, std::string
      >(tmesh, cmesh, f2fmap_, "f:color");
      copy_property<
        ffg_face_descriptor, face_descriptor, double
      >(tmesh, cmesh, f2fmap_, "f:scalar");

      xptrs(c) = Rcpp::XPtr<EMesh3>(new EMesh3(cmesh), false);
    }

    return xptrs;
  }


  // ----------------------------------------------------------------------- //
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
  
  
  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  bool doesBoundVolume() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    return PMP::does_bound_a_volume(mesh);
  }


  // ----------------------------------------------------------------------- //
  bool doesSelfIntersect() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    return PMP::does_self_intersect(mesh);
  }


  // ----------------------------------------------------------------------- //
  void DooSabin(unsigned int iterations) {
    removeProperties(
      mesh, {"v:normal", "v:scalar", "v:color", "f:scalar", "f:color"}
    );
    CGAL::Subdivision_method_3::DooSabin_subdivision(
      mesh, CGAL::parameters::number_of_iterations(iterations)
    );
    mesh.collect_garbage();
  }


  // ----------------------------------------------------------------------- //
  Rcpp::XPtr<EMesh3> dual() {
    EMesh3 dualmesh = dualMesh(mesh);
    return Rcpp::XPtr<EMesh3>(new EMesh3(dualmesh), false);
  }


  // ----------------------------------------------------------------------- //
  Rcpp::DataFrame edges() {
    return getEdges<EK, EMesh3, EPoint3>(mesh);
  }


  // ----------------------------------------------------------------------- //
  Rcpp::IntegerVector facesAroundVertex(int v) {
    if(v >= mesh.number_of_vertices()) {
      Rcpp::stop("Too high vertex index.");
    }
    vertex_descriptor vd = CGAL::SM_Vertex_index(v);
    Rcpp::IntegerVector faces;
    for(face_descriptor fd : 
          CGAL::faces_around_target(mesh.halfedge(vd), mesh)) {
      if(fd != EMesh3::null_face()) {
        faces.push_back(int(fd) + 1);
      }
    }
    return faces;
  }


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  Rcpp::XPtr<EMesh3> fillBoundaryHole(int border, bool fairhole) {
    std::vector<halfedge_descriptor> border_cycles;
    PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
    const int nborders = border_cycles.size();
    if(nborders == 0) {
      Rcpp::stop("There's no border in this mesh.");
    } 
    if(border >= nborders) {
      std::string msg;
      if(nborders == 1) {
        msg = "There's only one border in this mesh.";
      } else {
        msg = "There are only " + std::to_string(nborders) + " borders.";
      }
      Rcpp::stop(msg);
    }
    halfedge_descriptor h = border_cycles[border];
    std::vector<face_descriptor>   patch_faces;
    std::vector<vertex_descriptor> patch_vertices;
    if(fairhole) {
      const bool success = std::get<0>(
        PMP::triangulate_refine_and_fair_hole(
          mesh, h, 
          std::back_inserter(patch_faces), std::back_inserter(patch_vertices)
        )
      );
      if(!success) {
        Message("Fairing failed.");
      }
    } else {
      PMP::triangulate_and_refine_hole(
        mesh, h, 
        std::back_inserter(patch_faces), std::back_inserter(patch_vertices)
      );
    }
    Face_index_map fimap = mesh.add_property_map<face_descriptor, std::size_t>(
      "f:i", 0
    ).first;
    for(int i = 0; i < patch_faces.size(); i++) {
      fimap[patch_faces[i]] = 2;
    }
    EMesh3 hole;
    Filtered_graph ffg(mesh, 2, fimap);
    CGAL::copy_face_graph(ffg, hole);
    mesh.remove_property_map(fimap);
    return Rcpp::XPtr<EMesh3>(new EMesh3(hole), false);
  }


  // ----------------------------------------------------------------------- //
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
      const bool valid = ffg.is_selection_valid();
      if(!valid) {
        Rcpp::warning("Selection is possibly invalid.");
      }
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
      const bool valid = ffg.is_selection_valid();
      if(!valid) {
        Rcpp::warning("Selection is possibly invalid.");
      }
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


  // ----------------------------------------------------------------------- //
  void fixManifoldness() {
    std::size_t nmv = PMP::duplicate_non_manifold_vertices(mesh);
    if(nmv > 0) {
      std::string msg;
      if(nmv == 1) {
        msg = "One non-manifold vertex duplicated.";
      } else {
        msg = "Duplicated " + std::to_string(nmv) + " non-manifold vertices.";
      }
      Message(msg);
    } else {
      Message("No non-manifold vertex has been found.");
    }
  }


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  Rcpp::NumericMatrix getFacesInfo() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    Rcpp::CharacterVector rownames = {"cx", "cy", "cz", "area"};
    Rcpp::NumericMatrix FacesInfo(4, mesh.number_of_faces());
    int i = 0;
    for(face_descriptor fd : mesh.faces()) {
      auto vs = vertices_around_face(mesh.halfedge(fd), mesh).begin();
      EPoint3 p1 = mesh.point(*(vs++));
      EPoint3 p2 = mesh.point(*(vs++));
      EPoint3 p3 = mesh.point(*vs);
      EPoint3 centroid = CGAL::centroid(p1, p2, p3);
      EK::FT sarea = CGAL::squared_area(p1, p2, p3);
      Rcpp::NumericVector col_i = {
        CGAL::to_double<EK::FT>(centroid.x()),
        CGAL::to_double<EK::FT>(centroid.y()),
        CGAL::to_double<EK::FT>(centroid.z()),
        sqrt(CGAL::to_double<EK::FT>(sarea))
      };
      FacesInfo(Rcpp::_, i++) = col_i;
    }
    Rcpp::rownames(FacesInfo) = rownames;
    return Rcpp::transpose(FacesInfo);
  }


 // ----------------------------------------------------------------------- //
 Rcpp::List getFacesList() {
    return getFaces<EMesh3>(mesh);
  }


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  Rcpp::NumericMatrix getVertices() {
    return getVertices_EK(mesh);
  }  


  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  double HausdorffApproximate(Rcpp::XPtr<EMesh3> mesh2XPtr, bool symmetric) {
    if(CGAL::is_empty(mesh)) {
      Rcpp::stop("The reference mesh is empty.");
    }
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The reference mesh is not triangle.");
    }
    EMesh3 mesh2 = *(mesh2XPtr.get());
    if(CGAL::is_empty(mesh2)) {
      Rcpp::stop("The second mesh is empty.");
    }
    if(!CGAL::is_triangle_mesh(mesh2)) {
      Rcpp::stop("The second mesh is not triangle.");
    }
    double d;
    if(symmetric) {
      d = PMP::approximate_symmetric_Hausdorff_distance<PIA_TAG>(mesh, mesh2);
    } else {
      d = PMP::approximate_Hausdorff_distance<PIA_TAG>(mesh, mesh2);
    }
    return d;
  }
  
  
  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  double HausdorffEstimate(
      Rcpp::XPtr<EMesh3> mesh2XPtr, double errorBound, bool symmetric
  ) {
    if(CGAL::is_empty(mesh)) {
      Rcpp::stop("The reference mesh is empty.");
    }
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The reference mesh is not triangle.");
    }
    EMesh3 mesh2 = *(mesh2XPtr.get());
    if(CGAL::is_empty(mesh2)) {
      Rcpp::stop("The second mesh is empty.");
    }
    if(!CGAL::is_triangle_mesh(mesh2)) {
      Rcpp::stop("The second mesh is not triangle.");
    }
    double d;
    if(symmetric) {
      d = PMP::bounded_error_symmetric_Hausdorff_distance<PIA_TAG>(
        mesh, mesh2, errorBound
      );
    } else {
      d = PMP::bounded_error_Hausdorff_distance<PIA_TAG>(mesh, mesh2, errorBound);
    }
    return d;
  }
  
  
  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
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


  // ----------------------------------------------------------------------- //
  bool isClosed() {
    return CGAL::is_closed(mesh);
  }


  // ----------------------------------------------------------------------- //
  void isotropicRemeshing(
    const double targetEdgeLength, 
    const unsigned niters, const unsigned nrelaxsteps 
  ) {
    std::vector<halfedge_descriptor> borderHalfedges;
    PMP::border_halfedges(
      mesh.faces(), mesh, std::back_inserter(borderHalfedges)
    );
    std::vector<edge_descriptor> border;
    int nhborder = borderHalfedges.size();
    border.reserve(nhborder);
    for(int i = 0; i < nhborder; i++) {
      border.emplace_back(mesh.edge(borderHalfedges[i]));
    }
    PMP::split_long_edges(border, targetEdgeLength, mesh);
    PMP::isotropic_remeshing(
      mesh.faces(), targetEdgeLength, mesh,
      PMP::parameters::number_of_iterations(niters)
                      .number_of_relaxation_steps(nrelaxsteps)
                      .protect_constraints(true)
    );
    mesh.collect_garbage();
  }


  // ----------------------------------------------------------------------- //
  bool isOutwardOriented() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    return PMP::is_outward_oriented(mesh);
  }


  // ----------------------------------------------------------------------- //
  bool isQuad() {
    return CGAL::is_quad_mesh(mesh);
  }

  // ----------------------------------------------------------------------- //
  bool isTriangle() {
    return CGAL::is_triangle_mesh(mesh);
  }


  // ----------------------------------------------------------------------- //
  bool isValid() {
    return mesh.is_valid();
  }


  // ----------------------------------------------------------------------- //
  bool isValidFaceGraph() {
    return CGAL::is_valid_face_graph(mesh, true);
  }


  // ----------------------------------------------------------------------- //
  bool isValidHalfedgeGraph() {
    return CGAL::is_valid_halfedge_graph(mesh, true);
  }


  // ----------------------------------------------------------------------- //
  bool isValidPolygonMesh() {
    return CGAL::is_valid_polygon_mesh(mesh, true);
  }


  // ----------------------------------------------------------------------- //
  void LoopSubdivision(unsigned int iterations) {
    removeProperties(
      mesh, {"v:normal", "v:scalar", "v:color", "f:scalar", "f:color"}
    );
    CGAL::Subdivision_method_3::Loop_subdivision(
      mesh, CGAL::parameters::number_of_iterations(iterations)
    );
  }


  // ----------------------------------------------------------------------- //
  void merge(Rcpp::XPtr<EMesh3> mesh2XPtr) {
    EMesh3 mesh2 = *(mesh2XPtr.get());
    mesh += mesh2;
  }


  // ----------------------------------------------------------------------- //
  void orientToBoundVolume() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    PMP::orient_to_bound_a_volume(mesh);
    // faut-il updater normals?
  }


  // ----------------------------------------------------------------------- //
  void print() {
    Rcpp::Rcout << "Mesh with " << mesh.number_of_vertices() 
                << " vertices and " << mesh.number_of_faces() << " faces.\n";
    if(CGAL::is_triangle_mesh(mesh)) {
      Rcpp::Rcout << "The mesh is triangle.\n";
    } else {
      Rcpp::Rcout << "The mesh is not triangle.\n";
    }
  }


  // ----------------------------------------------------------------------- //
  void removeSelfIntersections() {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    PMP::experimental::remove_self_intersections(mesh);
    mesh.collect_garbage();
  }


  // ----------------------------------------------------------------------- //
  void reverseFaceOrientations() {
    PMP::reverse_face_orientations(mesh);
    // update normals ?
  }

  
  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  Rcpp::NumericMatrix sampleMesh(const unsigned nsims) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    std::vector<EPoint3> sims;
    PMP::sample_triangle_mesh(
      mesh, std::back_inserter(sims),
      PMP::parameters::number_of_points_on_faces(nsims)
    );
    Rcpp::NumericMatrix rsims = points3_to_matrix(sims);
    return Rcpp::transpose(rsims);
  }  
  

  // ----------------------------------------------------------------------- //
  // ----------------------------------------------------------------------- //
  Rcpp::IntegerMatrix sharpEdges(double angleBound) {
    EK::FT angle_bound(angleBound);
    EMesh3::Property_map<edge_descriptor, bool> edmap = 
      mesh.add_property_map<edge_descriptor, bool>("e:bool", false).first;
    PMP::detect_sharp_edges(mesh, angle_bound, edmap);
    std::vector<int> vedges;
    int n = 0;
    for(edge_descriptor ed : mesh.edges()) {
      if(edmap[ed]) {
        n++;
        vedges.push_back(int(ed) + 1);
        vedges.push_back(int(source(ed, mesh)) + 1);
        vedges.push_back(int(target(ed, mesh)) + 1);
      }
    }
    mesh.remove_property_map(edmap);
    Rcpp::IntegerVector Vedges(vedges.begin(), vedges.end());
    Rcpp::CharacterVector colnames = {"edge", "v1", "v2"};
    Vedges.attr("dim") = Rcpp::Dimension(3, n);
    Rcpp::IntegerMatrix Edges = 
      Rcpp::transpose(Rcpp::as<Rcpp::IntegerMatrix>(Vedges));
    Rcpp::colnames(Edges) = colnames;
    return Edges;
  }


  // ----------------------------------------------------------------------- //
  void Sqrt3Subdivision(unsigned int iterations) {
    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    removeProperties(
      mesh, {"v:normal", "v:scalar", "v:color", "f:scalar", "f:color"}
    );
    CGAL::Subdivision_method_3::Sqrt3_subdivision(
      mesh, CGAL::parameters::number_of_iterations(iterations)
    );
  }


  // ----------------------------------------------------------------------- //
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

    Face_index_map fimap_tmp = 
      mesh.add_property_map<face_descriptor, std::size_t>("f:i", 0).first;
    DifferenceVisitor vis;
    EMesh3 dmesh;
    const bool success = PMP::corefine_and_compute_difference(
      mesh, mesh2, dmesh,
      CGAL::parameters::visitor(vis)
    );
    if(!success) {
      Rcpp::stop("Difference computation has failed.");
    }
    mesh.remove_property_map(fimap_tmp);

    MapBetweenFaces fmap_mesh1      = *(vis.fmap_mesh1);
    MapBetweenFaces fmap_mesh2      = *(vis.fmap_mesh2);
    MapBetweenFaces fmap_difference = *(vis.fmap_difference);
    int nfaces_dmesh1 = *(vis.nfaces_dmesh1);

    const bool disjoint = fmap_mesh1.size() == 0 && fmap_mesh2.size() == 0;
    if(disjoint) {
      return Rcpp::List::create();
    }

    Fcolors_map fcolor;
    Fscalars_map fscalar;
    if(hasColors) {
      fcolor = dmesh.add_property_map<face_descriptor, std::string>(
        "f:color", ""
      ).first;
    }
    if(hasScalars) {
      fscalar = dmesh.add_property_map<face_descriptor, double>(
        "f:scalar", nan("")
      ).first;
    }

    Face_index_map fwhich = 
      dmesh.add_property_map<face_descriptor, std::size_t>(
        "f:which", 2
      ).first;
    if(hasColors || hasScalars) {
      for(int i = 0; i < nfaces_dmesh1; i++) {
        face_descriptor fi = CGAL::SM_Face_index(i);
        face_descriptor fd = fmap_difference[fi];
        face_descriptor fd1 = int(fd) < nfaces1 ? fd : fmap_mesh1[fd];
        if(hasColors) {
          fcolor[fi] = fcolor1[fd1];
        }
        if(hasScalars) {
          fscalar[fi] = fscalar1[fd1];
        }
        fwhich[fi] = 1;
      }
      for(int i = nfaces_dmesh1; i < dmesh.number_of_faces(); i++) {
        face_descriptor fi = CGAL::SM_Face_index(i);
        face_descriptor fd = fmap_difference[fi];
        face_descriptor fd2 = int(fd) < nfaces2 ? fd : fmap_mesh2[fd];
        if(hasColors) {
          fcolor[fi] = fcolor2[fd2];
        }
        if(hasScalars) {
          fscalar[fi] = fscalar2[fd2];
        }
      }
    } else {
      for(int i = 0; i < nfaces_dmesh1; i++) {
        fwhich[CGAL::SM_Face_index(i)] = 1;
      }
    }

    EMesh3 dmesh1;
    {
      Filtered_graph ffg(dmesh, 1, fwhich);
      //Rcpp::Rcout << "valid selection: " << ffg.is_selection_valid() << "\n";
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> 
        f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, dmesh1, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(
        dmesh, dmesh1, f2fmap_, "f:color"
      );
      copy_property<ffg_face_descriptor, face_descriptor, double>(
        dmesh, dmesh1, f2fmap_, "f:scalar"
      );
    }

    EMesh3 dmesh2;
    {
      Filtered_graph ffg(dmesh, 2, fwhich);
      //Rcpp::Rcout << "valid selection: " << ffg.is_selection_valid() << "\n";
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> 
        f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, dmesh2, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(
        dmesh, dmesh2, f2fmap_, "f:color"
      );
      copy_property<ffg_face_descriptor, face_descriptor, double>(
        dmesh, dmesh2, f2fmap_, "f:scalar"
      );
    }

    dmesh.remove_property_map(fwhich);

    return Rcpp::List::create(
      Rcpp::Named("dmesh") = Rcpp::XPtr<EMesh3>(new EMesh3(dmesh), false),
      Rcpp::Named("mesh1") = Rcpp::XPtr<EMesh3>(new EMesh3(dmesh1), false),
      Rcpp::Named("mesh2") = Rcpp::XPtr<EMesh3>(new EMesh3(dmesh2), false)
    );

  }


  // ----------------------------------------------------------------------- //
  void triangulate() {
    triangulateMesh(mesh);
  }


  // ----------------------------------------------------------------------- //
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

    Face_index_map fimap_tmp = 
      mesh.add_property_map<face_descriptor, std::size_t>("f:i", 0).first;
    UnionVisitor vis;
    EMesh3 umesh;
    const bool success = PMP::corefine_and_compute_union(
      mesh, mesh2, umesh,
      PMP::parameters::visitor(vis)
    );
    if(!success) {
      Rcpp::stop("Union computation has failed.");
    }
    mesh.remove_property_map(fimap_tmp);

    MapBetweenFaces fmap_mesh1 = *(vis.fmap_mesh1);
    MapBetweenFaces fmap_mesh2 = *(vis.fmap_mesh2);
    MapBetweenFaces fmap_union = *(vis.fmap_union);
    const bool disjoint = fmap_mesh1.size() == 0 && fmap_mesh2.size() == 0;

    if(disjoint) {
      std::map<vertex_descriptor, vertex_descriptor> vmap_union = 
        *(vis.vmap_union);
      // vertex normals
      std::pair<Normals_map, bool> vnormal1_ = 
        mesh.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
      std::pair<Normals_map, bool> vnormal2_ = 
        mesh2.property_map<vertex_descriptor, Rcpp::NumericVector>("v:normal");
      int nverts_mesh = mesh.number_of_vertices();
      int nverts_umesh = umesh.number_of_vertices();
      if(vnormal1_.second && vnormal2_.second) {
        Normals_map vnormal1 = vnormal1_.first;
        Normals_map vnormal2 = vnormal2_.first;
        Normals_map vnormal = 
          umesh.add_property_map<vertex_descriptor, Rcpp::NumericVector>(
          "v:normal", defaultNormal()
        ).first;
        for(int i = 0; i < nverts_mesh; i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vnormal[vi] = vnormal1[vd];
        }
        for(int i = nverts_mesh; i < nverts_umesh; i++) {
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
        Vcolors_map vcolor = 
          umesh.add_property_map<vertex_descriptor, std::string>(
            "v:color", ""
          ).first;
        for(int i = 0; i < nverts_mesh; i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vcolor[vi] = vcolor1[vd];
        }
        for(int i = nverts_mesh; i < nverts_umesh; i++) {
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
        for(int i = 0; i < nverts_mesh; i++) {
          vertex_descriptor vi = CGAL::SM_Vertex_index(i);
          vertex_descriptor vd = vmap_union[vi];
          vscalar[vi] = vscalar1[vd];
        }
        for(int i = nverts_mesh; i < nverts_umesh; i++) {
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

    int nfaces_umesh1 = *(vis.nfaces_umesh1);
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

    Face_index_map fwhich = 
      umesh.add_property_map<face_descriptor, std::size_t>(
        "f:which", 0
      ).first;
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

    EMesh3 umesh1;
    {
      Filtered_graph ffg(umesh, 1, fwhich);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> 
        f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, umesh1, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(
        umesh, umesh1, f2fmap_, "f:color"
      );
      copy_property<ffg_face_descriptor, face_descriptor, double>(
        umesh, umesh1, f2fmap_, "f:scalar"
      );
    }

    EMesh3 umesh2;
    {
      Filtered_graph ffg(umesh, 2, fwhich);
      MapBetweenFaceDescriptors f2fmap_;
      boost::associative_property_map<MapBetweenFaceDescriptors> 
        f2fmap(f2fmap_);
      CGAL::copy_face_graph(
        ffg, umesh2, CGAL::parameters::face_to_face_map(f2fmap)
      );
      copy_property<ffg_face_descriptor, face_descriptor, std::string>(
        umesh, umesh2, f2fmap_, "f:color"
      );
      copy_property<ffg_face_descriptor, face_descriptor, double>(
        umesh, umesh2, f2fmap_, "f:scalar"
      );
    }

    umesh.remove_property_map(fwhich);

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


  // ----------------------------------------------------------------------- //
  Rcpp::IntegerVector whereIs(Rcpp::NumericMatrix points) {

    if(!CGAL::is_triangle_mesh(mesh)) {
      Rcpp::stop("The mesh is not triangle.");
    }
    if(!CGAL::is_closed(mesh)) {
      Rcpp::stop("The mesh is not closed.");
    }

    CGAL::Side_of_triangle_mesh<EMesh3, EK> Where(mesh);

    int npoints = points.ncol();
    Rcpp::IntegerVector results(npoints);

    for(int i = 0; i < npoints; i++) {
      Rcpp::NumericVector point = points(Rcpp::_, i);
      EPoint3 pt(point(0), point(1), point(2));
      CGAL::Bounded_side side = Where(pt);
      int result = -1;
      if(side == CGAL::ON_BOUNDED_SIDE) {
        result = 1;
      } else if(side == CGAL::ON_BOUNDARY) {
        result = 0;
      }
      results(i) = result;
    }

    return results;
  }


  // ----------------------------------------------------------------------- //
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
  
};