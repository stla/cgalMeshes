#ifndef _HEADER_
#include "cgalMesh.h"
#endif

class CGALmesh {
public:
  EMesh3 mesh;
  MyMesh mymesh;
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
        MYMESH((
          xxx.mesh = mesh,
          xxx.normals = normals_,
          xxx.vcolors = vcolors_,
          xxx.fcolors = fcolors_
        ))
    ), 
      xptr(&mesh, false),
      normals(normals_),
      vcolors(vcolors_),
      fcolors(fcolors_) {}
  CGALmesh(Rcpp::XPtr<MyMesh> xptr_)
    : mymesh(*(xptr_.get())), 
      mesh(mymesh.mesh),
      xptr(Rcpp::XPtr<EMesh3>(&(mymesh.mesh), false)),
      normals(mymesh.normals),
      vcolors(mymesh.vcolors),
      fcolors(mymesh.fcolors) {}
  CGALmesh(const std::string filename, const bool binary)
    : mesh(readMeshFile(filename)), 
      mymesh(MYMESH((xxx.mesh = mesh))),
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
  
  Rcpp::XPtr<MyMesh> clone() {
    EMesh3 copy;
    CGAL::copy_face_graph(mesh, copy);
    MyMesh mycopy = MYMESH((
      xxx.mesh = copy,
      xxx.normals = normals,
      xxx.vcolors = vcolors,
      xxx.fcolors = fcolors
    ));
    return Rcpp::XPtr<MyMesh>(new MyMesh(mycopy), false);
  }
  
  Rcpp::List connectedComponents(const bool triangulate) {
    std::vector<EMesh3> cc_meshes;
    PMP::split_connected_components(mesh, cc_meshes);
    const size_t ncc = cc_meshes.size();
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
      MyMesh mycc = MYMESH((
        xxx.mesh = *cc
      ));
      xptrs(i) = Rcpp::XPtr<MyMesh>(new MyMesh(mycc), false);
      i++;
    }
    return xptrs;
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
      MyMesh mycmesh = MYMESH((
        xxx.mesh = cmesh
      ));
      out(i) = Rcpp::XPtr<MyMesh>(new MyMesh(mycmesh), false);
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
  
  Rcpp::List getRmesh(const bool normals) {
    if(CGAL::is_triangle_mesh(mesh)) {
      Rcpp::List rmesh = RSurfEKMesh2(mesh, normals, 3);
      rmesh["fcolors"] = fcolors;
      return rmesh;
    }
    if(CGAL::is_quad_mesh(mesh)) {
      return RSurfEKMesh2(mesh, normals, 4);
    }
    return RSurfEKMesh(mesh, normals);
  }
  
  Rcpp::XPtr<MyMesh> intersection(Rcpp::XPtr<EMesh3> mesh2XPtr) {
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
    MyMesh myimesh = MYMESH((
      xxx.mesh = imesh
    ));
    return Rcpp::XPtr<MyMesh>(new MyMesh(myimesh), false);
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
  
  Rcpp::XPtr<MyMesh> subtract(Rcpp::XPtr<EMesh3> mesh2XPtr) {
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
    MyMesh myimesh = MYMESH((
      xxx.mesh = imesh
    ));
    return Rcpp::XPtr<MyMesh>(new MyMesh(myimesh), false);
  }
  
  void triangulate() {
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
    normals = R_NilValue;
    fcolors = R_NilValue;
  }

  Rcpp::XPtr<MyMesh> Union(Rcpp::XPtr<EMesh3> mesh2XPtr) {
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
    MyMesh myimesh = MYMESH((
      xxx.mesh = imesh
    ));
    return Rcpp::XPtr<MyMesh>(new MyMesh(myimesh), false);
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