#ifndef _HEADER_
#include "cgalMesh.h"
#endif

void Message(std::string msg) {
  SEXP rmsg = Rcpp::wrap(msg);
  Rcpp::message(rmsg);
}

std::string q2str(CGAL::Gmpq r) {
  CGAL::Gmpz numer = r.numerator();
  CGAL::Gmpz denom = r.denominator();
  size_t n = mpz_sizeinbase(numer.mpz(), 10) + 2;
  size_t d = mpz_sizeinbase(denom.mpz(), 10) + 2;
  char* cnumer = new char[n];
  char* cdenom = new char[d];
  cnumer = mpz_get_str(cnumer, 10, numer.mpz());
  cdenom = mpz_get_str(cdenom, 10, denom.mpz());
  std::string snumer = cnumer;
  std::string sdenom = cdenom;
  delete[] cnumer;
  delete[] cdenom;
  return snumer + "/" + sdenom;
}

template <typename PointT>
std::vector<PointT> matrix_to_points3(const Rcpp::NumericMatrix M) {
  const size_t npoints = M.ncol();
  std::vector<PointT> points;
  points.reserve(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::NumericVector pt = M(Rcpp::_, i);
    points.emplace_back(PointT(pt(0), pt(1), pt(2)));
  }
  return points;
}

template std::vector<EPoint3> matrix_to_points3<EPoint3>(
    const Rcpp::NumericMatrix);

std::vector<QPoint3> matrix_to_qpoints3(const Rcpp::CharacterMatrix M) {
  const size_t npoints = M.ncol();
  std::vector<QPoint3> points;
  points.reserve(npoints);
  for(size_t i = 0; i < npoints; i++) {
    const Rcpp::CharacterVector pt = M(Rcpp::_, i);
    CGAL::Gmpq qpt0(CGAL::Gmpq(Rcpp::as<std::string>(pt(0))));
    CGAL::Gmpq qpt1(CGAL::Gmpq(Rcpp::as<std::string>(pt(1))));
    CGAL::Gmpq qpt2(CGAL::Gmpq(Rcpp::as<std::string>(pt(2))));
    points.emplace_back(QPoint3(qpt0, qpt1, qpt2));
  }
  return points;
}

std::vector<std::vector<int>> matrix_to_Tfaces(
    const Rcpp::IntegerMatrix Faces) {
  const size_t nfaces = Faces.ncol();
  std::vector<std::vector<int>> faces;
  faces.reserve(nfaces);
  for(size_t i = 0; i < nfaces; i++) {
    const Rcpp::IntegerVector face_rcpp = Faces(Rcpp::_, i);
    std::vector<int> face = {face_rcpp(0), face_rcpp(1), face_rcpp(2)};
    faces.emplace_back(face);
  }
  return faces;
}

std::vector<std::vector<int>> list_to_faces(const Rcpp::List L) {
  const size_t nfaces = L.size();
  std::vector<std::vector<int>> faces;
  faces.reserve(nfaces);
  for(size_t i = 0; i < nfaces; i++) {
    Rcpp::IntegerVector face_rcpp = Rcpp::as<Rcpp::IntegerVector>(L(i));
    std::vector<int> face(face_rcpp.begin(), face_rcpp.end());
    faces.emplace_back(face);
  }
  return faces;
}

template <typename MeshT, typename PointT>
MeshT soup2mesh(std::vector<PointT> points,
                std::vector<std::vector<int>> faces,
                const bool clean,
                const bool triangulate,
                const bool stopifnotclosed) {
  bool success = PMP::orient_polygon_soup(points, faces);
  if(success) {
    Message("Successful polygon orientation.");
  } else {
    Message("Polygon orientation failed.");
  }
  if(clean) {
    PMP::repair_polygon_soup(points, faces);
  }
  MeshT mesh;
  PMP::polygon_soup_to_polygon_mesh(points, faces, mesh);
  const bool valid = mesh.is_valid(false);
  if(!valid) {
    Message("The mesh is not valid.");
  }
  bool isTriangle;
  if(triangulate) {
    Message("Triangulation.");
    const bool success = PMP::triangulate_faces(mesh);
    if(!success) {
      Rcpp::stop("Triangulation has failed.");
    }
    isTriangle = true;
  } else {
    isTriangle = CGAL::is_triangle_mesh(mesh);
  }
  if(isTriangle) {
    Message("The mesh is triangle.");
  } else {
    Message(
        "The mesh is not triangle; no way to ensure it bounds a volume "
        "and whether it is outward oriented.");
  }
  if(CGAL::is_closed(mesh)) {
    Message("The mesh is closed.");
    if(isTriangle) {
      if(!PMP::is_outward_oriented(mesh)) {
        PMP::reverse_face_orientations(mesh);
      }
      const bool bv = PMP::does_bound_a_volume(mesh);
      std::string msg2;
      if(bv) {
        msg2 = "The mesh bounds a volume.";
      } else {
        msg2 = "The mesh does not bound a volume - reorienting.";
        PMP::orient_to_bound_a_volume(mesh);
      }
      Message(msg2);
    }
  } else {
    if(stopifnotclosed) {
      Rcpp::stop("The mesh is not closed.");
    } else {
      Message("The mesh is not closed.");
    }
  }
  return mesh;
}

template <typename MeshT, typename PointT>
MeshT makeSurfMesh(
  const Rcpp::List rmesh, const bool clean, const bool triangulate,
  const bool stopifnotclosed
) {
  const Rcpp::NumericMatrix vertices =
      Rcpp::as<Rcpp::NumericMatrix>(rmesh["vertices"]);
  const Rcpp::List rfaces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<PointT> points = matrix_to_points3<PointT>(vertices);
  std::vector<std::vector<int>> faces = list_to_faces(rfaces);
  return soup2mesh<MeshT, PointT>(
      points, faces, clean, triangulate, stopifnotclosed
  );
}

template EMesh3 makeSurfMesh<EMesh3, EPoint3>(
  const Rcpp::List, const bool, const bool, const bool
);

QMesh3 makeSurfQMesh(
  const Rcpp::List rmesh, const bool clean, const bool triangulate,
  const bool stopifnotclosed
) {
  const Rcpp::CharacterMatrix vertices =
      Rcpp::as<Rcpp::CharacterMatrix>(rmesh["vertices"]);
  const Rcpp::List rfaces = Rcpp::as<Rcpp::List>(rmesh["faces"]);
  std::vector<QPoint3> points = matrix_to_qpoints3(vertices);
  std::vector<std::vector<int>> faces = list_to_faces(rfaces);
  return soup2mesh<QMesh3, QPoint3>(
      points, faces, clean, triangulate, stopifnotclosed
  );
}
