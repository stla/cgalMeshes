#ifndef _HEADER_
#include "cgalMesh.h"
#endif

// [[Rcpp::export]]
Rcpp::XPtr<EMesh3> AFSreconstruction_cpp(
  const Rcpp::NumericMatrix pts, const unsigned nneighs
) {

  std::vector<Point3> points = matrix_to_points3<Point3>(pts);

  if(nneighs >= 2) {
    CGAL::jet_smooth_point_set<CGAL::Sequential_tag>(points, nneighs);
  }
  
  AFS_triangulation3 dt(points.begin(), points.end());
  AFS_reconstruction reconstruction(dt);
  reconstruction.run();
  const AFS_Tds2& tds = reconstruction.triangulation_data_structure_2();
  
  std::vector<EPoint3> vertices;
  vertices.reserve(pts.ncol());
  int counter = 0;
  for(
    AFS_Tds2::Face_iterator fit = tds.faces_begin(); 
    fit != tds.faces_end();++fit
  ) {
    if(reconstruction.has_on_surface(fit)) {
      counter++;
      AFS_triangulation3::Facet f = fit->facet();
      AFS_triangulation3::Cell_handle ch = f.first;
      int ci = f.second;
      Point3 vs[3];
      for(int i = 0, j = 0; i < 4; i++) {
        if(ci != i) {
          vs[j] = ch->vertex(i)->point();
          j++;
        }
      }
      for(size_t k = 0; k < 3; k++) {
        const Point3 p = vs[k];
        const EPoint3 v = EPoint3(p.x(), p.y(), p.z());
        vertices.push_back(v);
      }
    }
  }
  std::vector<std::vector<int>> triangles;
  triangles.reserve(counter);
  for(int i = 0; i < counter; i++) {
    const int k = 3*i;
    const std::vector<int> triangle = {k, k+1, k+2};
    triangles.emplace_back(triangle);
  }
  
  size_t x = PMP::merge_duplicate_points_in_polygon_soup(vertices, triangles);
  EMesh3 mesh = csoup2mesh<EMesh3, EPoint3>(vertices, triangles, false);
  return Rcpp::XPtr<EMesh3>(new EMesh3(mesh), false);
}
