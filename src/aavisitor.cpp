#ifndef _HEADER_
#include "cgalMesh.h"
#endif

//template <class MeshT>
//template <typename>
// void Visitor::new_vertex_added(std::size_t i_id, vertex_descriptor v, const EMesh3 & tm) {
//   Rcpp::Rcout << v << "\n";
// }

void new_vertex_added(std::size_t i_id, vertex_descriptor v, const EMesh3 & tm) {
  Rcpp::Rcout << v << "\n";
}
//template void Visitor::new_vertex_added(std::size_t, vertex_descriptor, const EMesh3&);

//template void Visitor<EMesh3>::new_vertex_added(std::size_t, vertex_descriptor, const EMesh3&);
// void Visitor<EMesh3>::new_vertex_added(std::size_t i_id, vertex_descriptor v, const EMesh3 & tm){
//   Rcpp::Rcout << v << "\n";  
// }
