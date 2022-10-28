#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include "MODULE.h"

RCPP_MODULE(class_CGALmesh) {
  using namespace Rcpp;
  class_<CGALmesh>("CGALmesh")
    .constructor<const NumericMatrix, const List, const bool>()
    .constructor<XPtr<EMesh3>>()
    .field("xptr", &CGALmesh::xptr)
    .method("centroid", &CGALmesh::centroid)
    .method("getRmesh", &CGALmesh::getRmesh)
    .method("isTriangle", &CGALmesh::isTriangle);
}