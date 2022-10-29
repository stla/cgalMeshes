#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include "MODULE.h"

RCPP_MODULE(class_CGALmesh) {
  using namespace Rcpp;
  class_<CGALmesh>("CGALmesh")
    .constructor<const NumericMatrix, const List, const bool>()
    .constructor<XPtr<EMesh3>>()
    .constructor<Rcpp::String, const bool>()
    .field("xptr", &CGALmesh::xptr)
    .method("centroid", &CGALmesh::centroid)
    .method("clone", &CGALmesh::clone)
    .method("doesBoundVolume", &CGALmesh::doesBoundVolume)
    .method("doesSelfIntersect", &CGALmesh::doesSelfIntersect)
    .method("edges", &CGALmesh::edges)
    .method("getRmesh", &CGALmesh::getRmesh)
    .method("isClosed", &CGALmesh::isClosed)
    .method("isTriangle", &CGALmesh::isTriangle)
    .method("print", &CGALmesh::print)
    .method("reverseFaceOrientations", &CGALmesh::reverseFaceOrientations)
    .method("triangulate", &CGALmesh::triangulate)
    .method("vertices", &CGALmesh::vertices)
    .method("writeFile", &CGALmesh::writeFile);
}