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
    .method("area", &CGALmesh::area)
    .method("centroid", &CGALmesh::centroid)
    .method("clipMesh", &CGALmesh::clipMesh)
    .method("clone", &CGALmesh::clone)
    .method("connectedComponents", &CGALmesh::connectedComponents)
    .method("convexParts", &CGALmesh::convexParts)
    .method("distance", &CGALmesh::distance)
    .method("doesBoundVolume", &CGALmesh::doesBoundVolume)
    .method("doesSelfIntersect", &CGALmesh::doesSelfIntersect)
    .method("dual", &CGALmesh::dual)
    .method("edges", &CGALmesh::edges)
    .method("fair", &CGALmesh::fair)
    .method("geoDists", &CGALmesh::geoDists)
    .method("getRmesh", &CGALmesh::getRmesh)
    .method("intersection", &CGALmesh::intersection)
    .method("isClosed", &CGALmesh::isClosed)
    .method("isOutwardOriented", &CGALmesh::isOutwardOriented)
    .method("isTriangle", &CGALmesh::isTriangle)
    .method("isValid", &CGALmesh::isValid)
    .method("orientToBoundVolume", &CGALmesh::orientToBoundVolume)
    .method("print", &CGALmesh::print)
    .method("removeSelfIntersections", &CGALmesh::removeSelfIntersections)
    .method("reverseFaceOrientations", &CGALmesh::reverseFaceOrientations)
    .method("subtract", &CGALmesh::subtract)
    .method("triangulate", &CGALmesh::triangulate)
    .method("Union", &CGALmesh::Union)
    .method("vertices", &CGALmesh::vertices)
    .method("volume", &CGALmesh::volume)
    .method("writeFile", &CGALmesh::writeFile);
}
