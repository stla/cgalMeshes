#ifndef _HEADER_
#include "cgalMesh.h"
#endif

#include "MODULE.h"

RCPP_MODULE(class_CGALmesh) {
  using namespace Rcpp;
  class_<CGALmesh>("CGALmesh")
    .constructor<
      const NumericMatrix, const List, bool,
      Nullable<NumericMatrix>, Nullable<StringVector>, Nullable<StringVector>
    >()
    .constructor<XPtr<EMesh3>>()
    .constructor<Rcpp::String, const bool>()
    .field("xptr", &CGALmesh::xptr)
    .method("area", &CGALmesh::area)
    .method("assignFaceColors", &CGALmesh::assignFaceColors)
    .method("assignVertexColors", &CGALmesh::assignVertexColors)
    .method("assignFaceScalars", &CGALmesh::assignFaceScalars)
    .method("assignVertexScalars", &CGALmesh::assignVertexScalars)
    .method("centroid", &CGALmesh::centroid)
    .method("clipMesh", &CGALmesh::clipMesh)
    .method("doubleclip", &CGALmesh::doubleclip)
    .method("clone", &CGALmesh::clone)
    .method("computeNormals", &CGALmesh::computeNormals)
    .method("connectedComponents", &CGALmesh::connectedComponents)
    .method("convexParts", &CGALmesh::convexParts)
    .method("distance", &CGALmesh::distance)
    .method("doesBoundVolume", &CGALmesh::doesBoundVolume)
    .method("doesSelfIntersect", &CGALmesh::doesSelfIntersect)
    .method("dual", &CGALmesh::dual)
    .method("edges", &CGALmesh::edges)
    .method("facesAroundVertex", &CGALmesh::facesAroundVertex)
    .method("fair", &CGALmesh::fair)
    .method("filterGraph", &CGALmesh::filterGraph)
    .method("filterMesh", &CGALmesh::filterMesh)
    .method("fixManifoldness", &CGALmesh::fixManifoldness)
    .method("geoDists", &CGALmesh::geoDists)
    .method("getFacesList", &CGALmesh::getFacesList)
    .method("getFacesMatrix", &CGALmesh::getFacesMatrix)
    .method("getFcolors", &CGALmesh::getFcolors)
    .method("getVcolors", &CGALmesh::getVcolors)
    .method("getFscalars", &CGALmesh::getFscalars)
    .method("getHalfedges", &CGALmesh::getHalfedges)
    .method("getVscalars", &CGALmesh::getVscalars)
    .method("getVertices", &CGALmesh::getVertices)
    .method("getNormals", &CGALmesh::getNormals)
    .method("getRmesh", &CGALmesh::getRmesh)
    .method("intersection", &CGALmesh::intersection)
    .method("isClosed", &CGALmesh::isClosed)
    .method("isOutwardOriented", &CGALmesh::isOutwardOriented)
    .method("isQuad", &CGALmesh::isQuad)
    .method("isTriangle", &CGALmesh::isTriangle)
    .method("isValid", &CGALmesh::isValid)
    .method("isValidFaceGraph", &CGALmesh::isValidFaceGraph)
    .method("isValidHalfedgeGraph", &CGALmesh::isValidHalfedgeGraph)
    .method("isValidPolygonMesh", &CGALmesh::isValidPolygonMesh)
    .method("orientToBoundVolume", &CGALmesh::orientToBoundVolume)
    .method("print", &CGALmesh::print)
    .method("removeSelfIntersections", &CGALmesh::removeSelfIntersections)
    .method("reverseFaceOrientations", &CGALmesh::reverseFaceOrientations)
    .method("subtract", &CGALmesh::subtract)
    .method("triangulate", &CGALmesh::triangulate)
    .method("Union", &CGALmesh::Union)
    .method("volume", &CGALmesh::volume)
    .method("writeFile", &CGALmesh::writeFile)
    .method("testsplit", &CGALmesh::testsplit);
}
