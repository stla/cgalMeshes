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
    .constructor<Rcpp::String, bool, bool>()
    .field("xptr", &CGALmesh::xptr)
    .method("area", &CGALmesh::area)
    .method("assignFaceColors", &CGALmesh::assignFaceColors)
    .method("assignNormals", &CGALmesh::assignNormals)
    .method("assignVertexColors", &CGALmesh::assignVertexColors)
    .method("assignFaceScalars", &CGALmesh::assignFaceScalars)
    .method("assignVertexScalars", &CGALmesh::assignVertexScalars)
    .method("boundingBox", &CGALmesh::boundingBox)
    .method("CatmullClark", &CGALmesh::CatmullClark)
    .method("centroid", &CGALmesh::centroid)
    .method("clipMesh", &CGALmesh::clipMesh)
    .method("clipToIsoCuboid", &CGALmesh::clipToIsoCuboid)
    .method("clipToPlane", &CGALmesh::clipToPlane)
    .method("clone", &CGALmesh::clone)
    .method("collectGarbage", &CGALmesh::collectGarbage)
    .method("computeNormals", &CGALmesh::computeNormals)
    .method("connectedComponents", &CGALmesh::connectedComponents)
    .method("convexParts", &CGALmesh::convexParts)
    .method("distance", &CGALmesh::distance)
    .method("doesBoundVolume", &CGALmesh::doesBoundVolume)
    .method("doesSelfIntersect", &CGALmesh::doesSelfIntersect)
    .method("DooSabin", &CGALmesh::DooSabin)
    .method("dual", &CGALmesh::dual)
    .method("edges", &CGALmesh::edges)
    .method("facesAroundVertex", &CGALmesh::facesAroundVertex)
    .method("fair", &CGALmesh::fair)
    .method("fillBoundaryHole", &CGALmesh::fillBoundaryHole)
    .method("filterMesh", &CGALmesh::filterMesh)
    .method("fixManifoldness", &CGALmesh::fixManifoldness)
    .method("geoDists", &CGALmesh::geoDists)
    .method("getBorders", &CGALmesh::getBorders)
    .method("getFacesInfo", &CGALmesh::getFacesInfo)
    .method("getFacesList", &CGALmesh::getFacesList)
    .method("getFacesMatrix", &CGALmesh::getFacesMatrix)
    .method("getFcolors", &CGALmesh::getFcolors)
    .method("getVcolors", &CGALmesh::getVcolors)
    .method("getFscalars", &CGALmesh::getFscalars)
    .method("getVscalars", &CGALmesh::getVscalars)
    .method("getVertices", &CGALmesh::getVertices)
    .method("getNormals", &CGALmesh::getNormals)
    .method("getRmesh", &CGALmesh::getRmesh)
    .method("HausdorffApproximate", &CGALmesh::HausdorffApproximate)
    .method("HausdorffEstimate", &CGALmesh::HausdorffEstimate)
    .method("intersection", &CGALmesh::intersection)
    .method("isClosed", &CGALmesh::isClosed)
    .method("isotropicRemeshing", &CGALmesh::isotropicRemeshing)
    .method("isOutwardOriented", &CGALmesh::isOutwardOriented)
    .method("isQuad", &CGALmesh::isQuad)
    .method("isTriangle", &CGALmesh::isTriangle)
    .method("isValid", &CGALmesh::isValid)
    .method("isValidFaceGraph", &CGALmesh::isValidFaceGraph)
    .method("isValidHalfedgeGraph", &CGALmesh::isValidHalfedgeGraph)
    .method("isValidPolygonMesh", &CGALmesh::isValidPolygonMesh)
    .method("LoopSubdivision", &CGALmesh::LoopSubdivision)
    .method("merge", &CGALmesh::merge)
    .method("optimalBoundingBox", &CGALmesh::optimalBoundingBox)
    .method("orientToBoundVolume", &CGALmesh::orientToBoundVolume)
    .method("print", &CGALmesh::print)
    .method("removeSelfIntersections", &CGALmesh::removeSelfIntersections)
    .method("reverseFaceOrientations", &CGALmesh::reverseFaceOrientations)
    .method("sampleInMesh", &CGALmesh::sampleInMesh)
    .method("sampleOnMesh", &CGALmesh::sampleOnMesh)
    .method("sharpEdges", &CGALmesh::sharpEdges)
    .method("Sqrt3Subdivision", &CGALmesh::Sqrt3Subdivision)
    .method("subtract", &CGALmesh::subtract)
    .method("triangulate", &CGALmesh::triangulate)
    .method("Union", &CGALmesh::Union)
    .method("volume", &CGALmesh::volume)
    .method("whereIs", &CGALmesh::whereIs)
    .method("writeFile", &CGALmesh::writeFile);
}
