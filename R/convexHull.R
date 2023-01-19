#' @title Convex hull
#' @description Mesh of the convex hull of a set of points.
#'
#' @param points numeric matrix with three columns, one point per row
#'
#' @return A \code{cgalMesh} object. The mesh is triangle.
#' @export
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' cube <- rbind(
#'   c(-1, -1, -1),
#'   c(-1, -1,  1),
#'   c(-1,  1, -1),
#'   c(-1,  1,  1),
#'   c( 1, -1, -1),
#'   c( 1, -1,  1),
#'   c( 1,  1, -1),
#'   c( 1,  1,  1),
#'   c( 0,  0,  0)
#' )
#' mesh <- convexHull(cube)
#' # plot
#' rmesh <- mesh$getMesh()
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(20, 20)
#' shade3d(rmesh, color = "orange")
#' wire3d(rmesh)
convexHull <- function(points) {
  stopifnot(is.matrix(points))
  stopifnot(ncol(points) == 3L)
  storage.mode(points) <- "double"
  if(anyNA(points)) {
    stop("Found missing values in the `points` matrix.")
  }
  xptr <- cxhull(t(points))
  cgalMesh$new(clean = xptr)
}

#' @title Intersection of convex hulls.
#' @description Mesh of the intersection of some convex hulls.
#'
#' @param Points a list of numeric matrices with three columns, each matrix 
#'   represents the points for which a convex hull will be computed
#' @param origin either \code{NULL} or a numeric vector of length 3 
#'   corresponding to a point in the intersection; set \code{NULL} if 
#'   you don't know such a point, otherwise give such a point to 
#'   gain speed
#'
#' @return A \code{cgalMesh} object.
#' @export
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' # intersection of the compound of five tetrahedra
#' mesh <- convexHullsIntersection(
#'   tetrahedraCompound[["Vertices"]],
#'   origin = c(0, 0, 0)
#' )
#' mesh$triangulate()
#' # plot
#' rmesh <- mesh$getMesh()
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(20, 20)
#' shade3d(rmesh, color = "orange")
#' wire3d(rmesh)
convexHullsIntersection <- function(Points, origin = NULL) {
  stopifnot(is.null(origin) || isVector3(origin))
  stopifnot(is.list(Points))
  stopifnot(length(Points) >= 2L)
  for(i in seq_along(Points)) {
    points <- Points[[i]]
    stopifnot(is.matrix(points))
    stopifnot(ncol(points) == 3L)
    storage.mode(points) <- "double"
    if(anyNA(points)) {
      stop("Found missing values.")
    }
    Points[[i]] <- t(points)
  }
  xptr <- cxhullsIntersection(Points, origin)
  cgalMesh$new(clean = xptr)
}
