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
