#' @title Hyperbolic triangle
#' @description Mesh of a hyperbolic triangle.
#'
#' @param A,B,C the vertices of the triangle
#' @param s hyperbolic curvature, a positive number
#' @param iterations number of iterations used to construct the mesh, 
#'   at least 1
#'
#' @return A \code{cgalMesh} object.
#' @export
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' mesh <- gyroTriangle(
#'   c(0, 0, 1), c(1, 0, 0), c(0, 1, 0), s = 0.7
#' )
#' mesh$computeNormals()
#' rmesh <- mesh$getMesh()
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(20, 20)
#' shade3d(rmesh, color = "green")
#' wire3d(rmesh)
gyroTriangle <- function(A, B, C, s, iterations = 3) {
  stopifnot(isVector3(A))
  stopifnot(isVector3(B))
  stopifnot(isVector3(C))
  stopifnot(isPositiveNumber(s))
  stopifnot(isStrictPositiveInteger(iterations))
  xptr <- gTriangle(
    A, B, C, s, as.integer(iterations) + 1L
  )
  cgalMesh$new(clean = xptr)
}