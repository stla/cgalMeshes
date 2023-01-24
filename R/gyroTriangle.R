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
#'
#' # hyperbolic icosahedron ####
#' library(cgalMeshes)
#' library(rgl)
#' icosahedron <- icosahedron3d()
#' vertices    <- icosahedron[["vb"]][-4L, ]
#' faces       <- icosahedron[["it"]]
#' colors <- hcl.colors(ncol(faces), palette = "Plasma")
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(15, -60, zoom = 0.7)
#' for(i in 1L:ncol(faces)) {
#'   triangle <- faces[, i]
#'   A <- vertices[, triangle[1L]]
#'   B <- vertices[, triangle[2L]]
#'   C <- vertices[, triangle[3L]]
#'   mesh <- gyroTriangle(A, B, C, s = 0.6)
#'   mesh$computeNormals()
#'   rmesh <- mesh$getMesh()
#'   shade3d(rmesh, color = colors[i])
#'   wire3d(rmesh)
#' }
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