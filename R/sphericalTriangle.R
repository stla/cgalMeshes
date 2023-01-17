#' @title Spherical triangle
#' @description Mesh of a spherical triangle.
#'
#' @param A,B,C the three vertices of the triangle, each given as a pair of 
#'   numbers: the polar angle (between 0 and pi) and the azimuthal angle 
#'   (between 0 and 2pi) 
#' @param center center of the sphere
#' @param radius radius of the sphere
#' @param iterations number of iterations used to construct the mesh of 
#'   the sphere
#'
#' @return A \code{cgalMesh} object. The mesh has normals.
#' @export
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' # orthant
#' A <- c(0, pi/2)
#' B <- c(pi/2, pi/2)
#' C <- c(0, 0)
#' mesh <- sphericalTriangle(A, B, C)
#' rmesh <- mesh$getMesh()
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(30, 30)
#' shade3d(rmesh, color = "red")
#' wire3d(rmesh)
sphericalTriangle <- function(
  A, B, C, center = c(0, 0, 0), radius = 1, iterations = 4    
) {
  A <- sph2cart(radius, A[1L], A[2L])
  B <- sph2cart(radius, B[1L], B[2L])
  C <- sph2cart(radius, C[1L], C[2L])
  stopifnot(isVector3(center))
  stopifnot(isPositiveNumber(radius))
  stopifnot(isStrictPositiveInteger(iterations))
  xptr <- sTriangle(A, B, C, center, radius, as.integer(iterations))
  cgalMesh$new(clean = xptr)
}