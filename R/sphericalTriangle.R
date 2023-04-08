#' @title Spherical triangle
#' @description Mesh of a spherical triangle.
#'
#' @param A,B,C the three vertices of the triangle, each given either as a pair 
#'   of numbers: the polar angle (between 0 and pi) and the azimuthal angle 
#'   (between 0 and 2pi), or as Cartesian coordinates 
#' @param center center of the sphere
#' @param radius radius of the sphere, ignored if the three vertices 
#'   are given as Cartesian coordinates
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
#'
#' # spherical icosahedron ####
#' library(cgalMeshes)
#' library(rgl)
#' icosahedron <- icosahedron3d()
#' vertices    <- icosahedron[["vb"]][-4L, ]
#' faces       <- icosahedron[["it"]]
#' colors <- rainbow(ncol(faces))
#' \donttest{open3d(windowRect = 50 + c(0, 0, 512, 512))
#' for(i in 1L:ncol(faces)) {
#'   triangle <- faces[, i]
#'   A <- vertices[, triangle[1L]]
#'   B <- vertices[, triangle[2L]]
#'   C <- vertices[, triangle[3L]]
#'   mesh <- sphericalTriangle(A, B, C)
#'   rmesh <- mesh$getMesh()
#'   shade3d(rmesh, color = colors[i])
#'   wire3d(rmesh)
#' }}
sphericalTriangle <- function(
  A, B, C, center = c(0, 0, 0), radius = 1, iterations = 4    
) {
  stopifnot(isVector3(center))
  if(length(A) == 2L) {
    A <- sph2cart(radius, A[1L], A[2L])
    B <- sph2cart(radius, B[1L], B[2L])
    C <- sph2cart(radius, C[1L], C[2L])
  } else {
    stopifnot(isVector3(A))
    stopifnot(isVector3(B))
    stopifnot(isVector3(C))
    r1 <- c(crossprod(A - center))
    r2 <- c(crossprod(B - center))
    r3 <- c(crossprod(C - center))
    if(!isTRUE(all.equal(c(r2-r1, r3-r1), c(0, 0)))) {
      stop("The three points are not on the sphere centered at `center`.")
    }
    radius <- sqrt(r1)
  }
  if(CGALversion() < 5.5) {
    x <- sqrt(1 + (1 + sqrt(5)) / 4) # bug make_icosahedron
    radius <- radius / x
  }
  stopifnot(isPositiveNumber(radius))
  stopifnot(isStrictPositiveInteger(iterations))
  xptr <- sTriangle(A, B, C, center, radius, as.integer(iterations))
  cgalMesh$new(clean = xptr)
}