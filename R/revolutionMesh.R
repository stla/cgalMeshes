#' @title Surface of revolution
#' @description Mesh of a surface of revolution. The axis of revolution is 
#'   the z-axis.
#'
#' @param x,y two numeric vectors of the same length defining the section to 
#'   be revoluted
#' @param n integer, the number of subdivisions used to construct the mesh
#'
#' @return A \strong{rgl} triangle mesh (class \code{mesh3d}).
#' @export
#' @importFrom rgl tmesh3d
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' t <- seq(0, 2*pi, length.out = 90)
#' x <- 4 + cos(t)/2
#' y <- sin(t)
#' rmesh <- revolutionMesh(x, y, n = 120)
#' rmesh <- addNormals(rmesh)
#' shade3d(rmesh, color = "red")
revolutionMesh <- function(x, y, n = 100) {
  stopifnot(isPositiveInteger(n), n >= 3)
  stopifnot(is.numeric(x), is.numeric(y), length(x) == length(y))
  if(anyNA(x) || anyNA(y)) {
    stop("Found missing values.")
  }
  nu <- as.integer(n)
  uperiodic <- TRUE
  u_ <- seq(0, 2, length.out = nu+1L)[-1L]
  l <- length(x)
  stopifnot(l >= 3L)
  vperiodic <- isTRUE(all.equal(c(x[1L], y[1L]), c(x[l], y[l])))
  nv <- if(vperiodic) { l - 1L } else { l }
  v_ <- 1L:nv
  f <- Vectorize(function(u, j) {
    c(x[j] * cospi(u), x[j] * sinpi(u), y[j])
  })
  Grid <- expand.grid(U = u_, V = v_)
  varray  <- with(Grid, array(f(U, V), dim = c(3L, nu, nv)))
  varray2 <- aperm(varray, c(1L, 3L, 2L))
  vs <- matrix(varray2, nrow = 3L, ncol = nu*nv)
  tris <- meshTopology(nu, nv, uperiodic, vperiodic)
  tmesh3d(
    vertices    = vs,
    indices     = tris,
    homogeneous = FALSE
  )
}
