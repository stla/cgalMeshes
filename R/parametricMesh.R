#' @title Mesh of a parametric surface
#' @description Mesh of a parametric surface.
#'
#' @param f vectorized function of two variables returning a three-dimensional 
#'   point 
#' @param urange,vrange the ranges of the two variables
#' @param periodic two Boolean values, whether \code{f} is periodic in the 
#'   first variable and in the second variable
#' @param nu,nv numbers of subdivisions
#'
#' @return A \strong{rgl} mesh (\code{mesh3d} object).
#' @export
#' @importFrom rgl tmesh3d
#' @examples 
#' library(cgalMeshes)
#' library(rgl)
#' # Richmond surface
#' Richmond <- function(r, t){
#'   radius <- 0.5
#'   r <- r + 1/4
#'   exprho <- exp(r*(1 + 3*radius) - 2 - radius)
#'   u <- exprho * cospi(2*t)
#'   v <- exprho * sinpi(2*t)
#'   rbind(
#'     -v/(u*u+v*v) - u*u*v + v*v*v/3,
#'     3*u,
#'      u/(u*u+v*v) - u*v*v + u*u*u/3
#'   ) 
#' }
#' rmesh <- parametricMesh(
#'   Richmond, urange = c(0, 1), vrange = c(0, 1),
#'   periodic = c(FALSE, TRUE), nu = 100, nv = 100
#' )
#' rmesh <- addNormals(rmesh)
#' # plot
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(50, 10, zoom = 0.8)
#' shade3d(rmesh, color = "gold")
parametricMesh <- function(
    f, urange, vrange, periodic = c(FALSE, FALSE), nu, nv
) {
  stopifnot(isPositiveInteger(nu), nu >= 3)
  stopifnot(isPositiveInteger(nv), nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  uperiodic <- periodic[1L]
  vperiodic <- periodic[2L]
  if(uperiodic) {
    u_ <- seq(urange[1L], urange[2L], length.out = nu+1L)[-1L]
  } else {
    u_ <- seq(urange[1L], urange[2L], length.out = nu)
  }
  if(vperiodic) {
    v_ <- seq(vrange[1L], vrange[2L], length.out = nv+1L)[-1L]
  } else {
    v_ <- seq(vrange[1L], vrange[2L], length.out = nv)
  }
  Grid <- expand.grid(U = u_, V = v_)
  varray <- with(Grid, array(f(U, V), dim = c(3L, nu, nv)))
  varray2 <- aperm(varray, c(1L, 3L, 2L))
  vs <- matrix(varray2, nrow = 3L, ncol = nu*nv)
  tris <- meshTopology(nu, nv, uperiodic, vperiodic)
  tmesh3d(
    vertices    = vs,
    indices     = tris,
    homogeneous = FALSE
  )
}

