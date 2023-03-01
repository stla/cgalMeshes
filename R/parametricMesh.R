#' @title Mesh of a parametric surface
#' @description Mesh of a parametric surface.
#'
#' @param f vectorized function of two variables returning a three-dimensional 
#'   point 
#' @param urange,vrange the ranges of the two variables
#' @param periodic two Boolean values, whether \code{f} is periodic in the 
#'   first variable and in the second variable
#' @param nu,nv numbers of subdivisions
#' @param fnormal if given (i.e. not \code{NULL}), a vectorized function 
#'   returning at \code{(u,v)} the normal for the vertex \code{f(u,v)}
#'
#' @return A \strong{rgl} mesh (\code{mesh3d} object), with normals if 
#'   \code{fnormal} is given.
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
#' 
#' # an example with normals ####
#' library(cgalMeshes)
#' library(rgl)
#' # twisted sphere
#' fx <- function(u, v) cos(u) * cos(v)
#' fy <- function(u, v) sin(u) * cos(v)
#' fz <- function(u, v) sin(v) + u
#' f <- function(u, v) {
#'   rbind(fx(u, v), fy(u, v), fz(u, v))
#' }
#' if(require("dfdr")) {
#'   # compute the gradient with the 'dfdr' package
#'   dfx <- gradient(fx, FALSE, u, v)
#'   dfy <- gradient(fy, FALSE, u, v)
#'   dfz <- gradient(fz, FALSE, u, v)
#'   fnormal <- Vectorize(function(u, v) {
#'     dx <- dfx(u, v)
#'     dy <- dfy(u, v)
#'     dz <- dfz(u, v)
#'     v1 <- c(dx[1L], dy[1L], dz[1L])
#'     v2 <- c(dx[2L], dy[2L], dz[2L])
#'     - crossProduct(v1, v2)
#'   })
#' } else {
#'   fnormal <- NULL
#' }
#' # compute mesh of the parametric surface
#' rmesh <- cgalMeshes::parametricMesh(
#'   f, c(0, 2*pi), c(0, 2*pi), c(FALSE, TRUE), nu = 60L, nv = 30L, fnormal
#' )
#' # plot
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(50, 10)
#' shade3d(rmesh, color = "midnightblue", specular = "black")
parametricMesh <- function(
    f, urange, vrange, periodic = c(FALSE, FALSE), nu, nv, fnormal = NULL
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
  if(!is.null(fnormal)) {
    narray <- with(Grid, array(fnormal(U, V), dim = c(3L, nu, nv)))
    narray2 <- aperm(narray, c(1L, 3L, 2L))
    normals <- t(matrix(narray2, nrow = 3L, ncol = nu*nv))
  } else {
    normals <- NULL
  }
  tmesh3d(
    vertices    = vs,
    indices     = tris,
    normals     = normals,
    homogeneous = FALSE
  )
}
