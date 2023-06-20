#' @title Mesh of a parametric surface
#' @description Mesh of a parametric surface.
#'
#' @param f vectorized function of two variables returning a three-dimensional 
#'   point 
#' @param urange,vrange the ranges of the two variables
#' @param periodic two Boolean values, whether \code{f} is periodic in the 
#'   first variable and in the second variable
#' @param nu,nv numbers of subdivisions
#' @param clean Boolean, whether to merge the duplicated vertices in the 
#'   output mesh; use only if necessary because it costs time
#' @param fcolor if given (i.e. not \code{NULL}), a vectorized function 
#'   returning at \code{(u,v)} a color for the vertex \code{f(u,v)}; 
#'   the behavior is undefined if some vertices are merged by \code{clean=TRUE}
#' @param fnormal if given (i.e. not \code{NULL}), a vectorized function 
#'   returning at \code{(u,v)} the normal for the vertex \code{f(u,v)}; 
#'   the behavior is undefined if some vertices are merged by \code{clean=TRUE}
#'
#' @return A \strong{rgl} mesh (\code{mesh3d} object), with normals if 
#'   \code{fnormal} is given.
#' @export
#' @importFrom rgl tmesh3d
#' @importFrom data.table uniqueN
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
#'   f, c(0, 2*pi), c(0, 2*pi), c(FALSE, TRUE), nu = 60L, nv = 30L, 
#'   fnormal = fnormal
#' )
#' # plot
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(50, 10) 
#' shade3d(rmesh, color = "midnightblue", specular = "black")
#' 
#' ##| example with colors: torus figure 8 ####
#' library(rgl)
#' library(cgalMeshes)
#' library(colorspace) # to use the  lighten  function
#' # function mapping (u,v) to a color
#' clrs <- c(
#'   hcl.colors(128L, "Rocket"), 
#'   hcl.colors(128L, "Rocket", rev = TRUE)
#' )
#' framp <- colorRamp(clrs)
#' fcolor <- function(u, v) {
#'   cols <- framp(v/(2*pi))
#'   cols <- rgb(cols[, 1L], cols[, 2L], cols[, 3L], maxColorValue = 255)
#'   lighten(cols, 0.3*cos(u))
#' }
#' # parameterization
#' f <- function(u, v, c = 1){
#'   h <- c + sin(v) * cos(u) - sin(2*v) * sin(u) / 2
#'   x <- h * cos(u)
#'   y <- h * sin(u)
#'   z <- sin(u) * sin(v) + cos(u) * sin(2*v) / 2
#'   rbind(x, y, z)  
#' }
#' # make the mesh and plot it
#' \donttest{rmesh <- parametricMesh(
#'   f, c(0, 2*pi), c(0, 2*pi), periodic = c(TRUE, TRUE),
#'   nu = 100L, nv = 100L, fcolor = fcolor)
#' rmesh <- addNormals(rmesh)
#' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
#' shade3d(rmesh, meshColor = "vertices")}
parametricMesh <- function(
    f, urange, vrange, periodic = c(FALSE, FALSE), 
    nu, nv, clean = FALSE, fcolor = NULL, fnormal = NULL
) {
  stopifnot(isPositiveInteger(nu), nu >= 3)
  stopifnot(isPositiveInteger(nv), nv >= 3)
  stopifnot(isBoolean(clean))
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
  if(!is.null(fcolor)) {
    colors <- with(Grid, fcolor(U, V))
  }
  if(!clean && !is.null(fnormal)) {
    narray <- with(Grid, array(fnormal(U, V), dim = c(3L, nu, nv)))
    narray2 <- aperm(narray, c(1L, 3L, 2L))
    normals <- t(matrix(narray2, nrow = 3L, ncol = nu*nv))
  } else {
    normals <- NULL
  }
  if(clean) {
    gather <- gatherVertices(vs, tris)
    validFaces <- apply(gather[["faces"]], 2L, uniqueN) == 3L
    tmesh3d(
      vertices    = gather[["vertices"]],
      indices     = gather[["faces"]][, validFaces],
      normals     = NULL,
      material    = if(!is.null(fcolor)) list(color = colors),
      homogeneous = FALSE
    )
  } else {
    tmesh3d(
      vertices    = vs,
      indices     = tris,
      normals     = normals,
      material    = if(!is.null(fcolor)) list(color = colors),
      homogeneous = FALSE
    )
  }
}
