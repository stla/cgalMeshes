f <- Vectorize(function(u, v) {
  c(0, u, v) # or rbind(0, u, v)
})


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
#'
#' @examples
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

