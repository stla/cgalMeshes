#' @title Sphere mesh
#' @description Triangle mesh of a sphere.
#'
#' @param x,y,z coordinates of the center
#' @param r radius
#' @param iterations number of iterations (the mesh is obtained by iteratively 
#'   subdividing the faces of an icosahedron)
#'
#' @return A \strong{rgl} mesh (class \code{mesh3d}).
#' @export
#'
#' @importFrom rgl subdivision3d icosahedron3d
sphereMesh <- function(x = 0, y = 0, z = 0, r = 1, iterations = 3L) {
  stopifnot(isNumber(x), isNumber(y), isNumber(z))
  stopifnot(isPositiveNumber(r))
  stopifnot(isStrictPositiveInteger(iterations))
  sphere <- subdivision3d(icosahedron3d(), depth = iterations)
  vs <- sphere[["vb"]][1L:3L, ]
  h <- sqrt(apply(vs, 2L, function(x) sum(x * x)))
  vs[1L, ] <- vs[1L, ] / h
  vs[2L, ] <- vs[2L, ] / h
  vs[3L, ] <- vs[3L, ] / h
  sphere[["vb"]] <- rbind(r*vs + c(x, y, z), 1) 
  sphere[["normals"]] <- vs
  sphere
}

#' @title Torus mesh
#' @description Triangle mesh of a torus.
#'
#' @param R,r major and minor radii, positive numbers
#' @param nu,nv numbers of subdivisions, integers (at least 3)
#'
#' @return A triangle \strong{rgl} mesh (class \code{mesh3d}).
#' @export
#'
#' @importFrom rgl tmesh3d
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' mesh <- torusMesh(R = 3, r = 1)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "green")
#' wire3d(mesh)
torusMesh <- function(R, r, nu = 50, nv = 30){
  stopifnot(isPositiveNumber(R), isPositiveNumber(r))
  stopifnot(R > r)
  stopifnot(nu >= 3, nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  nunv <- nu*nv
  vs      <- matrix(NA_real_, nrow = 3L, ncol = nunv)
  normals <- matrix(NA_real_, nrow = nunv, ncol = 3L)
  tris1   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  tris2   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  u_ <- seq(0, 2*pi, length.out = nu + 1L)[-1L]
  cosu_ <- cos(u_)
  sinu_ <- sin(u_)
  v_ <- seq(0, 2*pi, length.out = nv + 1L)[-1L]
  cosv_ <- cos(v_)
  sinv_ <- sin(v_)
  Rrcosv_ <- R + r*cosv_
  rsinv_ <- r*sinv_
  jp1_ <- c(2L:nv, 1L)
  j_ <- 1L:nv
  for(i in 1L:(nu-1L)){
    i_nv <- i*nv
    k1 <- i_nv - nv
    rg <- (k1 + 1L):i_nv
    cosu_i <- cosu_[i]
    sinu_i <- sinu_[i]
    vs[, rg] <- rbind(
      cosu_i * Rrcosv_,
      sinu_i * Rrcosv_,
      rsinv_
    )
    normals[rg, ] <- cbind(
      cosu_i * cosv_,
      sinu_i * cosv_,
      sinv_
    )
    k_ <- k1 + j_
    l_ <- k1 + jp1_
    m_ <- i_nv + j_
    tris1[, k_] <- rbind(m_, l_, k_)
    tris2[, k_] <- rbind(m_, i_nv + jp1_, l_)
  }
  i_nv <- nunv
  k1 <- i_nv - nv
  rg <- (k1 + 1L):i_nv
  vs[, rg] <- rbind(
    Rrcosv_,
    0,
    rsinv_
  )
  normals[rg, ] <- cbind(
    cosv_,
    0,
    sinv_
  )
  l_ <- k1 + jp1_
  k_ <- k1 + j_
  tris1[, k_] <- rbind(j_, l_, k_)
  tris2[, k_] <- rbind(j_, jp1_, l_)
  tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    normals = normals,
    homogeneous = FALSE
  )
}

#' @title Cyclide mesh
#' @description Triangle mesh of a Dupin cyclide.
#'
#' @param a,c,mu cyclide parameters, positive numbers such that
#'   \code{c < mu < a}
#' @param nu,nv numbers of subdivisions, integers (at least 3)
#'
#' @return A triangle \strong{rgl} mesh (class \code{mesh3d}).
#'
#' @details The Dupin cyclide in the plane \emph{z=0}:
#' \if{html}{
#'   \figure{cyclide.png}{options: style="max-width:75\%;" alt="cyclide"}
#' }
#' \if{latex}{
#'   \out{\begin{center}}\figure{cyclide.png}\out{\end{center}}
#' }
#'
#' @export
#'
#' @importFrom rgl tmesh3d
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' mesh <- cyclideMesh(a = 97, c = 32, mu = 57)
#' sphere <- sphereMesh(x = 32, y = 0, z = 0, r = 40)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "chartreuse")
#' wire3d(mesh)
#' shade3d(sphere, color = "red")
#' wire3d(sphere)
cyclideMesh <- function(a, c, mu, nu = 90L, nv = 40L){
  stopifnot(c > 0, a > mu, mu > c)
  stopifnot(nu >= 3, nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  vertices <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
  normals  <- matrix(NA_real_, nrow = nu*nv, ncol = 3L)
  b2 <- a * a - c * c;
  bb <- sqrt(b2 * (mu * mu - c * c))
  omega <- (a * mu + bb) / c
  Omega0 <- c(omega, 0, 0)
  inversion <- function(M) {
    OmegaM <- M - Omega0
    k <- c(crossprod(OmegaM))
    OmegaM / k + Omega0
  }
  h <- (c * c) / ((a - c) * (mu - c) + bb)
  r <- (h * (mu - c)) / ((a + c) * (mu - c) + bb)
  R <- (h * (a - c)) / ((a - c) * (mu + c) + bb)
  bb2 <- b2 * (mu * mu - c * c)
  denb1 <- c * (a*c - mu*c + c*c - a*mu - bb)
  b1 <- (a*mu*(c-mu)*(a+c) - bb2 + c*c + bb*(c*(a-mu+c) - 2*a*mu))/denb1
  denb2 <- c * (a*c - mu*c - c*c + a*mu + bb)
  b2 <- (a*mu*(c+mu)*(a-c) + bb2 - c*c + bb*(c*(a-mu-c) + 2*a*mu))/denb2
  omegaT <- (b1 + b2)/2
  OmegaT <- c(omegaT, 0, 0)
  tormesh <- torusMesh(R, r, nu, nv)
  rtnormals <- r * tormesh[["normals"]][1L:3L, ]
  xvertices <- tormesh[["vb"]][1L:3L, ] + OmegaT
  for(i in 1L:nu){
    k0 <- i * nv - nv
    for(j in 1L:nv){
      k <- k0 + j
      rtnormal <- rtnormals[, k]
      xvertex <- xvertices[, k]
      vertex <- inversion(xvertex)
      vertices[, k] <- vertex
      foo <- vertex - inversion(rtnormal + xvertex)
      normals[k, ] <- foo / sqrt(c(crossprod(foo)))
    }
  }
  tmesh3d(
    vertices    = vertices,
    indices     = tormesh[["it"]],
    normals     = normals,
    homogeneous = FALSE
  )
}

HopfTorusMeshHelper <- function(u, cos_v, sin_v, nlobes, A, alpha){
  B <- pi/2 - (pi/2 - A)*cos(u*nlobes)
  C <- u + A*sin(2*u*nlobes)
  y1 <- 1 + cos(B)
  y23 <- sin(B) * c(cos(C), sin(C))
  y2 <- y23[1L]
  y3 <- y23[2L]
  x1 <- cos_v*y3 + sin_v*y2
  x2 <- cos_v*y2 - sin_v*y3
  x3 <- sin_v*y1
  x4 <- cos_v*y1
  yden <- sqrt(2*y1)
  if(is.null(alpha)){
    t(cbind(x1, x2, x3) / (yden-x4))
  }else{
    t(acos(x4/yden)/(yden^alpha-abs(x4)^alpha)^(1/alpha) * cbind(x1, x2, x3))
  }
}

#' @title Hopf torus mesh
#' @description Triangle mesh of a Hopf torus.
#'
#' @param nlobes number of lobes of the Hopf torus, a positive integr
#' @param A parameter of the Hopf torus, number strictly between
#'   \code{0} and \code{pi/2}
#' @param alpha if not \code{NULL}, this is the exponent of a modified
#'   stereographic projection, a positive number; otherwise the ordinary
#'   stereographic projection is used
#' @param nu,nv numbers of subdivisions, integers (at least 3)
#'
#' @return A triangle \strong{rgl} mesh (class \code{mesh3d}).
#' @export
#'
#' @importFrom rgl tmesh3d addNormals
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' mesh <- HopfTorusMesh(nu = 90, nv = 90)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "forestgreen")
#' wire3d(mesh)
#' mesh <- HopfTorusMesh(nu = 90, nv = 90, alpha = 1.5)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(0, 0, zoom = 0.75)
#' shade3d(mesh, color = "yellowgreen")
#' wire3d(mesh)
HopfTorusMesh <- function(
    nlobes = 3, A = 0.44, alpha = NULL, nu, nv
){
  stopifnot(isStrictPositiveInteger(nlobes))
  stopifnot(isPositiveNumber(A))
  stopifnot(A < pi/2)
  stopifnot(is.null(alpha) || isPositiveNumber(alpha))
  stopifnot(nu >= 3, nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  vs    <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
  tris1 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
  tris2 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
  u_ <- seq(0, 2*pi, length.out = nu + 1L)[-1L]
  v_ <- seq(0, 2*pi, length.out = nv + 1L)[-1L]
  cos_v <- cos(v_)
  sin_v <- sin(v_)
  jp1_ <- c(2L:nv, 1L)
  j_ <- 1L:nv
  for(i in 1L:(nu-1L)){
    i_nv <- i*nv
    vs[, (i_nv - nv + 1L):i_nv] <-
      HopfTorusMeshHelper(u_[i], cos_v, sin_v, nlobes, A, alpha)
    k1 <- i_nv - nv
    k_ <- k1 + j_
    l_ <- k1 + jp1_
    m_ <- i_nv + j_
    tris1[, k_] <- rbind(k_, l_, m_)
    tris2[, k_] <- rbind(l_, i_nv + jp1_, m_)
  }
  i_nv <- nu*nv
  vs[, (i_nv - nv + 1L):i_nv] <-
    HopfTorusMeshHelper(2*pi, cos_v, sin_v, nlobes, A, alpha)
  k1 <- i_nv - nv
  k_ <- k1 + j_
  l_ <- k1 + jp1_
  tris1[, k_] <- rbind(k_, l_, j_)
  tris2[, k_] <- rbind(l_, jp1_, j_)
  addNormals(tmesh3d(
    vertices    = vs,
    indices     = cbind(tris1, tris2),
    normals     = NULL,
    homogeneous = FALSE
  ))
}
