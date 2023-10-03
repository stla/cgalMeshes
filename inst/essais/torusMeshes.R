torusMesh1 <- function(R, r, nu, nv, nrmls) {
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  nunv <- nu*nv
  vs      <- matrix(NA_real_, nrow = 3L, ncol = nunv)
  normals <- if(nrmls) matrix(NA_real_, nrow = nunv, ncol = 3L)
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
    if(nrmls) {
      normals[rg, ] <- cbind(
        cosu_i * cosv_,
        sinu_i * cosv_,
        sinv_
      )
    }
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
  if(nrmls) {
    normals[rg, ] <- cbind(
      cosv_,
      0,
      sinv_
    )
  }
  l_ <- k1 + jp1_
  k_ <- k1 + j_
  tris1[, k_] <- rbind(j_, l_, k_)
  tris2[, k_] <- rbind(j_, jp1_, l_)
  tmesh3d(
    vertices    = vs,
    indices     = cbind(tris1, tris2),
    normals     = normals,
    homogeneous = FALSE
  )
}

torusMesh2 <- function(R, r, nu, nv) {
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  nunv <- nu*nv
  vs    <- matrix(NA_real_, nrow = 3L, ncol = nunv)
  tris1 <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  tris2 <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  u_ <- seq(0, 2*pi, length.out = nu + 1L)[-1L]
  cosu_ <- cos(u_)
  sinu_ <- sin(u_)
  v_ <- seq(0, 2*pi, length.out = nv + 1L)[-1L]
  cosv_ <- cos(v_)
  sinv_ <- sin(v_)
  kxy <- R*R - r*r
  kz <- sqrt(kxy) * r
  kcosu_ <- kxy * cosu_
  ksinu_ <- kxy * sinu_
  w_ <- R - r * cosv_
  h_ <- kz * sinv_ / w_
  jp1_ <- c(2L:nv, 1L)
  j_ <- 1L:nv
  for(i in 1L:(nu-1L)){
    i_nv <- i*nv
    k1 <- i_nv - nv
    rg <- (k1 + 1L):i_nv
    kcosu_i <- kcosu_[i]
    ksinu_i <- ksinu_[i]
    vs[, rg] <- rbind(
      kcosu_i / w_,
      ksinu_i / w_,
      h_
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
    kxy / w_,
    0,
    h_
  )
  l_ <- k1 + jp1_
  k_ <- k1 + j_
  tris1[, k_] <- rbind(j_, l_, k_)
  tris2[, k_] <- rbind(j_, jp1_, l_)
  tmesh3d(
    vertices    = vs,
    indices     = cbind(tris1, tris2),
    homogeneous = FALSE
  )
}


torusMesh <- function(
  R, r, p1 = NULL, p2 = NULL, p3 = NULL, nu = 50, nv = 30, 
  normals = TRUE, conformal = FALSE
) {
  transformation <- !is.null(p1) && !is.null(p2) && !is.null(p3)
  if(transformation) {
    ccircle <- circumcircle(p1, p2, p3)
    R <- ccircle[["radius"]]
    # !! we have to take the inverse matrix for rgl::rotate3d
    rotMatrix <- rotationFromTo(ccircle[["normal"]], c(0, 0, 1))
    center <- ccircle[["center"]]
  }
  stopifnot(isPositiveNumber(R), isPositiveNumber(r))
  stopifnot(R > r)
  stopifnot(nu >= 3, nv >= 3)
  stopifnot(isBoolean(normals))
  stopifnot(isBoolean(conformal))
  if(conformal) {
    rmesh <- torusMesh2(R, r, nu, nv)
    if(transformation) {
      rmesh <- translate3d(
        rotate3d(
          rmesh, matrix = rotMatrix
        ),
        x = center[1L], y = center[2L], z = center[3L]
      )
    }
    if(normals) {
      rmesh <- addNormals(rmesh)
    }
  } else {
    rmesh <- torusMesh1(R, r, nu, nv, normals)
    if(transformation) {
      rmesh <- translate3d(
        rotate3d(
          rmesh, matrix = rotMatrix
        ),
        x = center[1L], y = center[2L], z = center[3L]
      )
    }
  }
  rmesh
}
