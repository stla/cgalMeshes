library(rgl)
library(cgalMeshes)

N <- 1024

g <- function(phi) {
  r <- sin(phi) / sqrt(1 - sin(phi)^2)
  R <- cos(phi) / (1 - sin(phi)) - r
  h <- R/r
  sqrt(h*h - 1)
}

phi <- pi/4 # c'est la root de g=1
r <- sin(phi) / sqrt(1 - sin(phi)^2)
R <- cos(phi) / (1 - sin(phi)) - r

torus <- function(R, r) {
  h <- R/r
  s <- sqrt(h*h - 1)
  r <- 1/s/r
  f <- function(u, v) {
    w <- h - cospi(2*v)
    rbind(
      s * cospi(2*u/s) / w,
      s * sinpi(2*u/s) / w,
      sinpi(2*v) / w
    ) / r
  }
  parametricMesh(
    f, c(0, s), c(0, 1), c(TRUE, TRUE),
    nu = N, nv = N
  )
}

mesh <- torus(7, 3)
summary(t(mesh$vb))
# r = 1/sqrt((h-1)*(h+1))


u_ <- seq(0, 1, length.out = N)
v_ <- seq(0, 1, length.out = N)
UV <- as.matrix(
  expand.grid(V = v_, U = u_)
)

rot <- function(alpha, uv) {
  t(rbind(
    c(cos(alpha), -sin(alpha)),
    c(sin(alpha),  cos(alpha))
  ) %*% t(uv))
}

UVrot <- rot(pi/4, UV)

clrs1 <- ifelse(
  (floor(4*sqrt(2)*UVrot[, 1L]) %% 2) == (floor(4*sqrt(2)*UVrot[, 2L]) %% 2),
  "yellow", "navy"
)
plot(UV, col = clrs1, asp = 1, pch = ".")

mesh$material <- list(color = clrs1)
shade3d(mesh, polygon_offset = 1) 
#bbox3d()

Villarceau <- function(beta, theta0 = 0) {
  d <- (1 - sin(beta) * sin(phi))
  cbind(
    cos(theta0 + beta) * cos(phi) / d,
    sin(theta0 + beta) * cos(phi) / d,
    cos(beta) * sin(phi) / d
  ) 
}

beta_ <- seq(0, 2*pi, len = 400)
pts <- Villarceau(beta_)

points3d(pts, size = 7)