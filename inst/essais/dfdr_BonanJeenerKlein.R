f <- function(u, v)
  rbind(
    4 * cos(u) - cos(4*u) - 0.75 * (sin(3*u) + 1.5) * sin(2.5*u) * cos(v),
    (sin(3 * u) + 1.5) * sin(v),
    4 * sin(u) - sin(4*u) + 0.75 * (sin(3*u) + 1.5) * cos(2.5*u) * cos(v)
  )

fx <- function(u, v)
  4 * cos(u) - cos(4*u) - 0.75 * (sin(3*u) + 1.5) * sin(2.5*u) * cos(v)
fy <- function(u, v)
  (sin(3 * u) + 1.5) * sin(v)
fz <- function(u, v)
  4 * sin(u) - sin(4*u) + 0.75 * (sin(3*u) + 1.5) * cos(2.5*u) * cos(v)

library(dfdr)

dfx <- gradient(fx, FALSE, u, v)
dfy <- gradient(fy, FALSE, u, v)
dfz <- gradient(fz, FALSE, u, v)

normal <- Vectorize(function(u, v) {
  dx <- dfx(u, v)
  dy <- dfy(u, v)
  dz <- dfz(u, v)
  v1 <- c(dx[1L], dy[1L], dz[1L])
  v2 <- c(dx[2L], dy[2L], dz[2L])
  cgalMeshes:::crossProduct(v1, v2)
})

nu <- 100L
nv <- 10L

u_ <- seq(0, 2*pi, length.out = nu)#[-1L]
v_ <- seq(0, 2*pi, length.out = nv+1L)[-1L]
Grid <- expand.grid(U = u_, V = v_)
varray <- with(Grid, array(normal(U, V), dim = c(3L, nu, nv)))
varray2 <- aperm(varray, c(1L, 3L, 2L))
normals <- matrix(varray2, nrow = 3L, ncol = nu*nv)

rmesh <- parametricMesh(
  f, c(0, 2*pi), c(0, 2*pi), c(FALSE, TRUE), nu, nv
)

rmesh$normals <- -normals


rgl::shade3d(rmesh, color = "red")

