library(cgalMeshes)
library(rgl)
# twisted sphere
fx <- function(u, v) cos(u) * cos(v)
fy <- function(u, v) sin(u) * cos(v)
fz <- function(u, v) sin(v) + u
f <- function(u, v) {
  rbind(fx(u, v), fy(u, v), fz(u, v))
}
if(require("dfdr")) {
  # compute the gradient with the 'dfdr' package
  dfx <- gradient(fx, FALSE, u, v)
  dfy <- gradient(fy, FALSE, u, v)
  dfz <- gradient(fz, FALSE, u, v)
  fnormal <- Vectorize(function(u, v) {
    dx <- dfx(u, v)
    dy <- dfy(u, v)
    dz <- dfz(u, v)
    v1 <- c(dx[1L], dy[1L], dz[1L])
    v2 <- c(dx[2L], dy[2L], dz[2L])
    - crossProduct(v1, v2)
  })
} else {
  fnormal <- NULL
}
# compute mesh of the parametric surface
rmesh <- cgalMeshes::parametricMesh(
  f, c(0, 2*pi), c(0, 2*pi), c(FALSE, TRUE), nu = 60L, nv = 30L, fnormal = fnormal
)
# plot
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(50, 10)
shade3d(rmesh, color = "midnightblue", specular = "black")


