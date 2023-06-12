library(rgl)
library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

# Barth sextic
# isosurface f=0
phi <- (1 + sqrt(5)) / 2
f <- function(x, y, z){
  4 * (phi^2*x^2 - y^2) * (phi^2*y^2 - z^2) * (phi^2*z^2 - x^2) - 
    (1 + 2*phi) * (x^2 + y^2 + z^2 - 1)^2
}
# make the isosurface
nx <- 220L; ny <- 220L; nz <- 220L
x <- seq(-1.8, 1.8, length.out = nx) 
y <- seq(-1.8, 1.8, length.out = ny)
z <- seq(-1.8, 1.8, length.out = nz) 
Grid <- expand.grid(X = x, Y = y, Z = z)
voxel <- array(with(Grid, f(X, Y, Z)), dim = c(nx, ny, nz))
mask  <- array(with(Grid, X^2 + Y^2 + Z^2 > 3), dim = c(nx, ny, nz))
voxel[mask] <- -1000
library(rmarchingcubes)
cont <- contour3d(voxel, level = 0, x = x, y = y, z = z)
# plot
library(rgl)
rmesh <- tmesh3d(
  vertices = t(cont[["vertices"]]),
  indices  = t(cont[["triangles"]]),
  normals  = cont[["normals"]],
  homogeneous = FALSE
)

#
mesh <- cgalMesh$new(rmesh)
mesh$writeMeshFile("barth.off")

vs <- mesh$getVertices()
dists <- apply(vs, 1L, crossprod)
which(dists == max(dists))

