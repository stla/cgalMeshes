library(rgl)
library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

# soccer star
# isosurface f=0
a1 <- a2 <- a3 <- a4 <- -100
f <- function(x, y, z){
  u <- x*x + y*y + z*z
  v <- -z*(2*x+z)*(x^4 - x*x*z*z + z^4 + 2*(x^3*z-x*z^3) + 5*(y^4-y*y*z*z) + 10*(x*y*y*z-x*x*y*y))
  w <- (4*x*x + z*z - 6*x*z) * 
    (z^4 - 2*z^3*x - x*x*z*z + 2*z*x^3 + x^4 - 25*y*y*z*z - 30*x*y*y*z - 10*x*x*y*y + 5*y^4) *
    (z^4 + 8*z^3*x + 14*x*x*z*z - 8*z*x^3 + x^4 - 10*y*y*z*z - 10*x*x*y*y + 5*y^4)
  1 + ((128565+115200*sqrt(5))/1295029 * a3 + (49231296000*sqrt(5)-93078919125)/15386239549 * a4 - a1 - 3*a2 - 3) * u +
    ((-230400*sqrt(5) - 257130)/1295029 * a3 + (238926989250-126373248000*sqrt(5))/15386239549 * a4 + 3*a1 + 8*a2 + 3) * u * u + 
    ((115200*sqrt(5)+128565)/1295029 * a3 + (91097280000*sqrt(5)-172232645625)/15386239549 * a4 - 3*a1 - 6*a2 - 1) * u*u*u + 
    (a3 + (121075-51200*sqrt(5))/11881 * a4) * v + ((102400*sqrt(5)-242150)/11881 - 2*a3) * u * v + 
    a1 * u^4 + a2 * u^5 + a3 * u*u*v + a4*w
}
# make the isosurface
nx <- 300L; ny <- 300L; nz <- 300L
x <- seq(-1.5, 1.5, length.out = nx) 
y <- seq(-1.5, 1.5, length.out = ny)
z <- seq(0, 1.5, length.out = nz) 
Grid <- expand.grid(X = x, Y = y, Z = z)
voxel <- array(with(Grid, f(X, Y, Z)), dim = c(nx, ny, nz))
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

view3d(0, 0)
shade3d(rmesh, color = "hotpink")

#
mesh <- cgalMesh$new(rmesh)
mesh$writeMeshFile("soccer.off")



UV <- cgalMeshes:::testo("soccer.off")
plot(UV, type="p", asp = 1, pch = ".", cex = 2)

vs <- mesh$getVertices()
dists <- apply(vs, 1L, crossprod)
ii <- which(dists > max(dists)-0.0005)

points3d(vs)
points3d(vs[ii,], col="red", size=10)

# 1 14 428129
points3d(vs[c(1, 14, 428129),], col="red", size=11)

open3d()
points3d(vs)
points3d(vs[ii,], col="red", size=10)
f <- select3d()
if (!is.null(f)) {
  keep <- f(x, y, z)
  pop3d()
  points3d(x[keep], y[keep], z[keep], color = 'green', size=10)
  points3d(x[!keep], y[!keep], z[!keep])
}