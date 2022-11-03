library(cgalMeshes)
library(rmarchingcubes)
# isosurface function (slice of a seven-dimensional toratope)
f <- function(x, y, z, a) {
  (sqrt(
    (sqrt((sqrt((x*sin(a))^2 + (z*cos(a))^2) - 5)^2 + (y*sin(a))^2) - 2.5)^2 + 
      (x*cos(a))^2) - 1.25
  )^2 + (sqrt((sqrt((z*sin(a))^2 + (y*cos(a))^2) - 2.5)^2) - 1.25)^2
}
# make grid
n <- 200L
x <- seq(-10, 10, len = n)
y <- seq(-10, 10, len = n)
z <- seq(-10, 10, len = n)
Grid <- expand.grid(X = x, Y = y, Z = z)
# compute isosurface
voxel <- array(with(Grid, f(X, Y, Z, a = pi/2)), dim = c(n, n, n))
isosurface <- contour3d(voxel, level = 0.25, x = x, y = y, z = z)
# make CGAL mesh
mesh <- cgalMesh$new(
  vertices = isosurface[["vertices"]], faces = isosurface[["triangles"]]
)
# connected components
components <- mesh$connectedComponents()
ncc <- length(components)
# plot
library(rgl)
colors <- rainbow(ncc)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(30, 50)
for(i in 1L:ncc) {
  rglMesh <- components[[i]]$getMesh()
  shade3d(rglMesh, color = colors[i])
}



mesh <- tmesh3d(
  vertices = t(cont[["vertices"]]),
  indices  = t(cont[["triangles"]]),
  normals  = cont[["normals"]],
  homogeneous = FALSE
)

open3d(windowRect = c(50, 50, 562, 562), zoom = 0.85)
view3d(30, 50)
shade3d(mesh, color = "turquoise")


