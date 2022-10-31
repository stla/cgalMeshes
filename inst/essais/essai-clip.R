# cube clipped to sphere ####
library(cgalMeshes)
library(rgl)
mesh    <- cgalMesh$new(cube3d())$triangulate()
clipper <- cgalMesh$new(sphereMesh(r= sqrt(2)))
mesh$clip(clipper, clipVolume = TRUE)
rglmesh <- mesh$getMesh(normals = FALSE)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(45, 45, zoom = 0.9)
shade3d(rglmesh, col = "darkorange")

# Togliatti surface clipped to a ball ####
library(rmarchingcubes)
library(rgl)
library(cgalMeshes)
# Togliatti surface equation: f(x,y,z) = 0
f <- function(x, y, z) {
  64*(x-1) *
    (x^4 - 4*x^3 - 10*x^2*y^2 - 4*x^2 + 16*x - 20*x*y^2 + 5*y^4 + 16 - 20*y^2) - 
    5*sqrt(5-sqrt(5))*(2*z - sqrt(5-sqrt(5))) * 
    (4*(x^2 + y^2 - z^2) + (1 + 3*sqrt(5)))^2
}
# grid
n <- 200L
x <- y <- seq(-5, 5, length.out = n)
z <- seq(-4, 4, length.out = n)
Grid <- expand.grid(X = x, Y = y, Z = z)
# calculate voxel
voxel <- array(with(Grid, f(X, Y, Z)), dim = c(n, n, n))
# calculate isosurface
contour_shape <- contour3d(
  griddata = voxel, level = 0, x = x, y = y, z = z
)
# make rgl mesh (plotted later)
rglMesh <- tmesh3d(
  vertices = t(contour_shape[["vertices"]]),
  indices  = t(contour_shape[["triangles"]]),
  normals  = contour_shape[["normals"]],
  homogeneous = FALSE
)
# make CGAL mesh
mesh <- cgalMesh$new(rglMesh)
# clip to sphere of radius 4.8
sphere <- sphereMesh(r = 4.8)
clipper <- cgalMesh$new(sphere)
mesh$clip(clipper, clipVolume = FALSE)
rglClippedMesh <- mesh$getMesh()
# plot
open3d(windowRect = 50 + c(0, 0, 900, 450))
mfrow3d(1L, 2L)
view3d(0, -70, zoom = 0.8)
shade3d(rglMesh, color = "firebrick")
next3d()
view3d(0, -70, zoom = 0.8)
shade3d(rglClippedMesh, color = "firebrick")
shade3d(sphere, color = "yellow", alpha = 0.2)
