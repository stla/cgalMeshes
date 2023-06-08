library(cgalMeshes)
library(rgl)
library(rmarchingcubes)

# the pseudo-gyroid is the isosurface f(x,y,z)=0
f <- function(x, y, z) {
  cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x)
}
# construct the isosurface f=0
ngrid <- 70L
x <- y <- z <- seq(-5, 5, length.out = ngrid)
Grid <- expand.grid(X = x, Y = y, Z = z)
voxel <- array(
  with(Grid, f(X, Y, Z)), dim = c(ngrid, ngrid, ngrid)
)
library(rmarchingcubes)
contour_shape <- contour3d(
  griddata = voxel, level = 0, x = x, y = y, z = z
)
# make mesh
library(rgl)
rglMesh <- tmesh3d(
  vertices = t(contour_shape[["vertices"]]),
  indices  = t(contour_shape[["triangles"]]),
  normals  = contour_shape[["normals"]]
)

open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-40, 35)
shade3d(rglMesh, color = "orange")
snapshot3d("figures/pseudogyroid_cube.png")


bdry <- getBoundary3d(rglMesh, color = "black", lwd = 3)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-40, 35)
shade3d(rglMesh, color = "orange")
shade3d(bdry)
snapshot3d("figures/pseudogyroid_cube_boundary.png")

# returns the squared norms of the vertices
sqnorm <- function(vertices) { 
  apply(vertices, 1L, function(row) crossprod(row))
}
# clipping
rglMesh2 <- clipMesh3d(rglMesh, sqnorm, bound = 25, greater = FALSE)
bdry <- getBoundary3d(rglMesh2, color = "black", lwd = 3)

open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.7)
shade3d(rglMesh2, color = "orangered")
shade3d(bdry)
snapshot3d("figures/pseudogyroid_sphere.png")


rglMesh2

vertices <- t(rglMesh2$vb[-4L, ])




# reconstruct the mesh from its vertices
wrapMesh <- alphaWrap(vertices, ralpha = 100, roffset = 1000) 
wrapMesh$computeNormals()
wrapMesh
# plot
wrapRglMesh <- wrapMesh$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.7)
shade3d(wrapRglMesh, color = "orangered")
snapshot3d("figures/pseudogyroid_wrapped.png")

g <- function(xyz) f(xyz[1L], xyz[2L], xyz[3L])
library(numDeriv)
normals <- apply(wrapRglMesh$vb[-4L, ], 2L, function(v) {
  grad(g, v)
})
vs <- wrapMesh$getVertices()
ff <- apply(vs, 1L, g)
eps <- 0
normals[, ff < eps] <- - normals[, ff < eps]

wrapRglMesh$normals <- normals

