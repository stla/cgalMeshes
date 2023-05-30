library(rgl)
library(cgalMeshes)

data("teapot", package = "misc3d")

pts <- teapot$vertices

xptr <- cgalMeshes:::alphaWrap_cpp(pts, 30, 600)

mesh <- cgalMesh$new(clean = xptr)

rmesh <- mesh$getMesh()

shade3d(rmesh, color = "red")

