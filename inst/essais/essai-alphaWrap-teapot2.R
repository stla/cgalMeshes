library(rgl)
library(cgalMeshes)

data("teapot", package = "misc3d")
pts <- teapot$vertices

teapotmesh <- readOBJ("teapot.obj")
mesh <- cgalMesh$new(teapotmesh)
mesh$LoopSubdivision(2)
pts <- t(mesh$getVertices())
dim(pts)

xptr <- cgalMeshes:::alphaWrap_cpp(pts, 40, 500)

mesh <- cgalMesh$new(clean = xptr)
#mesh$LoopSubdivision(3)
#mesh$smoothShape(time = 0.01, iterations = 4)
mesh$computeNormals()


rmesh <- mesh$getMesh()

shade3d(rmesh, color = "red")

