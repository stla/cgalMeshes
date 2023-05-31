library(cgalMeshes)
library(rgl)

frog <- readSTL("Tree_Frog.stl", ascii = FALSE, plot = FALSE)
pts <- t(frog)

xptr <- cgalMeshes:::alphaWrap_cpp(pts, 50, 1000)

mesh <- cgalMesh$new(clean = xptr)
#mesh$LoopSubdivision(3)
#mesh$smoothShape(time = 0.01, iterations = 4)
mesh$computeNormals()


rmesh <- mesh$getMesh()

shade3d(rmesh, color = "red")


pmesh <- PoissonReconstruction(frog)
