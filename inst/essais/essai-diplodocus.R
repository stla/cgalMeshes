library(cgalMeshes)
library(rgl)

mesh <- cgalMesh$new("diplodocus.off")
mesh$computeNormals()

rmesh <- mesh$getMesh()

shade3d(rmesh, color = "forestgreen")

pts <- t(mesh$getVertices())

xptr <- cgalMeshes:::alphaWrap_cpp(pts, 70, 3000)

wmesh <- cgalMesh$new(clean = xptr)
#wmesh$LoopSubdivision(1)
#mesh$smoothShape(time = 0.01, iterations = 1)
wmesh$computeNormals()

rmesh <- wmesh$getMesh()
rmesh
open3d()
shade3d(rmesh, color = "blue")


pmesh <- PoissonReconstruction(mesh$getVertices(), spacing = 0.005, sm_distance = 0.2)
pmesh$computeNormals()
rmesh <- pmesh$getMesh()
shade3d(rmesh, color = "red")

amesh <- AFSreconstruction(mesh$getVertices())
amesh$computeNormals()
rmesh <- amesh$getMesh()
shade3d(rmesh, color = "red")

library(AlphaHull3D)
ahull <- fullAhull3d(t(pts))
rmesh <- setAlpha(ahull, alpha = 0.005)

shade3d(rmesh, color = "green")
