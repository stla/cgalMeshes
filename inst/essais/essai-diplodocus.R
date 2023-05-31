library(cgalMeshes)
library(rgl)
# take the diplodocus mesh
off <- system.file("extdata", "diplodocus.off", package = "cgalMeshes")
diplodocusMesh <- cgalMesh$new(off)
diplodocusMesh$computeNormals()
# reconstruct the mesh from its vertices
pts <- diplodocusMesh$getVertices()
wrapMesh <- alphaWrap(pts, 70, 3000)
wrapMesh$computeNormals()
# plot
diplodocusRglMesh <- diplodocusMesh$getMesh()
wrapRglMesh <- wrapMesh$getMesh()
open3d(windowRect = 50 + c(0, 0, 800, 400))
mfrow3d(1, 2)
view3d(20, 0, zoom = 0.85)
shade3d(diplodocusRglMesh, color = "forestgreen")
next3d()
view3d(20, 0, zoom = 0.85)
shade3d(wrapRglMesh, color = "forestgreen")




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
