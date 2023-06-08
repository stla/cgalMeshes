library(cgalMeshes)
library(rgl)
# take the diplodocus mesh
off <- system.file("extdata", "diplodocus.off", package = "cgalMeshes")
diplodocusMesh <- cgalMesh$new(off)
diplodocusMesh$computeNormals()
diplodocusRglMesh <- diplodocusMesh$getMesh()

open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(20, 0, zoom = 0.85)
shade3d(diplodocusRglMesh, color = "forestgreen")

vertices <- diplodocusMesh$getVertices()
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(20, 0, zoom = 0.85)
points3d(vertices)

AFSmesh <- AFSreconstruction(vertices)
AFSmesh$computeNormals()
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(20, 0, zoom = 0.85)
shade3d(AFSmesh$getMesh(), color = "forestgreen")


# reconstruct the mesh from its vertices
pts <- diplodocusMesh$getVertices()
wrapMesh <- alphaWrap(pts, 100, 1000) -> wrapMesh0
# shade3d(Rvcg::vcgUpdateNormals(wrapMesh$getMesh(), type = 0), color = "orange")
# wrapMesh$smoothShape(time = 0.001, iterations = 10)
wrapMesh$computeNormals()
wrapMesh
# plot
wrapRglMesh <- wrapMesh$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512))
shade3d(wrapRglMesh, color = "forestgreen")


open3d(windowRect = 50 + c(0, 0, 800, 400))
mfrow3d(1, 2)
view3d(20, 0, zoom = 0.85)
shade3d(diplodocusRglMesh, color = "forestgreen")
next3d()
view3d(20, 0, zoom = 0.85)
shade3d(wrapRglMesh, color = "forestgreen")


pmesh <- PoissonReconstruction(vertices)
pmesh$computeNormals()
rmesh <- pmesh$getMesh()
shade3d(rmesh, color = "red")

library(cgalMeshes)
AFSmesh <- AFSreconstruction(vertices)
AFSmesh$computeNormals()
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.7)
shade3d(AFSmesh$getMesh(), color = "red")
snapshot3d("figures/pseudogyroid_AFS.png")

library(AlphaHull3D)
ahull <- fullAhull3d(vertices)
amesh <- setAlpha(ahull, alpha = 1.8)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.7)
shade3d(amesh, color = "green")
snapshot3d("figures/pseudogyroid_alphaHull.png")
