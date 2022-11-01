library(cgalMeshes)
rglHopf <- HopfTorusMesh(nu = 100, nv = 100)
hopf <- cgalMesh$new(rglhopf)
# squared norms of the vertices
normsq <- apply(hopf$vertices(), 1L, crossprod)
# fair the region where the squared norm is > 19
indices <- which(normsq > 19)
hopf$fair(indices)
rglHopf_faired <- hopf$getMesh()
# plot
library(rgl)
open3d(windowRect = 50 + c(0, 0, 900, 450))
mfrow3d(1L, 2L)
view3d(0, 0, zoom = 0.8)
shade3d(rglHopf, color = "orangered")
next3d()
view3d(0, 0, zoom = 0.8)
shade3d(rglHopf_faired, color = "orangered")
