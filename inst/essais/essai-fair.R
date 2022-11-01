rglhopf <- MeshesTools::HopfTorusMesh(nu = 100, nv = 100)
#rgl::shade3d(rglhopf, color = "orangered")

library(cgalMeshes)
hopf <- cgalMesh$new(rglhopf)
vs <- hopf$vertices()
summary(vs)

normsq <- apply(vs, 1L, crossprod)
indices <- which(normsq > 19)
hopf$fair(indices)

rglm <- hopf$getMesh()

library(rgl)
open3d(windowRect = c(50, 50, 950, 500))
mfrow3d(1L, 2L)
view3d(0, 0, zoom = 0.8)
shade3d(rglhopf, color = "orangered")
next3d()
view3d(0, 0, zoom = 0.8)
shade3d(rglm, color = "orangered")


open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, 0, zoom = 0.8)
shade3d(rglm, color = "orangered")


