library(cgalMeshes)
rmesh <- sphereMesh(iterations = 8)
mesh <- cgalMesh$new(rmesh)
mesh$clipToPlane(c(0,0,0), c(1,0,0), FALSE)
rmesh <- Rvcg::vcgClean(mesh$getMesh(), 0:7)
mesh <- cgalMesh$new(rmesh)
mesh$writeMeshFile("torus.off")
M <- cgalMeshes:::testparam()
clrs <- rep("yellow", nrow(M))
clrs[M[, 1] < 0.5 & M[, 2] < 0.5] <- "navy"
clrs[M[, 1] > 0.5 & M[, 2] > 0.5] <- "navy"
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

rgl::shade3d(rmesh, meshColor = "vertices")
rgl::snapshot3d("halftorus2_colored.png")
