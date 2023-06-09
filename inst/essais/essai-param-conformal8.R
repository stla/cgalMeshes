library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")
library(rgl)

source("~/Documents/R/rgl/parametricSurfaces/HopfTorusSmallLobes.R")

rmesh <- HopfTorusMesh(200, 200)
mesh <- cgalMesh$new(rmesh)
mesh$clipToPlane(c(0,0,0), c(0,0,1), FALSE)
rmesh <- mesh$getMesh()
rmesh <- Rvcg::vcgClean(rmesh, sel = 0)
mesh <- cgalMesh$new(rmesh)
str(mesh$getBorders())
mesh$fillBoundaryHole(2, fair = FALSE)
summary(mesh$getEdges())
mesh$isotropicRemeshing(0.015, iterations = 3, relaxSteps = 2)
mesh$writeMeshFile("torus.off")

M <- cgalMeshes:::testparam(normalizePath("torus.off"), 1L)

clrs <- ifelse(
  (floor(10*M[, 1L]) %% 2) == (floor(10*M[, 2L]) %% 2), 
  "yellow", "navy"
)
plot(M, type = "p", asp = 1, pch = ".", col=clrs, xlab = "u", ylab = "v")

mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(180, 0, zoom = 0.85)
shade3d(rmesh, meshColor = "vertices")

movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "half-HopfTorus-DCP.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

