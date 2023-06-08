library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")
library(rgl)

alpha <- 0.15
beta <- 1
gamma <- 0.1
n <- 2
f <- function(u, v) {
  rbind(
    alpha * (1 - v/(2*pi)) * cos(n*v) * (1 + cos(u)) + gamma * cos(n * v),
    alpha * (1 - v/(2*pi)) * sin(n*v) * (1 + cos(u)) + gamma * sin(n * v),
    alpha * (1 - v/(2*pi)) * sin(u) + beta * v/(2*pi)
  )
}

rmesh <- parametricMesh(
  f, c(0, 2*pi), c(0, 2*pi), periodic = c(TRUE, FALSE), nu = 10, nv = 10
)

rmesh <- Rvcg::vcgClean(rmesh, sel = 0)
rmesh$normals <- NULL
# Morpho::plotNormals(rmesh)
# rmesh <- Rvcg::vcgUniformRemesh(rmesh, discretize = TRUE, voxelSize = 5e-3, offset = 0)

mesh <- cgalMesh$new(rmesh)
summary(mesh$getEdges())
#mesh$isotropicRemeshing(2.5e-3, iterations = 3, relaxSteps = 2)
mesh$writeMeshFile("torus.off")

M <- cgalMeshes:::testparam(normalizePath("torus.off"), 3L)

stop()

clrs <- ifelse(
  (floor(10*M[, 1L]) %% 2) == (floor(10*M[, 2L]) %% 2), 
  "yellow", "navy"
)
plot(M, type = "p", asp = 1, pch = ".", col=clrs, xlab = "u", ylab = "v")


mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
shade3d(rmesh, meshColor = "vertices")

movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "twisted-horn-DAP.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

snapshot3d("Enneper_colored.png")
