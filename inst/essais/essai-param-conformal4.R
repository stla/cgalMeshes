library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

n <- 3
Enneper <- function(phi, r) {
  rbind(
    r*cos(phi) - r^(2*n-1)*cos((2*n-1)*phi)/(2*n-1),
    r*sin(phi) + r^(2*n-1)*sin((2*n-1)*phi)/(2*n-1),
    2*r^n*cos(n*phi)/n
  )
}

rmesh <- parametricMesh(
  Enneper, urange = c(0, 2*pi), vrange = c(0, 1.3),
  periodic = c(TRUE, FALSE), nu = 400, nv = 200
)
rmesh <- Rvcg::vcgClean(rmesh, sel = 0)
# rmesh <- Rvcg::vcgIsotropicRemeshing(rmesh, TargetLen = 0.01)
# rmesh <- Rvcg::vcgUpdateNormals(rmesh)
# rgl::shade3d(rmesh, color = "red")

mesh <- cgalMesh$new(rmesh)
mesh$isotropicRemeshing(0.01, iterations = 3, relaxSteps = 2)
mesh$writeMeshFile("torus.off")

M <- cgalMeshes:::testparam()

dists <- sqrt(apply(2*M-1, 1L, crossprod))

clrs <- rep("yellow", nrow(M))
clrs[dists <= 0.2] <- "navy"
clrs[dists > 0.4 & dists <= 0.6] <- "navy"
clrs[dists > 0.8] <- "navy"
plot(M, type = "p", asp = 1, pch = ".", col=clrs, xlab = "u", ylab = "v")


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
  gif_file = "Enneper_colored_circular.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

snapshot3d("Enneper_colored.png")
