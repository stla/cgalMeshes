library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")
library(rgl)

# Enneper ####
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
  periodic = c(TRUE, FALSE), nu = 400, nv = 200, clean = FALSE
)
rmesh <- Rvcg::vcgClean(rmesh, sel = 0)
mesh <- cgalMesh$new(rmesh)
mesh$isotropicRemeshing(0.01, iterations = 3, relaxSteps = 2)


mesh$writeMeshFile("torus.off")

M <- cgalMeshes:::testparam(normalizePath("torus.off"), 1L)

#### radial checkerboard
M0 <- M
M <- 10 * (M0 - 0.5)
radii <- sqrt(apply(M, 1L, crossprod))
angles <- 10 * (1 + atan2(M[, 2L], M[, 1L])/pi)
clrs <- ifelse(
  floor(radii) %% 2 == 0,
  ifelse(
    floor(angles) %% 2 == 0, "navy", "yellow"
  ),
  ifelse(
    floor(angles) %% 2 == 0, "yellow", "navy"
  )
)

plot(M0, type = "p", asp = 1, pch = ".", col=clrs, xlab = "u", ylab = "v")


mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, -20, zoom = 0.8)
shade3d(rmesh, meshColor = "vertices")

# animation ####
angles0 <- angles
fclrs <- function(alpha) {
  tests <- floor(angles0 + alpha) %% 2 == 0
  ifelse(
    floor(radii) %% 2 == 0,
    ifelse(
      tests, "navy", "yellow"
    ),
    ifelse(
      tests, "yellow", "navy"
    )
  )
}

alpha_ <- seq(0, 2, length.out = 19L)[-1L]
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, -20, zoom = 0.8)
for(i in seq_along(alpha_)) {
  clrs <- fclrs(alpha_[i])
  rmesh$material <- list(color = clrs)
  shade3d(rmesh, meshColor = "vertices")
  snapshot3d(sprintf("zzpic%03d.png", i), webshot = FALSE)
  clear3d()
}


library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "Enneper-radialCheckerboard-DCP.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

appendGIFs("Enneper-nopsi-DCP.gif", "Enneper-psi-DCP.gif", delay = 8)
