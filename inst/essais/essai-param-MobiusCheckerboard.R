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

#### Möbius stuff
M_power_t <- function(gamma = 0.2 - 0.1i, t){
  h <- sqrt(1-Mod(gamma)^2)
  d2 <- h^t * (cospi(t/2) + 1i*sinpi(t/2))
  d1 <- Conj(d2)
  a <- Re(d1) - 1i*Im(d1)/h
  b <- gamma * Im(d2)/h
  c <- Conj(b)
  d <- Conj(a)
  c(a = a, b = b, c = c, d = d)
}

Mobius <- function(xy, abcd) {
  a <- abcd["a"]
  b <- abcd["b"]
  c <- abcd["c"]
  d <- abcd["d"] 
  z <- complex(real = xy[1L], imaginary = xy[2L])
  Mz <- (a*z + b) / (c*z + d)
  c(Re(Mz), Im(Mz))
}

#### radial checkerboard
M0 <- M
abcd <- M_power_t(t = 1.5)
M <- t(apply(10 * (M0 - 0.5), 1L, Mobius, abcd = abcd))
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
fclrs <- function(t) {
  abcd <- M_power_t(t = t)
  M <- t(apply(10 * (M0 - 0.5), 1L, Mobius, abcd = abcd))
  radii <- sqrt(apply(M, 1L, crossprod))
  angles <- 10 * (1 + atan2(M[, 2L], M[, 1L])/pi)
  tests <- floor(angles) %% 2 == 0
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

alpha_ <- seq(0, 2, length.out = 61L)[-1L]
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
  gif_file = "Enneper-MobiusCheckerboard-DCP.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

