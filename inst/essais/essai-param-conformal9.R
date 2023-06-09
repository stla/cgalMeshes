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
  f, c(0, 2*pi), c(0, 2*pi), periodic = c(TRUE, FALSE), nu = 100, nv = 100, clean = TRUE
)

mesh <- cgalMesh$new(rmesh)
summary(mesh$getEdges())
mesh$isotropicRemeshing(2.5e-3, iterations = 3, relaxSteps = 2)

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

#### conformal map circle -> square
library(Carlson) # to get the `elliptic_F` function
# define the `psi` function
w <- sqrt(1i)
D <- elliptic_F(1i * asinh(w), -1)
psi <- function(z) { # maps the circle to the square
  elliptic_F(1i * asinh(w * z), -1, minerror = 1e-10) / D
}
#psi <- identity
colormap <- function(xy) {
  z <- psi(complex(real = xy[1L], imaginary = xy[2L]))
  x <- 3 * Re(z)
  y <- 3 * Im(z)
  ifelse(
    x %% 2 <= 1,
    ifelse(
      y %% 2 <= 1, "navy", "yellow"
    ),
    ifelse(
      y %% 2 <= 1, "yellow", "navy"
    )
  )
}

clrs <- apply(1.999*(M-0.5), 1L, colormap)

par(mar = c(2, 2, 1, 1), mfrow = c(1, 2))
plot(M, type = "p", asp = 1, pch = ".", col=clrs1, xlab = "u", ylab = "v")
plot(M, type = "p", asp = 1, pch = ".", col=clrs, xlab = "u", ylab = "v")


mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
shade3d(rmesh, meshColor = "vertices")

movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "Enneper-nopsi-DCP.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

appendGIFs("Enneper-nopsi-DCP.gif", "Enneper-psi-DCP.gif", delay = 8)
