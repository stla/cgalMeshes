library(cgalMeshes)
library(rgl)
library(rmarchingcubes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

# Togliatti surface equation with spherical coordinates
f <- function(x,y,z){
  64*(x-1)*
    (x^4-4*x^3-10*x^2*y^2-4*x^2+16*x-20*x*y^2+5*y^4+16-20*y^2) - 
    5*sqrt(5-sqrt(5))*(2*z-sqrt(5-sqrt(5)))*(4*(x^2+y^2-z^2)+(1+3*sqrt(5)))^2
}
h <- function(ρ, θ, ϕ){
  x <- ρ * cos(θ) * sin(ϕ)
  y <- ρ * sin(θ) * sin(ϕ)
  z <- ρ * cos(ϕ)
  f(x, y, z)
}
# make grid
nρ <- 300L; nθ <- 300L; nϕ <- 300L
ρ <- seq(0, 4.8, length.out = nρ) # ρ runs from 0 to the desired radius
θ <- seq(0, 2*pi, length.out = nθ)
ϕ <- seq(0, pi, length.out = nϕ) 
G <- expand.grid(ρ=ρ, θ=θ, ϕ=ϕ)
# calculate voxel
voxel <- array(with(G, h(ρ, θ, ϕ)), dim = c(nρ, nθ, nϕ))
# calculate isosurface
cont <- contour3d(voxel, level = 0, x = ρ, y = θ, z = ϕ)
# transform to Cartesian coordinates
surf <- t(apply(cont$vertices, 1L, function(ρθϕ){
  ρ <- ρθϕ[1L]; θ <- ρθϕ[2L]; ϕ <- ρθϕ[3L] 
  c(
    ρ * cos(θ) * sin(ϕ),
    ρ * sin(θ) * sin(ϕ),
    ρ * cos(ϕ)
  )
}))

# make mesh
rmesh <- tmesh3d(
  vertices = t(surf),
  indices  = t(cont[["triangles"]]),
  normals  = NULL,
  homogeneous = FALSE
)
rmesh <- Rvcg::vcgClean(rmesh, tol = 5e-7, sel=6)

# plot
open3d(windowRect = c(50, 50, 562, 562), zoom = 0.95)
clear3d(type = "lights")
light3d(x = -50, y = 100, z = 100, ambient = "black")
shade3d(rmesh, color = "#ff00FF", specular = "black")

################################################################################

mesh <- cgalMesh$new(rmesh)
summary(mesh$getEdges())
mesh$isotropicRemeshing(0.02, iterations = 3, relaxSteps = 2)
summary(mesh$getEdges())
mesh
mesh$writeMeshFile("togliatti.off")

M <- cgalMeshes:::testparam(normalizePath("togliatti.off"), 1L)

#### radial checkerboard
M0 <- M
M <- 9.999 * (M0 - 0.5)
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
view3d(0, -70, zoom = 0.8)
shade3d(rmesh, meshColor = "vertices")

# animation ####
movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "Togliatti-DCP.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

