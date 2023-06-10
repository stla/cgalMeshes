library(cgalMeshes)
library(rgl)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash/vignette")

Enneper <- function(phi, r, n = 4) {
  rbind(
    r*cos(phi) - r^(2*n-1) * cos((2*n-1)*phi) / (2*n-1),
    r*sin(phi) + r^(2*n-1) * sin((2*n-1)*phi) / (2*n-1),
    2*r^n * cos(n*phi)/n
  )
}

rmesh <- parametricMesh(
  Enneper, urange = c(0, 2*pi), vrange = c(0, 1.1),
  periodic = c(TRUE, FALSE), nu = 512, nv = 512, clean = FALSE
)
rmesh <- Rvcg::vcgClean(rmesh, sel = 0)

open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-10, -35, zoom = 0.7)
shade3d(rmesh, color = "maroon", polygon_offset = 1)
shade3d(
  getBoundary3d(
    rmesh, sorted = TRUE, color = "black", 
    lwd = 3, line_antialias = TRUE
  )
)
contourLines3d(
  rmesh, fn = function(x, y, z) sqrt(x^2 + y^2), 
  levels = seq(0.1, 1.1, by = 0.2), 
  plot = TRUE, lwd = 2, line_antialias = TRUE
)

# animation ####
movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "Enneper-maroon.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))


# checkerboard ####
u <- v <- seq(0, 1, length.out = 513L)[-513L]
Grid <- expand.grid(U = u, V = v)
checkerboard <- ifelse(
  (floor(10*Grid$U) %% 2) == (floor(10*Grid$V) %% 2), 
  "yellow", "navy"
)
opar <- par(mar = c(2, 2, 1, 1))
plot(
  Grid$U, Grid$V, type = "p", asp = 1, pch = ".", 
  col = checkerboard, xlab = "u", ylab = "v", axes = FALSE
)
axis(1); axis(2)
par(opar)

# assign these colors to the Enneper mesh
rmesh$material <- list(color = checkerboard)

open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-10, -35, zoom = 0.7)
shade3d(rmesh, meshColor = "vertices")

snapshot3d(
  "Enneper-checkerboard.png", width = 512, height = 512, webshot = TRUE
)

# convert the rgl mesh to cgalMesh
mesh <- cgalMesh$new(rmesh)

# look at the edge lengths
summary(mesh$getEdges()[["length"]])

# do an isotropic remeshing to get smaller faces
mesh$isotropicRemeshing(5e-3, iterations = 3L, relaxSteps = 2L)
summary(mesh$getEdges()[["length"]])
mesh$writeMeshFile("Enneper-remeshed.off")

# compute the discrete conformal parameterization
UV <- mesh$parameterization("DCP", UVborder = "circle")
head(UV)

# compute the ARAP parameterization
UV <- mesh$parameterization("ARAP", lambda = 1000, UVborder = "square")
head(UV)
plot(UV, type = "p", asp = 1, pch = ".")
UV[, 1L] <- UV[, 1L] - min(UV[, 1L])
UV[, 2L] <- UV[, 2L] - min(UV[, 2L])
UV[, 1L] <- UV[, 1L] / max(UV[, 1L])
UV[, 2L] <- UV[, 2L] / max(UV[, 2L])

# compute the iterative authalic parameterization
UV <- mesh$parameterization("IAP", UVborder = "square", iterations = 1)
head(UV)

# compute the discrete authalic parameterization
UV <- mesh$parameterization("DAP", UVborder = "circle")
head(UV)

# make a checkerboard with these points
UVcheckerboard <- ifelse(
  (floor(10*UV[, 1L]) %% 2) == (floor(10*UV[, 2L]) %% 2), 
  "yellow", "navy"
)

# compute mesh normals and convert to rgl mesh
mesh$computeNormals()
rmesh <- mesh$getMesh()

# assign these colors to the mesh and plot
rmesh$material <- list(color = UVcheckerboard)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-10, -35, zoom = 0.7)
shade3d(rmesh, meshColor = "vertices")

snapshot3d(
  "Enneper-DAP-circleBorder.png", width = 512, height = 512, webshot = FALSE
)


##| mapping an image ####
mesh <- cgalMesh$new("Enneper-remeshed.off")
# compute the discrete conformal parameterization with square border
UV <- mesh$parameterization("DCP", UVborder = "square")

# now we extract the colors from the image
library(imager)
# load the image
img <- load.image("spiral512x512.png")
# take the r, g, b channels
r <- squeeze(R(img))
g <- squeeze(G(img))
b <- squeeze(B(img))

# make interpolation functions to get the colors of the UV points
library(cooltools) # to get the `approxfun2` function
x_ <- y_ <- seq(0, 1, length.out = 512L)
f_r <- approxfun2(x_, y_, r)
f_g <- approxfun2(x_, y_, g)
f_b <- approxfun2(x_, y_, b)

# now, interpolate the r, g, b values
UV_r <- f_r(UV[, 1L], UV[, 2L])
UV_g <- f_g(UV[, 1L], UV[, 2L])
UV_b <- f_b(UV[, 1L], UV[, 2L])

# convert rgb to hex codes
clrs <- rgb(UV_r, UV_g, UV_b)

# we're done; now compute mesh normals, convert to rgl mesh, and assign colors
mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

# plot
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-10, -35, zoom = 0.7)
shade3d(rmesh, meshColor = "vertices")
snapshot3d(
  "Enneper-spiral.png", width = 512, height = 512, webshot = FALSE
)


##| yin yang ####
sphericalCoordinates <- function(θ, ϕ){
  x <- cos(θ) * sin(ϕ)
  y <- sin(θ) * sin(ϕ)
  z <- cos(ϕ)
  rbind(x, y, z)
}

rmesh <- parametricMesh(
  sphericalCoordinates, urange = c(0, 2*pi), vrange = c(0, pi/4),
  periodic = c(TRUE, FALSE), nu = 512, nv = 512, clean = FALSE
)
rmesh <- Rvcg::vcgClean(rmesh, sel = 0)

mesh <- cgalMesh$new(rmesh)

# look at the edge lengths
summary(mesh$getEdges()[["length"]])

# do an isotropic remeshing to get smaller faces
mesh$isotropicRemeshing(5e-3, iterations = 3L, relaxSteps = 2L)
summary(mesh$getEdges()[["length"]])

# compute the discrete conformal parameterization
UV <- mesh$parameterization("DCP", UVborder = "circle")

yy <- function(uv) {
  u <- uv[1L]
  v <- uv[2L]
  if(v > 0.5) {
    if(u < 0.5) {
      "yellow"
    } else {
      ifelse((u-0.5)^2 + (v-0.75)^2 < 0.25^2, "yellow", "navy")
    }
  } else {
    if(u < 0.5) {
      ifelse((u-0.5)^2 + (v-0.25)^2 < 0.25^2, "navy", "yellow")
    } else {
      "navy"
    }
  }
}

clrs <- apply(UV, 1L, yy)
rmesh$material <- list(color = clrs)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-10, -15, zoom = 0.7)
shade3d(rmesh, color = "darkorange", polygon_offset = 1)
contourLines3d(
  rmesh, fn = function(x, y, z) sqrt(x^2 + y^2), 
  levels = seq(0.1, 0.6, by = 0.1), 
  plot = TRUE, lwd = 2, line_antialias = TRUE
)
contourLines3d(
  rmesh, fn = function(x, y, z) atan2(y, x), 
  levels = seq(-pi, by = pi/4, length.out = 8L), 
  plot = TRUE, lwd = 2, line_antialias = TRUE
)
snapshot3d("sphericalCap-orange.png", webshot = TRUE)

open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-10, -15, zoom = 0.7)
shade3d(rmesh, meshColor = "vertices", polygon_offset = 1)
contourLines3d(
  rmesh, fn = function(x, y, z) sqrt(x^2 + y^2), 
  levels = seq(0.1, 0.6, by = 0.1), 
  plot = TRUE, lwd = 2, line_antialias = TRUE
)
contourLines3d(
  rmesh, fn = function(x, y, z) atan2(y, x), 
  levels = seq(-pi, by = pi/4, length.out = 8L), 
  plot = TRUE, lwd = 2, line_antialias = TRUE
)
snapshot3d("sphericalCap-yinyang.png", webshot = TRUE)
