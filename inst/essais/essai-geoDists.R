library(cgalMeshes)
library(rgl)

data(bunny, package = "onion")
rglmesh <- SurfaceReconstruction::AFSreconstruction(bunny)




library(cgalMeshes)
library(rgl)
rglmesh <- torusMesh(R = 3, r = 2, nu = 90, nv = 60)
mesh <- cgalMesh$new(rglmesh)
# estimated geodesic distances
geodists <- mesh$geoDists(1L)
# normalization to (0, 1)
geodists <- (geodists - min(geodists)) / (max(geodists) - min(geodists))
# color each vertex according to its geodesic distance from the source
fcolor <- colorRamp(viridisLite::turbo(200L))
colors <- fcolor(geodists)
colors <- rgb(colors[, 1L], colors[, 2L], colors[, 3L], maxColorValue = 255)
rglmesh[["material"]] <- list("color" = colors)
# plot
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
shade3d(rglmesh)
wire3d(rglmesh, color = "black")
if(!rgl.useNULL()) {
  play3d(spin3d(axis = c(1, 1, 1), rpm = 5), duration = 20)  
}


# a trefoil knot (taken from `?rgl::cylinder3d`) ####
library(cgalMeshes)
library(rgl)
theta <- seq(0, 2*pi, length.out = 50L)
knot <- cylinder3d(
  center = cbind(
    sin(theta) + 2*sin(2*theta), 
    2*sin(3*theta), 
    cos(theta) - 2*cos(2*theta)),
  e1 = cbind(
    cos(theta) + 4*cos(2*theta), 
    6*cos(3*theta), 
    sin(theta) + 4*sin(2*theta)),
  radius = 0.8, 
  closed = TRUE)
knot <- subdivision3d(knot, depth = 2)
mesh <- cgalMesh$new(knot)$triangulate()
rglmesh <- mesh$getMesh()
# estimated geodesic distances
geodists <- mesh$geoDists(1L)
# normalization to (0, 1)
geodists <- (geodists - min(geodists)) / (max(geodists) - min(geodists))
# color each vertex according to its geodesic distance from the source
fcolor <- colorRamp(viridisLite::inferno(200L))
colors <- fcolor(geodists)
colors <- rgb(colors[, 1L], colors[, 2L], colors[, 3L], maxColorValue = 255)
rglmesh[["material"]] <- list("color" = colors)
# plot
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
shade3d(rglmesh)

movie3d(spin3d(axis = c(1, 1, 0), rpm = 10),
        duration = 6, fps = 20,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/20)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "trefoilKnot.gif",
  width = 512,
  height = 512,
  delay = 1/11
)




if(!rgl.useNULL()) {
  play3d(spin3d(axis = c(1, 1, 0), rpm = 5), duration = 20)  
}

