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
rglmesh0 <- rglmesh
vs <- mesh$vertices()
# estimated geodesic distances
index <- which.min(vs[, 1L])
geodists <- mesh$geoDists(index)
x <- seq(min(geodists), max(geodists), length.out = 90L)
# coloring function
fcolor <- colorRamp(viridisLite::inferno(200L))
cols <- fcolor(geodists / max(geodists))
colors0 <- rgb(cols[, 1L], cols[, 2L], cols[, 3L], maxColorValue = 255)

for(i in seq_along(x)) {
  rglmesh <- rglmesh0
  colors <- colors0
  red <- geodists >= x[i]
  colors[red] <- "white"
  rglmesh[["vb"]][, which(red)] <- NA_real_
  # colors <- rep("white", nrow(vs))
  # red <- geodists <= x[i]
  # cols <- fcolor(geodists[red] / max(geodists[red]))
  # colors[which(red)] <- rgb(cols[, 1L], cols[, 2L], cols[, 3L], maxColorValue = 255)
  rglmesh[["material"]] <- list("color" = colors)
  # plot
  open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
  shade3d(rglmesh0, alpha = 0)
  shade3d(rglmesh)
  snapshot3d(sprintf("zzpic%03d.png", i), webshot = FALSE)
  close3d()
}

command <- "convert -delay 1x11 -duplicate 1,-2-1 -layers OptimizePlus zzpic*.png trefoilKnot4.gif"
system(command)
file.remove(Sys.glob("zzpic*.png"))
