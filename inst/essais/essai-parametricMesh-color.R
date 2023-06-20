##| example with colors: torus figure 8 ####
library(rgl)
library(cgalMeshes)
library(colorspace) # to use the  lighten  function
# function mapping (u,v) to a color
clrs <- c(
  hcl.colors(128L, "Rocket"), 
  hcl.colors(128L, "Rocket", rev = TRUE)
)
framp <- colorRamp(clrs)
fcolor <- function(u, v) {
  cols <- framp(v/(2*pi))
  cols <- rgb(cols[, 1L], cols[, 2L], cols[, 3L], maxColorValue = 255)
  lighten(cols, 0.3*cos(u))
}
# parameterization
f <- function(u, v, c = 1){
  h <- c + sin(v) * cos(u) - sin(2*v) * sin(u) / 2
  x <- h * cos(u)
  y <- h * sin(u)
  z <- sin(u) * sin(v) + cos(u) * sin(2*v) / 2
  rbind(x, y, z)  
}
# make the mesh and plot it
rmesh <- parametricMesh(
  f, c(0, 2*pi), c(0, 2*pi), periodic = c(TRUE, TRUE),
  nu = 100L, nv = 100L, fcolor = fcolor)
rmesh <- addNormals(rmesh)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
shade3d(rmesh, meshColor = "vertices")

##| animation ####
M <- par3d("userMatrix")
movie3d(
  par3dinterp(
    time = seq(0, 1, len = 9),
    userMatrix = list(
      M,
      rotate3d(M, pi, 1, 0, 0),
      rotate3d(M, pi, 1, 1, 0),
      rotate3d(M, pi, 1, 1, 1),
      rotate3d(M, pi, 0, 1, 1),
      rotate3d(M, pi, 0, 1, 0),
      rotate3d(M, pi, 1, 0, 1),
      rotate3d(M, pi, 0, 0, 1),
      M
    )
  ),
  fps = 100,
  duration = 1,
  dir = ".",
  movie = "zzpic",
  convert = FALSE, webshot = FALSE
)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "figure8.gif",
  width = 512,
  height = 512,
  delay = 1/10
)
file.remove(Sys.glob("zzpic*.png"))