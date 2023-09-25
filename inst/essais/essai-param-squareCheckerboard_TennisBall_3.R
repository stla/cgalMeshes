setwd("C:/SL/MyPackages/cgalMeshes/inst/trash")


library(cgalMeshes)
library(rgl)

# https://laustep.github.io/stlahblog/posts/TennisBall.html ####
f <- function(x, y, z) { # Enneper surface: f=0
  64*z^9 - 128*z^7 + 64*z^5 - 702*x^2*y^2*z^3 - 18*x^2*y^2*z + 
    144*(y^2*z^6-x^2*z^6) + 162*(y^4*z^2-x^4*z^2) + 27*(y^6-x^6) +
    9*(x^4*z+y^4*z) + 48*(x^2*z^3+y^2*z^3) - 432*(x^2*z^5+y^2*z^5) +
    81*(x^4*y^2-x^2*y^4) + 240*(y^2*z^4-x^2*z^4) - 135*(x^4*z^3+y^4*z^3)
}

smesh <- sphereMesh(r = 0.5, iterations = 5L)
rmesh1 <- clipMesh3d(smesh, f, greater = TRUE, minVertices = 20000L)

# convert to CGAL mesh ####
mesh1 <- cgalMesh$new(rmesh1)

# add vertices in order that the checkeboard has regular lines ####
mesh1$isotropicRemeshing(0.001, iterations = 3, relaxSteps = 2)

# compute mesh1 parameterization ####
UV <- mesh1$parameterization(method = "DCP", UVborder = "circle")

# square checkerboard ####
clrs1 <- ifelse(
  (floor(10 * UV[, 1L]) %% 2) == (floor(20 * UV[, 2L]) %% 2), 
  "yellow", "navy"
)

# convert to 'rgl' mesh, and add normals and colors ####
rmesh1 <- mesh1$getMesh()
rmesh1$normals <- rmesh1$vb[-4L, ]
rmesh1$material <- list(color = clrs1)

# make the other part of the tennis ball ####
rmesh2 <- rotate3d(rotate3d(rmesh1, pi, 1, 0, 0), pi/2, 0, 0, 1)
clrs2 <- rmesh2$material$color
clrs2[clrs2 =="navy"] <- "white"
clrs2[clrs2 =="yellow"] <- "black"
rmesh2$material$color <- clrs2

# extract the boundary ####
b <- getBoundary3d(rmesh1, sorted = TRUE, color = "red")

# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
bg3d("#363940")
view3d(0, -20, zoom = 0.7)
shade3d(rmesh1, meshColor = "vertices", polygon_offset = 1)
shade3d(rmesh2, meshColor = "vertices", polygon_offset = 1)
shade3d(b, lwd = 2)

snapshot3d(sprintf("TennisBallWithCheckerboard_BW.png"), webshot = FALSE)


# animation ####
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

# mount animation
library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "TennisBallWithCheckerboard_BW.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))
