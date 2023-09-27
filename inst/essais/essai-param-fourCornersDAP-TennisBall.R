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
rmesh1 <- clipMesh3d(smesh, f, greater = TRUE, minVertices = 30000L)

# convert to CGAL mesh ####
mesh1 <- cgalMesh$new(rmesh1)

vs <- mesh1$getVertices()
vsx <- vs[, 1L]
head(order(vsx))
head(order(vsx, decreasing = TRUE))

v1 <- 20883L
v2 <- 21397L
v3 <- 21091L
v4 <- 21531L

open3d()
view3d(0, 0)
shade3d(rmesh1)
points3d(rbind(vs[v1, ]), col = "red", size = 12)
points3d(rbind(vs[v2, ]), col = "green", size = 12)
points3d(rbind(vs[v3, ]), col = "blue", size = 12)
points3d(rbind(vs[v4, ]), col = "black", size = 12)

gdists <- mesh1$geoDists(v2)

gdv1 <- gdists[v1]
gdv3 <- gdists[v3]

# add vertices in order that the checkeboard has regular lines ####
mesh1$isotropicRemeshing(0.002, iterations = 3, relaxSteps = 2)

# compute mesh1 parameterization ####
# first, find the four corners
vs <- mesh1$getVertices()
vsx <- vs[, 1L]
head(order(vsx))
head(order(vsx, decreasing = TRUE))

v1 <- 157552L
v2 <- 25446L
v3 <- 157549L
v4 <- 25440L

open3d()
view3d(0, 0)
shade3d(rmesh1)
points3d(rbind(vs[v1, ]), col = "red", size = 12)
points3d(rbind(vs[v2, ]), col = "green", size = 12)
points3d(rbind(vs[v3, ]), col = "blue", size = 12)
points3d(rbind(vs[v4, ]), col = "black", size = 12)

UV <- mesh1$parameterization(method = "DAP", corners = c(v1, v2, v3, v4))

# square checkerboard ####
clrs1 <- ifelse(
  (floor(5 * UV[, 1L]) %% 2) == (floor(5 * UV[, 2L]) %% 2), 
  "yellow", "navy"
)

clrs <- ifelse(
  (floor(5 * (UV[, 1L]+UV[, 2L])) %% 2) == (floor(5 * (UV[, 1L]-UV[, 2L])) %% 2), 
  "yellow", "navy"
)
plot(UV, col = clrs)

clrs <- ifelse(
  (floor(5 * (G$U+G$V)) %% 2) == (floor(5 * (G$U-G$V)) %% 2), 
  "yellow", "navy"
)
plot(G$U, G$V, col = clrs, asp = 1, pch = ".")

# convert to 'rgl' mesh, and add normals and colors ####
rmesh1 <- mesh1$getMesh()
rmesh1$normals <- rmesh1$vb[-4L, ]
rmesh1$material <- list(color = clrs1)

# make the other part of the tennis ball ####
rmesh2 <- rotate3d(rotate3d(rmesh1, pi, 1, 0, 0), pi/2, 0, 0, 1)
clrs2 <- rmesh2$material$color
clrs2[clrs2 =="navy"] <- "springgreen"
clrs2[clrs2 =="yellow"] <- "firebrick4"
rmesh2$material$color <- clrs2

# extract the boundary ####
b <- getBoundary3d(rmesh1, sorted = TRUE, color = "black")

# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
bg3d("#363940")
view3d(0, -20, zoom = 0.7)
shade3d(rmesh1, meshColor = "vertices", polygon_offset = 1)
shade3d(rmesh2, meshColor = "vertices", polygon_offset = 1)
shade3d(b, lwd = 2)

snapshot3d(sprintf("TennisBallWithCheckerboard_DAP_FourCorners.png"), webshot = FALSE)


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
  gif_file = "TennisBallWithCheckerboard_DAP_FourCorners.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))


