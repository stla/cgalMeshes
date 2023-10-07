setwd("C:/SL/MyPackages/cgalMeshes/inst/trash")

library(cgalMeshes)
library(rgl)

n <- 2
Enneper <- function(phi, r) {
  rbind(
    r*cos(phi) - r^(2*n-1)*cos((2*n-1)*phi)/(2*n-1),
    r*sin(phi) + r^(2*n-1)*sin((2*n-1)*phi)/(2*n-1),
    2*r^n*cos(n*phi)/n
  )
}

rmesh <- parametricMesh(
  Enneper, urange = c(0, 2*pi), vrange = c(0, 1.2),
  periodic = c(TRUE, FALSE), nu = 512L, nv = 128L, clean = TRUE
)

shade3d(rmesh, col = "deeppink4")
bbox3d()

# convert to CGAL mesh ####
mesh <- cgalMesh$new(rmesh)

# take a look at the edge lengths
edges <- mesh$getEdges()
summary(edges[["length"]])

# add vertices in order that the checkeboard has regular lines ####
mesh$isotropicRemeshing(0.008, iterations = 3, relaxSteps = 2)

# compute mesh parameterization ####
# first, find the four corners
vs <- mesh$getVertices()
vsx <- vs[, 1L]
oinc <- order(vsx)
odec <- order(vsx, decreasing = TRUE)
v1 <- oinc[1L]
v2 <- oinc[2L]
v3 <- odec[2L]
v4 <- odec[1L]

open3d()
view3d(0, 0)
shade3d(rmesh)
points3d(rbind(vs[v1, ]), col = "red", size = 12)
points3d(rbind(vs[v2, ]), col = "green", size = 12)
points3d(rbind(vs[v3, ]), col = "blue", size = 12)
points3d(rbind(vs[v4, ]), col = "black", size = 12)

UV <- mesh$parameterization(method = "DCP", corners = c(v1, v2, v3, v4))

# square checkerboard ####
checkerboard <- ifelse(
  (floor(5 * UV[, 1L]) %% 2) == (floor(5 * UV[, 2L]) %% 2), 
  "yellow", "navy"
)

# add normals, convert to 'rgl' mesh, and add colors ####
mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh[["material"]] <- list("color" = checkerboard)
b <- getBoundary3d(rmesh, sorted = TRUE, color = "black")

# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
bg3d("#363940")
par3d(userMatrix = um)
shade3d(rmesh, meshColor = "vertices", polygon_offset = 1)
shade3d(b, lwd = 4)

snapshot3d(sprintf("EnneperWithCheckerboard_FourCorners.png"), webshot = FALSE)
saveRDS(um, "userMatrix.rds")

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
  gif_file = "TennisBallWithCheckerboard_FourCorners.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))


