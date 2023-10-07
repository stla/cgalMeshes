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

bndry <- getBoundary3d(rmesh, sorted = TRUE, color = "black")

open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
par3d(userMatrix = um)
shade3d(rmesh, col = "deeppink4")
shade3d(bndry, lwd = 4)
snapshot3d("EnneperOrderTwo.png", webshot = FALSE)

# convert to CGAL mesh ####
mesh <- cgalMesh$new(rmesh)

# take a look at the edge lengths
edges <- mesh$getEdges()
summary(edges[["length"]])

# add vertices in order that the checkeboard has regular lines ####
mesh$isotropicRemeshing(0.008, iterations = 3, relaxSteps = 2)

# compute mesh parameterization ####
# without the four corners
UV <- mesh$parameterization(method = "DCP")
# the UV-space is the square [0,1]x[0,1]
# make square checkerboard with 5 squares x 5 squares ####
checkerboard <- ifelse(
  (floor(5 * UV[, 1L]) %% 2) == (floor(5 * UV[, 2L]) %% 2), 
  "gold", "magenta4"
)
# add normals, convert to 'rgl' mesh, and add colors ####
mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh[["material"]] <- list("color" = checkerboard)
# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
bg3d("#363940")
par3d(userMatrix = um)
shade3d(rmesh, meshColor = "vertices", polygon_offset = 1)
shade3d(b, lwd = 4)
snapshot3d(sprintf("EnneperCheckerboard_NoCorners.png"), webshot = FALSE)

# now we will fix four corners ####
# first, find the four corners
vs <- mesh$getVertices()
vsx <- vs[, 1L]
oinc <- order(vsx)
odec <- order(vsx, decreasing = TRUE)
v1 <- oinc[1L]
v2 <- oinc[2L]
v3 <- odec[2L]
v4 <- odec[1L]
# plot with the four corner vertices
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, -25, zoom = 0.8)
shade3d(rmesh, color = "deeppink4", polygon_offset = 1)
shade3d(bndry, lwd = 3)
points3d(rbind(vs[v1, ]), col = "red",    size = 12)
points3d(rbind(vs[v2, ]), col = "green",  size = 12)
points3d(rbind(vs[v3, ]), col = "blue",   size = 12)
points3d(rbind(vs[v4, ]), col = "yellow", size = 12)
snapshot3d("Enneper_fourCorners.png", webshot = FALSE)

# compute the DCP parameterization with the four given corners
UV <- mesh$parameterization(method = "DCP", corners = c(v1, v2, v3, v4))

# the UV-space is the square [0,1]x[0,1]
# make square checkerboard with 5 squares x 5 squares ####
checkerboard <- ifelse(
  (floor(5 * UV[, 1L]) %% 2) == (floor(5 * UV[, 2L]) %% 2), 
  "gold", "magenta4"
)

# add normals, convert to 'rgl' mesh, and add colors ####
mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh[["material"]] <- list("color" = checkerboard)
#b <- getBoundary3d(rmesh, sorted = TRUE, color = "black")

# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
bg3d("#363940")
par3d(userMatrix = um)
shade3d(rmesh, meshColor = "vertices", polygon_offset = 1)
shade3d(b, lwd = 4)
snapshot3d(sprintf("EnneperCheckerboard_FourCorners.png"), webshot = FALSE)


# animation ####
movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 20,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/20)

# mount animation
library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "EnneperCheckerboard_FourCorners.gif",
  width = 512,
  height = 512,
  delay = 1/10
)
file.remove(Sys.glob("zzpic*.png"))


