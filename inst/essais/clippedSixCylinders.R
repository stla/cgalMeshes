library(rgl)
library(cgalMeshes)

t_ <- seq(0, 1, length.out = 200)

A1 <- c(0, -1, -1)
B1 <- c(0, 1, 1)
v1 <- B1 - A1 
path1 <- t(sapply(t_, function(t) A1 + t*v1))
cyl1 <- as.tmesh3d(cylinder3d(path1, radius = 1, sides = 360, closed = -2))
A2 <- c(0, -1, 1)
B2 <- c(0, 1, -1)
v2 <- B2 - A2 
path2 <- t(sapply(t_, function(t) A2 + t*v2))
cyl2 <- as.tmesh3d(cylinder3d(path2, radius = 1, sides = 360, closed = -2))
A3 <- c(1, 0, 1)
B3 <- c(-1, 0, -1)
v3 <- B3 - A3 
path3 <- t(sapply(t_, function(t) A3 + t*v3))
cyl3 <- as.tmesh3d(cylinder3d(path3, radius = 1, sides = 360, closed = -2))
A4 <- c(1, 0, -1)
B4 <- c(-1, 0, 1)
v4 <- B4 - A4 
path4 <- t(sapply(t_, function(t) A4 + t*v4))
cyl4 <- as.tmesh3d(cylinder3d(path4, radius = 1, sides = 360, closed = -2))
A5 <- c(1, 1, 0)
B5 <- c(-1, -1, 0)
v5 <- B5 - A5 
path5 <- t(sapply(t_, function(t) A5 + t*v5))
cyl5 <- as.tmesh3d(cylinder3d(path5, radius = 1, sides = 360, closed = -2))
A6 <- c(1, -1, 0)
B6 <- c(-1, 1, 0)
v6 <- B6 - A6
path6 <- t(sapply(t_, function(t) A6 + t*v6))
cyl6 <- as.tmesh3d(cylinder3d(path6, radius = 1, sides = 360, closed = -2))

mesh1 <- cgalMesh$new(cyl1)
mesh2 <- cgalMesh$new(cyl2)
mesh3 <- cgalMesh$new(cyl3)
mesh4 <- cgalMesh$new(cyl4)
mesh5 <- cgalMesh$new(cyl5)
mesh6 <- cgalMesh$new(cyl6)
mesh1$assignFaceColors("darkred")
mesh2$assignFaceColors("darkgreen")
mesh3$assignFaceColors("darkblue")
mesh4$assignFaceColors("gold")
mesh5$assignFaceColors("darkorange")
mesh6$assignFaceColors("turquoise")

.meshes1 <- mesh1$clip(mesh2, TRUE)
.meshes2 <- mesh1$clip(mesh3, TRUE)
.meshes3 <- mesh1$clip(mesh4, TRUE)
.meshes4 <- mesh1$clip(mesh5, TRUE)
.meshes5 <- mesh1$clip(mesh6, TRUE)

rmesh1 <- mesh1$getMesh()
shade3d(rmesh1, meshColor = "faces")

gold <- which(mesh1$getFaceColors() == "gold")
mm1 <- mesh1$filterMesh(gold)
mgold <- mm1[[1]]
m2    <- mm1[[2]]

red <- which(m2$getFaceColors() == "darkred")
mm2 <- m2$filterMesh(red)
mred <- mm2[[1]]
m3   <- mm2[[2]]

green <- which(m3$getFaceColors() == "darkgreen")
mm3 <- m3$filterMesh(green)
mgreen <- mm3[[1]]
m4  <- mm3[[2]]

turquoise <- which(m4$getFaceColors() == "turquoise")
mm4 <- m4$filterMesh(turquoise)
mturquoise <- mm4[[1]]
m5  <- mm4[[2]]

blue <- which(m5$getFaceColors() == "darkblue")
mm5 <- m5$filterMesh(blue)
mblue <- mm5[[1]]
morange  <- mm5[[2]]

mgold$computeNormals()
mred$computeNormals()
mgreen$computeNormals()
mturquoise$computeNormals()
mblue$computeNormals()
morange$computeNormals()
rm1 <- mgold$getMesh()
rm2 <- mred$getMesh()
rm3 <- mgreen$getMesh()
rm4 <- mturquoise$getMesh()
rm5 <- mblue$getMesh()
rm6 <- morange$getMesh()

open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(30, 30, zoom = 0.6)
shade3d(rm1, meshColor = "faces")
shade3d(rm2, meshColor = "faces")
shade3d(rm3, meshColor = "faces")
shade3d(rm4, meshColor = "faces")
shade3d(rm5, meshColor = "faces")
shade3d(rm6, meshColor = "faces")

borders <- mgold$getBorders()
vertices <- mgold$getVertices()
for(border in borders) {
  plotEdges(
    vertices, border[, c("v1", "v2")], spheresRadius = 0.02, edgesAsTubes = FALSE
  )
}
borders <- mred$getBorders()
vertices <- mred$getVertices()
for(border in borders) {
  plotEdges(
    vertices, border[, c("v1", "v2")], spheresRadius = 0.02, edgesAsTubes = FALSE
  )
}
borders <- mgreen$getBorders()
vertices <- mgreen$getVertices()
for(border in borders) {
  plotEdges(
    vertices, border[, c("v1", "v2")], spheresRadius = 0.02, edgesAsTubes = FALSE
  )
}
borders <- mturquoise$getBorders()
vertices <- mturquoise$getVertices()
for(border in borders) {
  plotEdges(
    vertices, border[, c("v1", "v2")], spheresRadius = 0.02, edgesAsTubes = FALSE
  )
}
borders <- mblue$getBorders()
vertices <- mblue$getVertices()
for(border in borders) {
  plotEdges(
    vertices, border[, c("v1", "v2")], spheresRadius = 0.02, edgesAsTubes = FALSE
  )
}
borders <- morange$getBorders()
vertices <- morange$getVertices()
for(border in borders) {
  plotEdges(
    vertices, border[, c("v1", "v2")], spheresRadius = 0.02, edgesAsTubes = FALSE
  )
}

# # -- if you want an animation 
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
  movie = "cc",
  convert = FALSE,
  webshot = FALSE
)

library(gifski)
pngs <- list.files(pattern = "^cc.*png$")
gifski(
  pngs, "clippedSixCylinders.gif", width = 512, height = 512, delay = 1/10
)

file.remove(pngs)
