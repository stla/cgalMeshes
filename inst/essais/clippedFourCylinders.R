library(rgl)
library(cgalMeshes)

t_ <- seq(0, 1, length.out = 100)

A1 <- c(2, 2, 2)
B1 <- c(-2, -2, -2)
v1 <- B1 - A1 
path1 <- t(sapply(t_, function(t) A1 + t*v1))
cyl1 <- as.tmesh3d(cylinder3d(path1, radius = 1, sides = 120, closed = -2))
A2 <- c(2, 2, -2)
B2 <- c(-2, -2, 2)
v2 <- B2 - A2 
path2 <- t(sapply(t_, function(t) A2 + t*v2))
cyl2 <- as.tmesh3d(cylinder3d(path2, radius = 1, sides = 120, closed = -2))
A3 <- c(2, -2, 2)
B3 <- c(-2, 2, -2)
v3 <- B3 - A3 
path3 <- t(sapply(t_, function(t) A3 + t*v3))
cyl3 <- as.tmesh3d(cylinder3d(path3, radius = 1, sides = 120, closed = -2))
A4 <- c(2, -2, -2)
B4 <- c(-2, 2, 2)
v4 <- B4 - A4 
path4 <- t(sapply(t_, function(t) A4 + t*v4))
cyl4 <- as.tmesh3d(cylinder3d(path4, radius = 1, sides = 120, closed = -2))

mesh1 <- cgalMesh$new(cyl1)
mesh2 <- cgalMesh$new(cyl2)
mesh3 <- cgalMesh$new(cyl3)
mesh4 <- cgalMesh$new(cyl4)
mesh1$assignFaceColors("darkred")
mesh2$assignFaceColors("darkgreen")
mesh3$assignFaceColors("darkblue")
mesh4$assignFaceColors("gold")

.meshes1 <- mesh1$clip(mesh2, TRUE)
.meshes2 <- mesh1$clip(mesh3, TRUE)
.meshes3 <- mesh1$clip(mesh4, TRUE)

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
mblue  <- mm3[[2]]

mgold$computeNormals()
mred$computeNormals()
mgreen$computeNormals()
mblue$computeNormals()
rm1 <- mgold$getMesh()
rm2 <- mred$getMesh()
rm3 <- mgreen$getMesh()
rm4 <- mblue$getMesh()

open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(30, 30, zoom = 0.6)
shade3d(rm1, meshColor = "faces")
shade3d(rm2, meshColor = "faces")
shade3d(rm3, meshColor = "faces")
shade3d(rm4, meshColor = "faces")

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
  webshot = TRUE
)

library(gifski)
pngs <- list.files(pattern = "^cc.*png$")
gifski(
  pngs, "clippedFourCylinders.gif", width = 512, height = 512, delay = 1/10
)

file.remove(pngs)
