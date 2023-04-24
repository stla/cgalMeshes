library(rgl)
library(cgalMeshes)

t_ <- seq(0, 1, length.out = 150)

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

meshes1 = mesh1$clip(mesh2, TRUE)
cols <- mesh1$getFaceColors()
cols[cols==""] <- "darkred"
mesh1$assignFaceColors(cols)

# m1 <- meshes1[[1]]
# cols <- m1$getFaceColors()
# cols[cols==""] <- "white"
# m1$assignFaceColors(cols)
# 
# m2 <- meshes1[[2]]
# mmm = m1$union(m2)

meshes2 = mesh3$clip(mesh4, TRUE)
cols <- mesh3$getFaceColors()
cols[cols==""] <- "darkblue"
mesh3$assignFaceColors(cols)

mmesh1 <- mesh1$copy()
mmesh3 <- mesh3$copy()

# cols <- mmesh1$getFaceColors()
# rmesh1 <- Rvcg::vcgClean(mmesh1$getMesh(), sel = 1:7)
# ok <- rmesh1$remface == 0
# rmesh1$material$color <- cols[ok]
# 
# cols <- mmesh3$getFaceColors()
# rmesh3 <- Rvcg::vcgClean(mmesh3$getMesh(), sel = 0:7)
# ok <- rmesh3$remface == 0
# rmesh3$material$color <- cols[ok]
# 
# mmmesh1 <- cgalMesh$new(rmesh1)
# mmmesh3 <- cgalMesh$new(rmesh3)

. <- mesh3$clip(mesh1, FALSE)
# table(mesh3$getFaceColors())
# mesh3$computeNormals()
# rmesh3 <- mesh3$getMesh()
# shade3d(rmesh3, meshColor = "faces")

. <- mmesh1$clip(mmesh3, FALSE)
# table(mmesh1$getFaceColors())
# mmesh1$computeNormals()
# rmesh1 <- mmesh1$getMesh()
# shade3d(rmesh1, meshColor = "faces")




m1 = mesh3
m2 = mmesh1

gold <- which(m1$getFaceColors() == "gold")
mm1 <- m1$filterMesh(gold)
m1a <- mm1[[1]]
m1b <- mm1[[2]]

red <- which(m2$getFaceColors() == "darkred")
mm2 <- m2$filterMesh(red)
m2a <- mm2[[1]]
m2b <- mm2[[2]]

m1a$computeNormals()
m1b$computeNormals()
m2a$computeNormals()
m2b$computeNormals()
rm1 <- m1a$getMesh()
rm2 <- m1b$getMesh()
rm3 <- m2a$getMesh()
rm4 <- m2b$getMesh()
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
  webshot = FALSE
)

library(gifski)
pngs <- list.files(pattern = "^cc.*png$")
gifski(
  pngs, "clippedFourCylinders.gif", width = 512, height = 512, delay = 1/10
)

file.remove(pngs)
