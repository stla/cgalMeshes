library(rgl)
library(cgalMeshes)

t_ <- seq(0, 1, length.out = 50)

A1 <- c(-2, 0, 0)
B1 <- c(2, 0, 0)
v1 <- B1 - A1 
path1 <- t(sapply(t_, function(t) A1 + t*v1))
cyl1 <- as.tmesh3d(cylinder3d(path1, radius = 1, sides = 60, closed = -2))
A2 <- c(0, -2, 0)
B2 <- c(0, 2, 0)
v2 <- B2 - A2 
path2 <- t(sapply(t_, function(t) A2 + t*v2))
cyl2 <- as.tmesh3d(cylinder3d(path2, radius = 1, sides = 60, closed = -2))
A3 <- c(0, 0, -2)
B3 <- c(0, 0, 2)
v3 <- B3 - A3 
path3 <- t(sapply(t_, function(t) A3 + t*v3))
cyl3 <- as.tmesh3d(cylinder3d(path3, radius = 1, sides = 60, closed = -2))

mesh1 <- cgalMesh$new(cyl1)
mesh2 <- cgalMesh$new(cyl2)
mesh3 <- cgalMesh$new(cyl3)
mesh1$assignFaceColors("red")
mesh2$assignFaceColors("green")
mesh3$assignFaceColors("blue")

.meshes = mesh1$clip(mesh2, TRUE)

meshes = mesh1$clip(mesh3, TRUE)
m1 = meshes[[1]]
m2 = meshes[[2]]

red <- which(m1$getFaceColors() == "red")
mm <- m1$filterMesh(red)
m1a <- mm[[1]]
m1b <- mm[[2]]

m1a$computeNormals()
m1b$computeNormals()
m2$computeNormals()
rm1 <- m1a$getMesh()
rm2 <- m1b$getMesh()
rm3 <- m2$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512)); view3d(30, 30)
shade3d(rm1, color = "red")
shade3d(rm2, color = "green")
shade3d(rm3, color = "blue")
