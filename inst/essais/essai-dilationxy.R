library(cgalMeshes)
library(rgl)

n <- 30
vs <- t(sapply(1:n, function(i) {
  c(cospi(2*i/n), sinpi(2*i/n), 0)
}))
vs <- rbind(c(0, 0, 0), vs)

faces <- t(sapply(2:n, function(i) {
  c(1, i, i + 1)
}))
faces <- rbind(faces, c(1, n+1, 2))

circle <- cgalMesh$new(vertices = vs, faces = faces, clean = TRUE)
torus <- cgalMesh$new(torusMesh(30, 25, nu = 30, nv = 20))
torus <- cgalMesh$new(octahedron3d())
torus <- cgalMesh$new(sphereMesh(it=2))

msum <- MinkowskiSum(circle, torus)
msum$triangulate()
rmesh <- msum$getMesh()

shade3d(rmesh, color = "green")

