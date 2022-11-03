library(cgalMeshes)
library(Boov)
library(microbenchmark)

rgl1 <- HopfTorusMesh(nu = 30, nv = 30)
rgl2 <- sphereMesh(r = 5, iterations = 4)

mesh1 <- cgalMesh$new(rgl1)
mesh2 <- cgalMesh$new(rgl2)

microbenchmark(
  cgalMesh = mesh1$intersection(mesh2),
  Boov     = MeshesIntersection(list(rgl1, rgl2)),
  times = 5
)

microbenchmark(
  cgalMesh = mesh1$intersection(mesh2)$getMesh(normals = TRUE),
  Boov     = MeshesIntersection(list(rgl1, rgl2), normals = TRUE),
  times = 5
)
