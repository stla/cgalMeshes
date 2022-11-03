library(cgalMeshes)
library(Boov)
library(MeshesTools)
library(microbenchmark)

rgl1 <- HopfTorusMesh(nu = 100, nv = 100)
rgl2 <- sphereMesh(r = 5, iterations = 4)

f1 <- function() {
  mesh1 <- cgalMesh$new(rgl1)
  mesh2 <- cgalMesh$new(rgl2)
  inter <- mesh1$intersection(mesh2)
  vol <- inter$volume()
  area <- inter$area()
  # centroid <- inter$centroid()
  rglinter <- inter$getMesh()
}

# open3d(windowRect = c(50, 50, 562, 562))
# view3d(0, 0, zoom = 0.6)
# shade3d(rglinter, color = "yellow")
# wire3d(rglinter)

f2 <- function() {
  inter <- MeshesIntersection(list(rgl1, rgl2))
  vol <- meshVolume(inter)
  area <- meshArea(inter)
  # centroid <- meshCentroid(inter)
  rglinter <- toRGL(inter)
}

o <- capture.output(
  microbenchmark(
    cgalMesh = f1(),
    others   = f2(),
    times = 3
  )
)



mesh1 <- cgalMesh$new(rgl1)
mesh2 <- cgalMesh$new(rgl2)

capture.output(
  microbenchmark(
    cgalMesh = mesh1$intersection(mesh2),
    Boov     = MeshesIntersection(list(rgl1, rgl2)),
    times = 5
  )
)

microbenchmark(
  cgalMesh = mesh1$intersection(mesh2)$getMesh(normals = TRUE),
  Boov     = MeshesIntersection(list(rgl1, rgl2), normals = TRUE),
  times = 5
)
