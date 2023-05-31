library(cgalMeshes)
library(rgl)

rglMesh <- HopfTorusMesh(nu = 100, nv = 80)
rglMesh <- cyclideMesh(a = 97, c = 32, mu = 57, nu = 30L, nv = 20L)
mesh <- cgalMesh$new(rglMesh)

areas <- mesh$getFacesInfo()[, "area"]
bigFaces <- which(areas > 0.6)
length(bigFaces)

mesh$smoothAngle(bigFaces, iterations = 100, safety = FALSE)
mesh$computeNormals()

rglSmoothMesh <- mesh$getMesh()

open3d(windowRect = 50 + c(0, 0, 800, 400))
mfrow3d(1, 2)
view3d(0, 0, zoom = 0.85)
shade3d(rglMesh, color = "orange")
wire3d(rglMesh)
next3d()
view3d(0, 0, zoom = 0.85)
shade3d(rglSmoothMesh, color = "orange")
wire3d(rglSmoothMesh)

