library(cgalMeshes)
library(rgl)

rglMesh <- HopfTorusMesh(nu = 50, nv = 50)
mesh <- cgalMesh$new(rglMesh)
mesh$smoothAngle(iterations = 3, safety = TRUE)
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

