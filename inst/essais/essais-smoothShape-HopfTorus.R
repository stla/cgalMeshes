library(cgalMeshes)
library(rgl)

rglMesh <- HopfTorusMesh(nu = 30, nv = 20)
mesh <- cgalMesh$new(rglMesh)

areas <- mesh$getFacesInfo()[, "area"]
summary(areas)
bigFaces <- which(areas > 0.00088)
length(bigFaces)

summary(mesh$getFacesInfo())
cy <- mesh$getFacesInfo()[, "ccx"]
bigFaces <- which(cy > 0)
length(bigFaces)
#bigFaces <- 4801:9600
bigFaces <- NULL
mesh$smoothShape(bigFaces, time = 2, iterations = 2)
mesh$computeNormals()
rglSmoothMesh <- mesh$getMesh()

open3d(windowRect = 50 + c(0, 0, 800, 400))
mfrow3d(1, 2)
view3d(0, 0, zoom = 0.85)
shade3d(rglMesh, color = "orangered")
wire3d(rglMesh)
next3d()
view3d(0, 0, zoom = 0.85)
shade3d(rglSmoothMesh, color = "orangered")
wire3d(rglSmoothMesh)

