library(cgalMeshes)
library(rgl)

rglMesh <- HopfTorusMesh(nu = 80, nv = 60)
mesh <- cgalMesh$new(rglMesh)

areas <- mesh$getFacesInfo()[, "area"]
summary(areas)
bigFaces <- which(areas > 0.00088)
length(bigFaces)

# cy <- mesh$getFacesInfo()[, "cy"]
# bigFaces <- which(cy > -5.8)
# length(bigFaces)

mesh$smoothShape(bigFaces, time = 1, iterations = 4)
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

