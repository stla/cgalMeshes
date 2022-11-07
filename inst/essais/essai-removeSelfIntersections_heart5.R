library(cgalMeshes)
library(rgl)

ply <- "inst/essais/HEART_OBS5.ply"

mesh <- cgalMesh$new(ply)
rglMesh <- mesh$getMesh(normals = FALSE)

mesh$selfIntersects()
mesh$removeSelfIntersections()
mesh$isValid()
rglMesh <- mesh$getMesh(normals = FALSE)

