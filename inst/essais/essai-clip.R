library(cgalMeshes)
library(rgl)
mesh <- cgalMesh$new(cube3d())$triangulate()
clipper <- cgalMesh$new(sphereMesh(r= sqrt(2)))

mesh$clipMesh(clipper, TRUE)

rglmesh <- mesh$getMesh()
shade3d(rglmesh, col = "darkorange")