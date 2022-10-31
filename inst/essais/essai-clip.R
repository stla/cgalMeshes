library(cgalMeshes)
library(rgl)
mesh <- cgalMesh$new(cube3d())$triangulate()
clipper <- cgalMesh$new(sphereMesh(r= sqrt(2)))

mesh$clipMesh(clipper, TRUE)
