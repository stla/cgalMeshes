library(cgalMeshes)

mesh1 <- cgalMesh$new("inst/essais/mesh_1.ply")
mesh1$selfIntersects() # TRUE
mesh1$removeSelfIntersections()
mesh1$selfIntersects() # FALSE

mesh2 <- cgalMesh$new("inst/essais/mesh_2.ply")
mesh2$selfIntersects() # TRUE
mesh2$removeSelfIntersections()
mesh2$selfIntersects() # FALSE

# union
umesh <- mesh1$union(mesh2)
rglumesh <- umesh$getMesh(normals = FALSE)

rgl::shade3d(rglumesh, color = "orange")

