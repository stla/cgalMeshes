library(cgalMeshes)

mesh1 <- cgalMesh$new("mesh_1.ply")
mesh1$selfIntersects() # TRUE
mesh1$removeSelfIntersections()
mesh1$selfIntersects() # FALSE

mesh2 <- cgalMesh$new("mesh_2.ply")
mesh2$selfIntersects() # TRUE
mesh2$removeSelfIntersections()
mesh2$selfIntersects() # FALSE

# union
umesh <- mesh1$union(mesh2)
