library(cgalMeshes)
library(rgl)
# take two cubes
rglmesh1 <- cube3d()
rglmesh2 <- translate3d(cube3d(), 1, 1, 1)
mesh1 <- cgalMesh$new(rglmesh1)
mesh2 <- cgalMesh$new(rglmesh2)
# the two meshes must be triangle
mesh1$triangulate()
mesh2$triangulate()
# intersection
imesh <- mesh1$intersection(mesh2)
rglimesh <- imesh$getMesh(normals = FALSE)
# extract edges for plotting
extEdges <- exteriorEdges(imesh$edges())
# plot
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
shade3d(rglimesh, color = "red")
plotEdges(imesh$vertices(), extEdges)
shade3d(rglmesh1, color = "yellow", alpha = 0.2)
shade3d(rglmesh2, color = "cyan", alpha = 0.2)
