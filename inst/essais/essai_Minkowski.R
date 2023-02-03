library(cgalMeshes)
library(rgl)
# cube + icosahedron
mesh1 <- cgalMesh$new(cube3d())$triangulate()
mesh2 <- cgalMesh$new(icosahedron3d())
mesh <- MinkowskiSum(mesh1, mesh2)
# get the edges before triangulation for plotting
edges <- mesh$getEdges()
# triangulation
mesh$triangulate()
# plot
rmesh <- mesh$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(20, 40, zoom = 0.85)
shade3d(rmesh, color = "turquoise")
plotEdges(
  mesh$getVertices(), edges[, c("i1", "i2")], 
  color = "darkred", tubesRadius = 0.07, spheresRadius = 0.1
)
