library(cgalMeshes)
library(rgl)
# the dodecahedron is the dual of the icosahedron
icosahedron <- cgalMesh$new(icosahedron3d())
dodecahedron <- icosahedron$dual()
# take vertices and edges for plotting
vertices <- dodecahedron$vertices()
edges <- dodecahedron$edges()
# plot
rglMesh <- dodecahedron$triangulate()$getMesh(normals = FALSE)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(30, 40, zoom = 0.7)
shade3d(rglMesh, color = "navy")
plotEdges(
  vertices, edges, color = "yellow"
)

