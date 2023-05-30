library(cgalMeshes)
library(sphereTessellation)
library(rgl)

icosphere <- cgalMesh$new(icosphereMesh())
dicosphere <- icosphere$dual()
# take vertices and edges for plotting
vertices <- dicosphere$getVertices()
edges <- dicosphere$getEdges()
# plot
rglMesh <- dicosphere$triangulate()$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(30, 40, zoom = 0.7)
shade3d(rglMesh, color = "navy")
plotEdges(
  vertices, edges, color = "yellow", 
  tubesRadius = 0.01, spheresRadius = 0.02
)

vor <- VoronoiOnSphere(vertices)
plotVoronoiOnSphere(vor, colors = "random")

del <- DelaunayOnSphere(vertices)
plotDelaunayOnSphere(del)
