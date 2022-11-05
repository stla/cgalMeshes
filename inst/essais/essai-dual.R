library(cgalMeshes)

rgl_icosphere <- sphereMesh()
icosphere <- cgalMesh$new(rgl_icosphere)

dual <- icosphere$dual()
vertices <- dual$vertices()
edges <- dual$edges()
rgl_dual <- dual$triangulate()$getMesh()


library(rgl)
open3d(windowRect = 50 + c(0, 0, 900, 450))
mfrow3d(1L, 2L)
view3d(zoom = 0.6)
shade3d(rgl_icosphere, color = "navy")
wire3d(rgl_icosphere, color = "yellow", lwd = 2)
next3d()
view3d(zoom = 0.6)
shade3d(rgl_dual, color = "navy")
plotEdges(vertices, edges, edgesAsTubes = FALSE, color = "yellow")

snapshot3d("icosphere_dual.png", webshot = FALSE)

