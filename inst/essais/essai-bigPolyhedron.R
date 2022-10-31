library(cgalMeshes)

meshFile <- system.file(
  "extdata", "bigPolyhedron.off", package = "cgalMeshes"
)
mesh <- cgalMesh$new(meshFile)
vertices <- mesh$vertices()

library(cxhull)
chull <- cxhull(vertices)
chullmesh <- hullMesh(chull)

mesh2 <- cgalMesh$new(chullmesh)

mesh2 <- cgalMesh$new("bigPolyhedron.off")
rglmesh2 <- mesh2$getMesh(normals = FALSE)

library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
shade3d(rglmesh2, color = "tomato")
cgalMeshes::plotEdges(
   mesh2$vertices(), mesh2$edges(), color = "darkred"
)