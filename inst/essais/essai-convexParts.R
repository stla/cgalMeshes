library(cgalMeshes)
library(rgl)
mesh <- cgalMesh$new(pentagrammicPrism)$triangulate()
cxparts <- mesh$convexParts()
ncxparts <- length(cxparts)
colors <- hcl.colors(ncxparts, palette = "plasma")
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(20, -20, zoom = 0.8)
for(i in 1L:ncxparts) {
  cxmesh <- cxparts[[i]]$getMesh(normals = FALSE)
  shade3d(cxmesh, color = colors[i])
}

