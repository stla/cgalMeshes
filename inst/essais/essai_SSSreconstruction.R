library(cgalMeshes)
mesh <- SSSreconstruction(
  SolidMobiusStrip, scaleIterations = 4, 
  forceManifold = TRUE, neighbors = 30
)
mesh$computeNormals()
rglMesh <- mesh$getMesh()
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(20, -40, zoom = 0.85)
shade3d(rglMesh, color = "tomato")