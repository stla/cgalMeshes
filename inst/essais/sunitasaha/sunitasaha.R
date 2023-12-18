library(cgalMeshes)
library(rgl)

pts <- as.matrix(read.table("surface_points.txt"))

mesh <- PoissonReconstruction(pts)
mesh$computeNormals()
rmesh <- mesh$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(20, -40, zoom = 1)
shade3d(rmesh, color = "darkorange")
wire3d(rmesh)

mesh <- AFSreconstruction(pts)
mesh$computeNormals()
rglMesh <- mesh$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512))
shade3d(rglMesh, color = "firebrick")


mesh <- alphaWrap(pts, 20, 3000)
mesh$computeNormals()
rglMesh <- mesh$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512))
shade3d(rglMesh, color = "forestgreen", polygon_offset = 1)
wire3d(rglMesh)


mesh <- SSSreconstruction(
  pts, scaleIterations = 10, 
  forceManifold = FALSE, neighbors = 40, samples = 300,
  separateShells = FALSE, borderAngle = 20
)
mesh$computeNormals()
rglMesh <- mesh$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(20, -40, zoom = 0.85)
shade3d(rglMesh, color = "tomato", polygon_offset = 1)
wire3d(rglMesh)