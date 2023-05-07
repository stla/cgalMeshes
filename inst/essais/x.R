library(cgalMeshes)
library(rgl)
rmesh <- sphereMesh()
mesh <- cgalMesh$new(rmesh)
nfaces <- nrow(mesh$getFaces())
if(require("randomcoloR")) {
  colors <- 
    randomColor(nfaces, hue = "random", luminosity = "dark")
} else {
  colors <- rainbow(nfaces)
}
mesh$assignFaceColors(colors)
meshes <- mesh$clipToPlane(
  planePoint  = c(0, 0, 0), 
  planeNormal = c(0, 0, 1), 
  clipVolume = TRUE
)
mesh1 <- meshes[[1]]
mesh2 <- meshes[[2]]
mesh1$computeNormals()
rClippedMesh1 <- mesh1$getMesh()
rClippedMesh2 <- mesh2$getMesh()
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(70, 0)
shade3d(rClippedMesh1, meshColor = "faces")
shade3d(rClippedMesh2, color = "orange")