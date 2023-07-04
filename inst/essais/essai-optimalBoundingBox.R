library(cgalMeshes)
library(rgl)
cyclide <- cyclideMesh(a = 97, c = 32, mu = 57)
mesh <- cgalMesh$new(cyclide)
x <- mesh$optimalBoundingBox()
obb <- x[["rmesh"]]

open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
shade3d(cyclide, color = "green")
shade3d(obb, color = "yellow", alpha = 0.3)

# tetrahedra ####
faces <- cbind(
	c(1L, 2L, 3L),
	c(3L, 2L, 4L),
	c(4L, 2L, 1L),
	c(1L, 3L, 4L)
)

pts <- x[["hxVertices"]]
vs1 <- pts[, c(1L, 5L, 3L, 7L)]
vs2 <- pts[, c(7L, 1L, 2L, 3L)]
vs3 <- pts[, c(6L, 1L, 7L, 5L)]
vs4 <- pts[, c(8L, 7L, 3L, 5L)]
vs5 <- pts[, c(4L, 1L, 3L, 5L)]

tth1 <- tmesh3d(
	vertices = vs1,
	indices  = faces
)
tth2 <- tmesh3d(
	vertices = vs2,
	indices  = faces
)
tth3 <- tmesh3d(
	vertices = vs3,
	indices  = faces
)
tth4 <- tmesh3d(
	vertices = vs4,
	indices  = faces
)
tth5 <- tmesh3d(
	vertices = vs5,
	indices  = faces
)

open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
shade3d(tth1, color = "red")
shade3d(tth2, color = "green")
shade3d(tth3, color = "blue")
shade3d(tth4, color = "yellow")
shade3d(tth5, color = "gray")
wire3d(obb, color = "black")

