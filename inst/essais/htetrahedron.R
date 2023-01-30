library(cgalMeshes)
library(rgl)

r <- 1/1.25^2
r <- 1/sqrt(3)

tetrahedron <- tetrahedraCompound$rglMeshes[[1]]
vertices    <- tetrahedron[["vb"]][-4L, ]
faces       <- tetrahedron[["it"]]
mesh <- cgalMesh$new(tetrahedron)

for(i in 1:ncol(faces)) {
  triangle <- faces[, i]
  A <- vertices[, triangle[1L]]
  B <- vertices[, triangle[2L]]
  C <- vertices[, triangle[3L]]
  n <- rgl:::xprod(C-A, B-A)
  n <- n / sqrt(c(crossprod(n)))
  offset <- c(crossprod(n, B)) 
  v <- -2*offset*n
  clipper <- cgalMesh$new(sphereMesh(v[1], v[2], v[3], r, iterations = 4))
  mesh <- mesh$subtract(clipper)
}


rmesh <- mesh$getMesh()
shade3d(rmesh, col = "red")

sharpEdges <- mesh$sharpEdges(120)

mesh$computeNormals()
rmesh <- mesh$getMesh()

open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.7)
shade3d(rmesh, col = "springgreen")
plotEdges(mesh$getVertices(), sharpEdges[, c("v1", "v2")], 
          edgesAsTubes = TRUE, tubesRadius = 0.011,
          verticesAsSpheres = FALSE)


