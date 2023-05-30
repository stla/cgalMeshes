library(cgalMeshes)
library(rgl)

icosahedron <- icosahedron3d()
vertices    <- icosahedron[["vb"]][-4L, ]
faces       <- icosahedron[["it"]]
colors <- rainbow(ncol(faces))

open3d(windowRect = 50 + c(0, 0, 512, 512))
for(i in 1L:ncol(faces)) {
  triangle <- faces[, i]
  A <- vertices[, triangle[1L]]
  B <- vertices[, triangle[2L]]
  C <- vertices[, triangle[3L]]
  mesh <- sphericalTriangle(A, B, C)
  rmesh <- mesh$getMesh()
  shade3d(rmesh, color = colors[i])
  wire3d(rmesh)
}

phi <- (1+sqrt(5))/2
r <- sqrt((phi-1)^2 + 1)
open3d(windowRect = 50 + c(0, 0, 512, 512))
for(i in 1L:ncol(faces)) {
  triangle <- faces[, i]
  A <- vertices[, triangle[1L]]
  B <- vertices[, triangle[2L]]
  C <- vertices[, triangle[3L]]
  mesh <- cgalMesh$new(vertices = rbind(A, B, C)/r, faces = rbind(c(1L, 2L, 3L)))
  mesh$LoopSubdivision(iterations = 3)
  rmesh <- mesh$getMesh()
  shade3d(rmesh, color = colors[i])
  wire3d(rmesh)
}