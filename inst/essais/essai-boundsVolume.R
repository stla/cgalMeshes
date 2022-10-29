library(cgalMeshes)


vertices1 <- rbind(
  c(-1, -1, -1),
  c(1, 1, -1),
  c(1, -1, 1),
  c(-1, 1, 1)
)
faces1 <- rbind(
  c(3, 2, 1),
  c(3, 4, 2),
  c(1, 2, 4),
  c(4, 3, 1)
)

vertices2 <- vertices1 + 3
faces2 <- faces1[, 3:1] + 4

vertices <- rbind(vertices1, vertices2)
faces <- rbind(faces1, faces2)

mesh <- cgalMesh$new(vertices = vertices, faces = faces)

vertices <- rbind(
  c(0, 0, 0),
  c(2, 2, 0),
  c(2, 0, 2),
  c(0, 2, 2),
  c(3, 3, 3),
  c(5, 5, 3),
  c(5, 3, 5),
  c(3, 5, 5)
)
faces <- rbind(
  c(3, 2, 1),
  c(3, 4, 2),
  c(1, 2, 4),
  c(4, 3, 1),
  c(5, 6, 7),
  c(6, 8, 7),
  c(8, 6, 5),
  c(5, 7, 8)
)
mesh <- cgalMesh$new(vertices = vertices, faces = faces)
mesh$boundsVolume() # FALSE
mesh$orientToBoundVolume()
mesh$boundsVolume() # TRUE