library(cgalMeshes)

# a tetrahedron with ill-oriented faces ####
vertices <- rbind(
  c(-1, -1, -1),
  c(1, 1, -1),
  c(1, -1, 1),
  c(-1, 1, 1)
)
faces <- rbind(
  c(1, 2, 3),
  c(3, 4, 2),
  c(4, 2, 1),
  c(4, 3, 1)
)

# mesh <- cgalMesh$new(vertices = vertices, faces = faces)
# mesh$writeMeshFile("illOrientedTetrahedron.off")

mesh <- cgalMesh$new("illOrientedTetrahedron.off")
