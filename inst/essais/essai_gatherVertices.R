
M <- cbind(
  c(1, 2, 3),
  c(4, 5, 6),
  c(7, 8, 9),
  c(4, 5, 6),
  c(6, 5, 4),
  c(7, 8, 9)
)

cgalMeshes:::duplicated_vertices(M)