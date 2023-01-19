# the compound of five tetrahedra ####
library(rgl)

# all vertices
Vertices <- function(a, b, c, phi, O){
  t(cbind(
    c( a,  a,  a),
    c( a,  a, -a),
    c( a, -a,  a),
    c(-a, -a,  a),
    c(-a,  a, -a),
    c(-a,  a,  a),
    c( O,  b, -c),
    c( O, -b, -c),
    c( O, -b,  c),
    c( c,  O, -b),
    c(-c,  O, -b),
    c(-c,  O,  b),
    c( b,  c,  O),
    c( b, -c,  O),
    c(-b, -c,  O),
    c(-b,  c,  O),
    c( O,  b,  c),
    c( a, -a, -a),
    c( c,  O,  b),
    c(-a, -a, -a)
  ))
}

# 'double' vertex coordinates
phi <- (1 + sqrt(5)) / 2
a <- 1 / sqrt(3)
b <- a / phi
c <- a * phi
O <- 0
vertices <- Vertices(a, b, c, phi, O)

# the face indices are common to all tetrahedra
faces <- rbind(
  c(1L, 2L, 3L),
  c(3L, 2L, 4L),
  c(4L, 2L, 1L),
  c(1L, 3L, 4L)
)

# define the five tetrahedra meshes
mesh1 <- list(
  "vertices" = vertices[c(17, 14, 2, 11), ],
  "faces" = faces
)
mesh2 <- list(
  "vertices" = vertices[c(18, 1, 4, 5), ],
  "faces" = faces
)
mesh3 <- list(
  "vertices" = vertices[c(19, 6, 15, 7), ],
  "faces" = faces
)
mesh4 <- list(
  "vertices" = vertices[c(3, 13, 12, 8), ],
  "faces" = faces
)
mesh5 <- list(
  "vertices" = vertices[c(20, 16, 10, 9), ],
  "faces" = faces
)

meshes <- list(
  mesh1, mesh2, mesh3, mesh4, mesh5
)


tetrahedraCompound <- list(
  "Vertices" = list(
    mesh1[["vertices"]],
    mesh2[["vertices"]],
    mesh3[["vertices"]],
    mesh4[["vertices"]],
    mesh5[["vertices"]]
  ),
  
  "rglMeshes" = list(
    tmesh3d(
      "vertices"    = t(mesh1[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(mesh2[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(mesh3[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(mesh4[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    ),
    tmesh3d(
      "vertices"    = t(mesh5[["vertices"]]),
      "indices"     = t(faces),
      "homogeneous" = FALSE
    )
  )
)
