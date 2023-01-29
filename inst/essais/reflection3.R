library(cgalMeshes)
library(rgl)

icosahedron <- icosahedron3d()
vertices    <- icosahedron[["vb"]][-4L, ]
faces       <- icosahedron[["it"]]

mesh <- cgalMesh$new(icosahedron)

i <- 1L
triangle <- faces[, i]
A <- vertices[, triangle[1L]]
B <- vertices[, triangle[2L]]
C <- vertices[, triangle[3L]]

n <- rgl:::xprod(C-A, B-A)
n <- n / sqrt(c(crossprod(n)))
offset <- c(crossprod(n, B)) 
v <- -2*offset*n

phi <- (1+sqrt(5))/2
r <- sqrt(1 + (phi - 1)^2)

clipper <- cgalMesh$new(sphereMesh(v[1], v[2], v[3], r))

mesh$clip(clipper, FALSE)
rmesh <- mesh$getMesh()

shade3d(rmesh, col = "red")



