library(cgalMeshes)
library(rgl)

sphericalCoordinates <- function(θ, ϕ){
  x <- cos(θ) * sin(ϕ)
  y <- sin(θ) * sin(ϕ)
  z <- cos(ϕ)
  rbind(x, y, z)
}

rmesh <- parametricMesh(
  sphericalCoordinates, urange = c(0, 2*pi), vrange = c(0, pi/2),
  periodic = c(TRUE, FALSE), nu = 512, nv = 512, clean = TRUE
)

u <- v <- seq(0, 1, length.out = 513L)[-513L]
Grid <- expand.grid(U = u, V = v)
checkerboard <- ifelse(
  (floor(10*Grid$U) %% 2) == (floor(10*Grid$V) %% 2), 
  "yellow", "navy"
)
rmesh$material <- list(color = checkerboard)

shade3d(rmesh, meshColor = "vertices")
