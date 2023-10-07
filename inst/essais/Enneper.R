library(cgalMeshes)
library(rgl)

n <- 2
Enneper <- function(phi, r) {
  rbind(
    r*cos(phi) - r^(2*n-1)*cos((2*n-1)*phi)/(2*n-1),
    r*sin(phi) + r^(2*n-1)*sin((2*n-1)*phi)/(2*n-1),
    2*r^n*cos(n*phi)/n
  )
}

rmesh <- parametricMesh(
  Enneper, urange = c(0, 2*pi), vrange = c(0, 1.2),
  periodic = c(TRUE, FALSE), nu = 200, nv = 100, clean = TRUE
)
#rmesh <- Rvcg::vcgClean(rmesh, sel = 0)

shade3d(rmesh)