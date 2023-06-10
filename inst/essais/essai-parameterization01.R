library(cgalMeshes)
library(rgl)

n <- 3
Enneper <- function(phi, r) {
  rbind(
    r*cos(phi) - r^(2*n-1)*cos((2*n-1)*phi)/(2*n-1),
    r*sin(phi) + r^(2*n-1)*sin((2*n-1)*phi)/(2*n-1),
    2*r^n*cos(n*phi)/n
  )
}

rmesh <- parametricMesh(
  Enneper, urange = c(0, 2*pi), vrange = c(0, 1.3),
  periodic = c(TRUE, FALSE), nu = 400, nv = 200, clean = TRUE
)

mesh <- cgalMesh$new(rmesh)

UV <- mesh$parameterization("IAP", iterations = 2)

plot(UV, type = "p", asp = 1, pch = ".")