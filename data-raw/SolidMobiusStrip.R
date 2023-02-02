a = 0.4
b = 0.1
f <- function(x, y, z, a, b){
  ((x*x+y*y+1)*(a*x*x+b*y*y)+z*z*(b*x*x+a*y*y)-2*(a-b)*x*y*z-a*b*(x*x+y*y))^2 -
    4*(x*x+y*y)*(a*x*x+b*y*y-x*y*z*(a-b))^2
}

# run the marching cubes algorithm ####
nx <- 120; ny <- 120; nz <- 120
x <- seq(-1.4, 1.4, length.out = nx)
y <- seq(-1.7, 1.7, length.out = ny)
z <- seq(-0.7, 0.7, length.out = nz)
G <- expand.grid(x = x, y = y, z = z)
voxel <- array(with(G, f(x, y, z, a, b)), c(nx, ny, nz))

library(misc3d)
surf <- computeContour3d(
  voxel, maxvol = max(voxel), level = 0, x = x, y = y, z = z
)

set.seed(666L)
SolidMobiusStrip <- surf[sample.int(nrow(surf), 10000L), ]
Noise <- matrix(rnorm(30000L, sd = 0.005), nrow = 10000L, ncol = 3L)
SolidMobiusStrip <- SolidMobiusStrip + Noise
