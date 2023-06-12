library(cgalMeshes)
library(rgl)

Enneper <- function(phi, r, n = 4) {
  rbind(
    r*cos(phi) - r^(2*n-1) * cos((2*n-1)*phi) / (2*n-1),
    r*sin(phi) + r^(2*n-1) * sin((2*n-1)*phi) / (2*n-1),
    2*r^n * cos(n*phi)/n
  )
}

rmesh <- parametricMesh(
  Enneper, urange = c(0, 2*pi), vrange = c(0, 1.1),
  periodic = c(TRUE, FALSE), nu = 1024, nv = 1024, clean = FALSE
)

mesh <- cgalMesh$new(rmesh)
# compute the discrete conformal parameterization with square border
UV <- mesh$parameterization("DCP", UVborder = "square")

# now we extract the colors from the image
library(imager)
# load the image
img <- load.image("schmidt2_square.png")
# take the r, g, b channels
r <- squeeze(R(img))
g <- squeeze(G(img))
b <- squeeze(B(img))

# make interpolation functions to get the colors of the UV points
library(cooltools) # to get the `approxfun2` function
x_ <- y_ <- seq(0, 1, length.out = 1024L)
f_r <- approxfun2(x_, y_, r)
f_g <- approxfun2(x_, y_, g)
f_b <- approxfun2(x_, y_, b)

# now, interpolate the r, g, b values
UV_r <- f_r(UV[, 1L], UV[, 2L])
UV_g <- f_g(UV[, 1L], UV[, 2L])
UV_b <- f_b(UV[, 1L], UV[, 2L])

# convert rgb to hex codes
clrs <- rgb(UV_r, UV_g, UV_b)

# we're done; now compute mesh normals, convert to rgl mesh, and assign colors
mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

# plot
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-10, -35, zoom = 0.7)
shade3d(rmesh, meshColor = "vertices")

# animation ####
movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "Enneper-Schmidt.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))
