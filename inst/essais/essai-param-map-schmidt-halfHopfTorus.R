library(rgl)
library(cgalMeshes)

HT <- function(t, phi, h = 0.4, nlobes = 4, alpha = NULL){
  p2 <- cos(t) * cos(h * cos(nlobes*t))
  p3 <- sin(t) * cos(h * cos(nlobes*t))
  p1 <- sin(h * cos(nlobes*t))
  ## alternatively, 
  # den = sqrt(1+h^2*cos(n*t)^2);
  # p2 = cos(t)/den;
  # p3 = sin(t)/den;
  # p1 = h*cos(n*t)/den;
  ##
  yden <- sqrt(2 * (1 + p1))
  y1 <- (1 + p1) / yden
  y2 <- p2 / yden
  y3 <- p3 / yden
  cosphi <- cos(phi)
  sinphi <- sin(phi)
  x1 <- cosphi*y3 + sinphi*y2
  x2 <- cosphi*y2 - sinphi*y3
  x3 <- sinphi*y1
  x4 <- cosphi*y1
  if(is.null(alpha)) {
    t(cbind(x1, x2, x3)/(1 - x4))
  }else{
    t(acos(x4)/(1 - abs(x4)^alpha)^(1/alpha) * cbind(x1, x2, x3))
  }
}

rmesh <- parametricMesh(HT, c(-pi, pi), c(0, pi), c(TRUE, FALSE), 1024, 1024)

shade3d(rmesh, col = "red")

mesh <- cgalMesh$new(rmesh)
UV <- mesh$parameterization("DCP", "square")


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
view3d(0, 0, zoom = 0.7)
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
  gif_file = "halfHopfTorus-Schmidt.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))



