library(cgalMeshes)
library(rgl)
library(trekcolors)
fpalette <- colorRamp(
  viridisLite::inferno(255), bias = 1, interpolate = "spline"
)

t <- seq(0, 2, length.out = 1000)[-1L]
crcl1 <- cbind(cospi(t)+0.5, sinpi(t), 0)
crcl2 <- cbind(cospi(t)-0.5, 0, sinpi(t))
points <- rbind(crcl1, crcl2)

mesh <- convexHull(points)
#mesh$computeNormals()

rmesh <- mesh$getMesh()
d2 <- apply(rmesh$vb[-4L, ], 2L, crossprod)
d2 <- (d2 - min(d2)) / diff(range(d2))
RGB <- fpalette(d2)
rmesh[["material"]] <- 
  list(color = rgb(RGB[, 1L], RGB[, 2L], RGB[, 3L], maxColorValue = 255))

open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.85)
shade3d(rmesh)

# # -- if you want an animation 
M <- par3d("userMatrix")
movie3d(
  par3dinterp(
    time = seq(0, 1, len = 9),
    userMatrix = list(
      M,
      rotate3d(M, pi, 1, 0, 0),
      rotate3d(M, pi, 1, 1, 0),
      rotate3d(M, pi, 1, 1, 1),
      rotate3d(M, pi, 0, 1, 1),
      rotate3d(M, pi, 0, 1, 0),
      rotate3d(M, pi, 1, 0, 1),
      rotate3d(M, pi, 0, 0, 1),
      M
    )
  ),
  fps = 100,
  duration = 1,
  dir = ".",
  movie = "cc",
  convert = FALSE,
  webshot = FALSE
)

library(gifski)
pngs <- list.files(pattern = "^cc.*png$")
gifski(
  pngs, "oloid.gif", width = 512, height = 512, delay = 1/10
)

file.remove(pngs)
