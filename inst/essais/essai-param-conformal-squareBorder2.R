library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

alpha <- 0.35
beta <- 0.35
gamma <- 0.1
n <- 1
f <- function(u, v) {
  rbind(
    alpha * (1 - v/(2*pi)) * cos(n*v) * (1 + cos(u)) + gamma * cos(n * v),
    alpha * (1 - v/(2*pi)) * sin(n*v) * (1 + cos(u)) + gamma * sin(n * v),
    alpha * (1 - v/(2*pi)) * sin(u) + beta * v/(2*pi)
  )
}

rmesh <- parametricMesh(
  f, c(0, 2*pi), c(0, 2*pi), periodic = c(TRUE, FALSE), nu = 400, nv = 300, clean = FALSE
)
rmesh <- Rvcg::vcgClean(rmesh, 0)

mesh <- cgalMesh$new(rmesh)
summary(mesh$getEdges())
mesh$isotropicRemeshing(0.0025, iterations = 3, relaxSteps = 2)
summary(mesh$getEdges())
mesh
mesh$writeMeshFile("horn2.off")

M <- cgalMeshes:::testparam(normalizePath("horn2.off"), 5L)

plot(M, type = "p", pch = ".", asp = 1)


# now we extract the colors
library(imager)
# load the image
img <- load.image("missPsychoCat512x512.png")
# if the image is gray, add colors
#img <- add.color(img) 
# take the r, g, b channels
r <- squeeze(R(img))
g <- squeeze(G(img))
b <- squeeze(B(img))

# interpolation
library(cooltools) # to get the `approxfun2` function
x_ <- y_ <- seq(0, 1, length.out = 512L)
f_r <- approxfun2(x_, y_, r)
f_g <- approxfun2(x_, y_, g)
f_b <- approxfun2(x_, y_, b)

M_r <- f_r(M[, 1L], M[, 2L])
M_g <- f_g(M[, 1L], M[, 2L])
M_b <- f_b(M[, 1L], M[, 2L])
# convert to hex codes
clrs <- rgb(M_r, M_g, M_b)


plot(M, type = "p", asp = 1, pch = ".", col=clrs, xlab = "u", ylab = "v")

mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(90, 90, zoom = 0.9)
shade3d(rmesh, meshColor = "vertices")

movie3d(spin3d(axis = c(1, 0, 0), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "horn_missPsychoCat.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

snapshot3d("Enneper_colored.png")
