library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

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
summary(mesh$getEdges())
mesh$isotropicRemeshing(0.01, iterations = 3, relaxSteps = 2)
mesh$writeMeshFile("enneper.off")

M <- cgalMeshes:::testparam(normalizePath("enneper.off"), 5L)

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


rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
shade3d(rmesh, meshColor = "vertices")

movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "Enneper_missPsychoCat.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

snapshot3d("Enneper_colored.png")
