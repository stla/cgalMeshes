library(cgalMeshes)

n <- 3
Enneper <- function(phi, r) {
  cbind(
    r*cos(phi) - r^(2*n-1)*cos((2*n-1)*phi)/(2*n-1),
    r*sin(phi) + r^(2*n-1)*sin((2*n-1)*phi)/(2*n-1),
    2*r^n*cos(n*phi)/n
  )
}

library(rgl)
parametricMesh3d <- function(
    Fxyz, umin, umax, vmin, vmax, nu, nv 
){
  u <- seq(umin, umax, length.out = nu)
  v <- seq(vmin, vmax, length.out = nv)
  tg <- misc3d:::expandTriangleGrid(u, v)
  f <- function(uv) Fxyz(uv[, 1L], uv[, 2L])
  v1 <- f(tg$v1)
  v2 <- f(tg$v2)
  v3 <- f(tg$v3)
  tris <- misc3d::makeTriangles(v1, v2, v3)
  mesh0 <- misc3d:::t2ve(tris)
  Rvcg::vcgUpdateNormals(
    tmesh3d(
      vertices = mesh0$vb,
      indices = mesh0$ib
    )
  )
}

rmesh <- parametricMesh3d(Enneper, 0, 2*pi, 0, 1.3, 300, 200)

rmesh <- Rvcg::vcgClean(rmesh, sel = 0)

mesh <- cgalMesh$new(rmesh)
#mesh$Sqrt3Subdivision()
mesh$writeMeshFile("torus.off")

M <- cgalMeshes:::testparam()

clrs <- rep("yellow", nrow(M))
clrs[M[, 1] <= 0.25 & M[, 2] <= 0.25] <- "navy"
clrs[M[, 1] > 0.5 & M[, 1] <= 0.75 & M[, 2] <= 0.25] <- "navy"
clrs[M[, 1] > 0.25 & M[, 1] <= 0.5 & M[, 2] > 0.25 & M[, 2] <= 0.5] <- "navy"
clrs[M[, 1] > 0.75 & M[, 2] > 0.25 & M[, 2] <= 0.5] <- "navy"
clrs[M[, 1] <= 0.25 & M[, 2] > 0.5  & M[, 2] <= 0.75] <- "navy"
clrs[M[, 1] > 0.5 & M[, 1] <= 0.75 & M[, 2] > 0.5 & M[, 2] <= 0.75] <- "navy"
clrs[M[, 1] > 0.25 & M[, 1] <= 0.5 & M[, 2] > 0.75] <- "navy"
clrs[M[, 1] > 0.75 & M[, 2] > 0.75] <- "navy"
plot(M, type = "p", asp = 1, pch = 19, col=clrs, xlab = "u", ylab = "v")


rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
shade3d(rmesh, meshColor = "vertices")

movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "Enneper_colored.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

snapshot3d("Enneper_colored.png")
