library(cgalMeshes)
library(rgl)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

mobiusStrip <- function(c, w, k){
  n <- 200
  t <- seq(-1, 1, len=n) 
  s <- -1
  r <- (1-w) + w*s*sin(k*pi*t/2)
  x1 <- r * sin(c*pi*t)/c
  y1 <- r*(1-(1-cos(c*pi*t))/c)
  z1 <- w*s*cos(k*pi*t/2)
  s <- 1
  r <- (1-w) + w*s*sin(k*pi*t/2)
  x2 <- r * sin(c*pi*t)/c
  y2 <- r*(1-(1-cos(c*pi*t))/c)
  z2 <- -z1
  a <- array(NA, dim=c(4,3,n-1))
  for(i in 1:(n-1)){
    a[,,i] <- rbind(c(x1[i],y1[i],z1[i]),
                    c(x1[i+1],y1[i+1],z1[i+1]),
                    c(x2[i+1],y2[i+1],z2[i+1]),
                    c(x2[i],y2[i],z2[i]))
  }
  c1 <- cbind(x1,y1,z1)
  c2 <- cbind(x2,y2,z2)
  return(list(strip=a, curve1=c1, curve2=c2))
}

strip <- mobiusStrip(c=1, w=1/4, k=2)
a <- aperm(strip$strip, c(2, 1, 3))
dim(a) <- c(3, 199*4)
rmesh <- Rvcg::vcgClean(as.tmesh3d(qmesh3d(
  vertices = a,
  indices = matrix(1:(199*4), nrow = 4)
)), sel = 0)


# rmesh <- Rvcg::vcgIsotropicRemeshing(rmesh, TargetLen = 0.01)
# rmesh <- Rvcg::vcgUpdateNormals(rmesh)
# rgl::shade3d(rmesh, color = "red")

mesh <- cgalMesh$new(rmesh)
mesh$isotropicRemeshing(0.01, iterations = 3, relaxSteps = 2)
mesh$writeMeshFile("torus.off")

M <- cgalMeshes:::testparam("torus.off", 3L)

clrs <- rep("yellow", nrow(M))
clrs[M[, 1] <= 0.25 & M[, 2] <= 0.25] <- "navy"
clrs[M[, 1] > 0.5 & M[, 1] <= 0.75 & M[, 2] <= 0.25] <- "navy"
clrs[M[, 1] > 0.25 & M[, 1] <= 0.5 & M[, 2] > 0.25 & M[, 2] <= 0.5] <- "navy"
clrs[M[, 1] > 0.75 & M[, 2] > 0.25 & M[, 2] <= 0.5] <- "navy"
clrs[M[, 1] <= 0.25 & M[, 2] > 0.5  & M[, 2] <= 0.75] <- "navy"
clrs[M[, 1] > 0.5 & M[, 1] <= 0.75 & M[, 2] > 0.5 & M[, 2] <= 0.75] <- "navy"
clrs[M[, 1] > 0.25 & M[, 1] <= 0.5 & M[, 2] > 0.75] <- "navy"
clrs[M[, 1] > 0.75 & M[, 2] > 0.75] <- "navy"
plot(M, type = "p", asp = 1, pch = "+", col=clrs, xlab = "u", ylab = "v")


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
  gif_file = "Enneper_colored_isotropicRemeshing.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

snapshot3d("Enneper_colored.png")
