library(cgalMeshes)
library(rgl)

refl <- function(A, B, C) {
  n <- rgl:::xprod(C-A, B-A)
  n <- n / sqrt(c(crossprod(n)))
  R <- diag(3) - 2 * tcrossprod(n)
  offset <- c(crossprod(n, B)) 
  M <- t(rbind(
    cbind(R, -2*offset*n), 
    c(0, 0, 0, 1)
  ))
  mesh1 <- sphericalTriangle(A, B, C)
  rmesh1 <- mesh1$getMesh()
  rmesh2 <- transform3d(rmesh1, M)
  rmesh2$normals <- -rmesh2$normals
  rmesh2
}

icosahedron <- icosahedron3d()
vertices    <- icosahedron[["vb"]][-4L, ]
faces       <- icosahedron[["it"]]
colors <- randomcoloR::distinctColorPalette(ncol(faces))
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.7)
for(i in 1L:ncol(faces)) {
  triangle <- faces[, i]
  A <- vertices[, triangle[1L]]
  B <- vertices[, triangle[2L]]
  C <- vertices[, triangle[3L]]
  rmesh <- refl(A, B, C)
  shade3d(rmesh, color = colors[i])
  wire3d(rmesh)
}


rmeshes <- vector("list", ncol(faces))
for(i in 1L:ncol(faces)) {
  triangle <- faces[, i]
  A <- vertices[, triangle[1L]]
  B <- vertices[, triangle[2L]]
  C <- vertices[, triangle[3L]]
  rmesh <- refl(A, B, C)
  rmeshes[[i]] <- Rvcg::vcgClean(rmesh, sel=0:7)
}

rmesh <- Rvcg::vcgClean(Morpho::mergeMeshes(rmeshes), sel=0:7)

mesh <- cgalMesh$new(rmesh)
#mesh$removeSelfIntersections()

ico <- cgalMesh$new(icosahedron3d())
mesh$clip(ico, FALSE)


i <- 1L
triangle <- faces[, i]
A <- vertices[, triangle[1L]]
B <- vertices[, triangle[2L]]
C <- vertices[, triangle[3L]]

mesh$clipToPlane(A, (B+C), FALSE)




movie3d(spin3d(axis = c(1, 1, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "gyrosphericalIcosahedron.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
