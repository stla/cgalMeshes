library(rgl)
library(cgalMeshes)

N <- 1024

g <- function(phi) {
  r <- sin(phi) / sqrt(1 - sin(phi)^2)
  R <- cos(phi) / (1 - sin(phi)) - r
  h <- R/r
  sqrt(h*h - 1)
}

phi <- pi/4 # c'est la root de g=1
r <- sin(phi) / sqrt(1 - sin(phi)^2) # 1
R <- cos(phi) / (1 - sin(phi)) - r   # sqrt(2)

torus <- function(R, r) {
  h <- R/r
  s <- sqrt(h*h - 1)
  r <- s*r
  f <- function(u, v) {
    w <- h - cospi(2*v)
    rbind(
      s * cospi(2*u/s) / w,
      s * sinpi(2*u/s) / w,
      sinpi(2*v) / w
    ) * r
  }
  parametricMesh(
    f, c(0, s), c(0, 1), c(TRUE, TRUE),
    nu = N, nv = N
  )
}

mesh <- addNormals(torus(sqrt(2), 1))
summary(t(mesh$vb))
# r = 1/sqrt((h-1)*(h+1))


x_ <- seq(0, 1, length.out = N)
UV <- as.matrix(
  expand.grid(U = x_, V = x_)
)

rotation <- function(alpha, uv) {
  t(rbind(
    c(cos(alpha), -sin(alpha)),
    c(sin(alpha),  cos(alpha))
  ) %*% t(uv))
}

UVrot <- rotation(pi/4, UV)

rotatedCheckerboardColors <- ifelse(
  (floor(4*sqrt(2)*UVrot[, 1L]) %% 2) == (floor(4*sqrt(2)*UVrot[, 2L]) %% 2),
  "yellow", "navy"
)

png("rotatedCheckerboard.png", width = 384, height = 384)
opar <- par(mar = c(0, 0, 0, 0))
plot(
  UV, col = rotatedCheckerboardColors, asp = 1, pch = ".", 
  xaxs = "i", yaxs = "i", axes = FALSE
)
dev.off()

mesh$material <- list(color = rotatedCheckerboardColors)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-15, -25, zoom = 0.7)
shade3d(mesh, polygon_offset = 1) 
#bbox3d()
snapshot3d("conformalTorus.png", webshot = FALSE)

movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 1.5, fps = 20,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/20)

pngs <- Sys.glob("zzpic*.png")
gifski::gifski(
  pngs,
  "conformalTorus.gif",
  width = 512, height = 512,
  delay = 1/10
)

Villarceau <- function(beta, theta0) {
  d <- (1 - sin(beta) * sin(phi))
  cbind(
    cos(theta0 - beta) * cos(phi) / d,
    sin(theta0 - beta) * cos(phi) / d,
    cos(beta) * sin(phi) / d
  ) 
}

beta_ <- seq(0, 2*pi, len = 400)
pts <- Villarceau(beta_, theta0 = 0)
tube <- addNormals(
  cylinder3d(pts, radius = 0.04, sides = 30)
)
shade3d(tube, color = "black")
pts <- Villarceau(beta_, theta0 = pi/2)
tube <- addNormals(
  cylinder3d(pts, radius = 0.04, sides = 30)
)
shade3d(tube, color = "black")
pts <- Villarceau(beta_, theta0 = pi)
tube <- addNormals(
  cylinder3d(pts, radius = 0.04, sides = 30)
)
shade3d(tube, color = "black")
pts <- Villarceau(beta_, theta0 = 3*pi/2)
tube <- addNormals(
  cylinder3d(pts, radius = 0.04, sides = 30)
)
shade3d(tube, color = "black")


Villarceau <- function(beta, theta0) {
  d <- (1 - sin(beta) * sin(phi))
  cbind(
    cos(theta0 + beta) * cos(phi) / d,
    sin(theta0 + beta) * cos(phi) / d,
    cos(beta) * sin(phi) / d
  ) 
}

beta_ <- seq(0, 2*pi, len = 400)
pts <- Villarceau(beta_, theta0 = 0)
tube <- addNormals(
  cylinder3d(pts, radius = 0.04, sides = 30)
)
shade3d(tube, color = "black")
pts <- Villarceau(beta_, theta0 = pi/2)
tube <- addNormals(
  cylinder3d(pts, radius = 0.04, sides = 30)
)
shade3d(tube, color = "black")
pts <- Villarceau(beta_, theta0 = pi)
tube <- addNormals(
  cylinder3d(pts, radius = 0.04, sides = 30)
)
shade3d(tube, color = "black")
pts <- Villarceau(beta_, theta0 = 3*pi/2)
tube <- addNormals(
  cylinder3d(pts, radius = 0.04, sides = 30)
)
shade3d(tube, color = "black")

snapshot3d("Villarceau.png", webshot = FALSE)



################################################################################

f <- function(beta, theta0, phi = pi/4) {
  d <- (1 - sin(beta) * sin(phi))
  rbind(
    cos(theta0 + beta) * cos(phi) / d,
    sin(theta0 + beta) * cos(phi) / d,
    cos(beta) * sin(phi) / d
  ) 
}

mesh <- parametricMesh(
  f, c(0, 2*pi), c(0, 2*pi), c(TRUE, TRUE),
  nu = 256, nv = 256
)

f <- function(u, v) { # clifford
  rbind(
    cospi(u) / (sqrt(2) - sinpi(v)),
    sinpi(u) / (sqrt(2) - sinpi(v)),
    cospi(v) / (sqrt(2) - sinpi(v))
  ) 
}

N <- 256

mesh <- parametricMesh(
  f, c(0, 2), c(0, 2), c(TRUE, TRUE),
  nu = 256, nv = 256
)

x_ <- seq(0, 1, length.out = N)
UV <- as.matrix(
  expand.grid(U = x_, V = x_)
)

rotation <- function(alpha, uv) {
  t(rbind(
    c(cos(alpha), -sin(alpha)),
    c(sin(alpha),  cos(alpha))
  ) %*% t(uv))
}

UVrot <- rotation(pi/4, UV)

K <- 4*sqrt(2)
rotatedCheckerboardColors <- ifelse(
  (floor(K*UVrot[, 1L]) %% 2) == (floor(K*UVrot[, 2L]) %% 2),
  "yellow", "navy"
)

mesh$material <- list(color = rotatedCheckerboardColors)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-15, -25, zoom = 0.7)
shade3d(mesh, polygon_offset = 1) 
