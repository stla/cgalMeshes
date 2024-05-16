R <- 1
r <- 1/sqrt(2)
R^2/r^2-1
x1 <- R+r + .8
x2 <- R-r + .8
x3 <- -R+r + .8
x4 <- -R-r + .8

ar <- 1
R <- 2
r <- sqrt(R^2/(ar^2 + 1))
2*R-2*r
( x2 <- 2*R-2*r - 0.01 )
# x2 - 2*R < - 1
# x2-2*R < -2*r
#R <- (x2+1)/2 + 1
x1 <- x2+2*r
(x3 <- x1-2*R)
(x4 <- x2-2*R)
R^2/r^2-1

R <- (x1-x3)/2
r <- (x1-x2)/2
ratio <- (1-x3/x1)/(1-x2/x1)


invx1 <- 0.3 # = radius^2/x1
invx2 <- 0.8
radius <- 0.81
x1 <- radius^2/invx1
x2 <- radius^2/invx2
ar <- 1
# => ratio = sqrt(ar^2+1)
ratio <- sqrt(ar^2+1)
(x3 <- x1 -  ratio * (x1-x2))
( x3 <- x3_over_x1 / invx1 )
# shift <- 0.5
# x3 <- x3 - shift
r <- (1/invx1 - 1/invx2) / 2
R <- r * ratio
(x4 <- x3 - 2*r)
x1 <- 1/invx1
x2 <- 1/invx2




ar <- 1
d <- 0.5 # petit diamètre
r <- 2/(ar^2*d) 
R <- sqrt(r^2*(1+ar^2))
# R-r = r*(sqrt(1+ar^2)-1) = 2*(sqrt(1+ar^2)-1)/(ar^2*d)
# R-r < 1 pour que x3 > -1, mais x2=R-r >1  
R_minus_r <- R - r
# r <- r / (0.5*R_minus_r)
# R <- R / (0.5*R_minus_r)
R - r
( x2 <- 1 ) #2*R-2*r - (R-r) )
# x2 - 2*R < - 1
# x2-2*R < -2*r
#R <- (x2+1)/2 + 1
x1 <- x2+2*r
(x3 <- x1-2*R)
(x4 <- x2-2*R)
1/x2 - 1/x1 # = 2*r / (x2*(x2+2*r)) = d => d*x2*(x2+2*r) = 2*r

a <- (1/x1+1/x2-1/x3-1/x4)/4
c <- (1/x1-1/x2-1/x3+1/x4)/4
mu <- (-1/x1+1/x2-1/x3+1/x4)/4
a;c;mu

a = 97/100; c = 32/100; mu = 57/100
b2 <- a * a - c * c
bb <- sqrt(b2 * (mu * mu - c * c))
h <- (c * c) / ((a - c) * (mu - c) + bb)
r <- (h * (mu - c)) / ((a + c) * (mu - c) + bb)
R <- (h * (a - c)) / ((a - c) * (mu + c) + bb)
R^2/r^2-1

library(rgl)
library(cgalMeshes)
cyclide <- function(a, c, mu, nu = 90L, nv = 40L){
  stopifnot(c > 0, a > mu, mu > c)
  stopifnot(nu >= 3, nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  vertices <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
  normals  <- matrix(NA_real_, nrow = nu*nv, ncol = 3L)
  b2 <- a * a - c * c
  bb <- sqrt(b2 * (mu * mu - c * c))
  omega <- (a * mu + bb) / c
  Omega0 <- c(omega, 0, 0)
  inversion <- function(M) {
    OmegaM <- M - Omega0
    k <- c(crossprod(OmegaM))
    OmegaM / k + Omega0
  }
  h <- (c * c) / ((a - c) * (mu - c) + bb)
  r <- (h * (mu - c)) / ((a + c) * (mu - c) + bb)
  R <- (h * (a - c)) / ((a - c) * (mu + c) + bb)
  bb2 <- b2 * (mu * mu - c * c)
  denb1 <- c * (a*c - mu*c + c*c - a*mu - bb)
  b1 <- (a*mu*(c-mu)*(a+c) - bb2 + c*c + bb*(c*(a-mu+c) - 2*a*mu))/denb1
  denb2 <- c * (a*c - mu*c - c*c + a*mu + bb)
  b2 <- (a*mu*(c+mu)*(a-c) + bb2 - c*c + bb*(c*(a-mu-c) + 2*a*mu))/denb2
  omegaT <- (b1 + b2)/2
  OmegaT <- c(omegaT, 0, 0)
  tormesh <- torusMesh(R, r, nu = nu, nv = nv, conformal = TRUE)
  rtnormals <- r * tormesh[["normals"]][1L:3L, ]
  xvertices <- tormesh[["vb"]][1L:3L, ] + OmegaT
  for(i in 1L:nu){
    k0 <- i * nv - nv
    for(j in 1L:nv){
      k <- k0 + j
      rtnormal <- rtnormals[, k]
      xvertex <- xvertices[, k]
      vertex <- inversion(xvertex)
      vertices[, k] <- vertex
      foo <- vertex - inversion(rtnormal + xvertex)
      normals[k, ] <- foo / sqrt(c(crossprod(foo)))
    }
  }
  tmesh3d(
    vertices    = vertices,
    indices     = tormesh[["it"]],
    normals     = normals,
    homogeneous = FALSE
  )
}

mesh <- cyclide(a, c, mu, nu=256, nv=256)

library(jacobi)
library(RcppColors)
# generate the values of wp on a grid
x <- y <- seq(-4, 0, length.out = 256L)
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, omega = c(0.5, 0.5 + 0.5i))
}
Z <- outer(x, y, f)
# map them to colors
img <- colorMap5(Z)
clrs <- c(img)

mesh[["material"]] <- list("color" = clrs)
# plot
open3d(windowRect = 50 + c(0, 0, 512, 512))
clear3d(type = "lights") # destroy current lights
light3d(x = -50, y = 100, z = 100)
bg3d("#363940")
view3d(-25, -25, zoom = 0.75)
shade3d(mesh, specular = "gold")


##########"

conformalCyclideMesh <- function(
  a, c, aspectRatio, normals = TRUE, nu = 90L, nv = 40L
){
  stopifnot(c > 0, a > c)
  stopifnot(nu >= 3, nv >= 3)
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  vertices <- matrix(NA_real_, nrow = 3L, ncol = nu*nv)
  b2 <- a*a - c*c
  ratio <- sqrt(aspectRatio^2 + 1)
  mu <- sqrt(c^2 + b2/(aspectRatio^2 + 1))
  bb <- b2/sqrt(aspectRatio^2 + 1)
  bb2 <- b2 * b2 / (aspectRatio^2 + 1)
  omega <- (a * mu + bb) / c
  Omega0 <- c(omega, 0, 0)
  inversion <- function(M) {
    OmegaM <- M - Omega0
    k <- c(crossprod(OmegaM))
    OmegaM / k + Omega0
  }
  h <- (c * c) / ((a - c) * (mu - c) + bb)
  r <- (h * (mu - c)) / ((a + c) * (mu - c) + bb)
  R <- ratio * r
  # a' = a/c; mu' = mu/c
  # R/r = (a'-1)/(mu'-1) * ((a'+1)*(mu'-1) + sqrt((a'^2-1)*(mu'^2-1))) / ((a'-1)*(mu'+1)+sqrt((a'^2-1)*(mu'^2-1)))
  # = (a'-1)/(mu'-1) * ((a'+1)*(mu'-1)/(a'-1) + sqrt((mu'^2-1)*(a'+1)/(a'-1))) / ((mu'+1) + sqrt((mu'^2-1)*(a'+1)/(a'-1)))
  # = ((aa+1) * (1 + sqrt((aa-1)/(aa+1))*sqrt((muu+1)/(muu-1))) / ((muu+1) * (1  + sqrt((aa+1)/(aa-1))*sqrt((muu-1)/(muu+1))))) 
  #(Wolfram) muu = sqrt(1 + (aa^2-1)/ratio^2)
  denb1 <- c * (a*c - mu*c + c*c - a*mu - bb)
  b1 <- (a*mu*(c-mu)*(a+c) - bb2 + c*c + bb*(c*(a-mu+c) - 2*a*mu))/denb1
  denb2 <- c * (a*c - mu*c - c*c + a*mu + bb)
  b2 <- (a*mu*(c+mu)*(a-c) + bb2 - c*c + bb*(c*(a-mu-c) + 2*a*mu))/denb2
  omegaT <- (b1 + b2)/2
  OmegaT <- c(omegaT, 0, 0)
  tormesh <- cgalMeshes:::torusMesh2(R, r, nu = nu, nv = nv)
  xvertices <- tormesh[["vb"]][1L:3L, ] + OmegaT
  for(i in 1L:nu){
    k0 <- i * nv - nv
    for(j in 1L:nv){
      k <- k0 + j
      vertices[, k] <- inversion(xvertices[, k])
    }
  }
  mesh <- tmesh3d(
    vertices    = vertices,
    indices     = tormesh[["it"]],
    homogeneous = FALSE
  )
  if(normals) {
    mesh <- addNormals(mesh)
  }
  mesh
}

n <- 512L
mesh <- conformalCyclideMesh(
  a = 97, c = 32, aspectRatio = 1, nu = n, nv = n, normals = FALSE
)

library(jacobi)
library(RcppColors)
# generate the values of wp on a grid
x <- y <- seq(-4, 0, length.out = n)
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, omega = c(0.5, 0.5 + 0.5i))
}
Z <- outer(x, y, f)
# map them to colors
img <- colorMap5(Z)
clrs <- c(img)

mesh[["material"]] <- list("color" = clrs)
# plot
open3d(windowRect = 50 + c(0, 0, 512, 512))
clear3d(type = "lights") # destroy current lights
light3d(x = -50, y = 100, z = 100)
bg3d("#363940")
view3d(-25, -25, zoom = 0.75)
shade3d(mesh, specular = "gold")
