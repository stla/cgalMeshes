library(rgl)

# cross product
xprod <- function(v, w){
  c(
    v[2] * w[3] - v[3] * w[2], 
    v[3] * w[1] - v[1] * w[3], 
    v[1] * w[2] - v[2] * w[1]
  )
}

# knot curve
knot <- function(t, p, q){
  r <- cos(q*t)+2
  c(
    r * cos(p*t),
    r * sin(p*t),
    -sin(q*t)
  )
}
# derivative (tangent)
dknot <- function(t, p, q){
  v <- c(
    -q*sin(q*t)*cos(p*t) - p*sin(p*t)*(cos(q*t)+2),
    -q*sin(q*t)*sin(p*t) + p*cos(p*t)*(cos(q*t)+2),
    -q*cos(q*t)
  )
  v / sqrt(c(crossprod(v)))
}
# second derivative (normal)
ddknot <- function(t, p, q){
  v <- c(
    -q*(q*cos(q*t)*cos(p*t)-p*sin(p*t)*sin(q*t)) - 
      p*(p*cos(p*t)*cos(q*t+2)-q*sin(q*t)*sin(p*t)),
    -q*(q*cos(q*t)*sin(p*t)+p*cos(p*t)*sin(q*t)) + 
      p*(-p*sin(p*t)*cos(q*t+2)-q*sin(q*t)*cos(p*t)),
    q*q*sin(q*t)
  )
  v / sqrt(c(crossprod(v)))
}
# binormal
bnrml <- function(t, p, q){
  v <- xprod(dknot(t, p, q), ddknot(t, p, q))
  v / sqrt(c(crossprod(v)))
}

# mesh: tubular knot
TubularKnotMesh <- function(p, q, a, nu, nv){
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  vs <- matrix(NA_real_, nrow=3L, ncol=nu*nv)
  u_ <- seq(0, 2*pi, length.out = nu+1L)[-1L]
  v_ <- seq(0, 2*pi, length.out = nv+1L)[-1L]
  for(i in 1:nu){
    u <- u_[i]
    for(j in 1:nv){
      v <- v_[j]
      h <- knot(u, p, q)
      vs[,(i-1)*nv+j] <- 
        h + 
        a*(cos(v) * 
             (cos(u)*ddknot(u, p, q) + 
                sin(u)*bnrml(u, p, q)) + 
             sin(v) * 
             (-sin(u)*ddknot(u, p, q) + 
                cos(u)*bnrml(u, p, q)))
    }
  }
  tris1 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
  tris2 <- matrix(NA_integer_, nrow = 3L, ncol = nu*nv)
  for(i in 1L:nu){
    ip1 <- ifelse(i == nu, 1L, i+1L)
    for(j in 1L:nv){
      jp1 <- ifelse(j==nv, 1L, j+1L)
      tris1[,(i-1)*nv+j] <- c((i-1L)*nv+j,(i-1L)*nv+jp1, (ip1-1L)*nv+j)   
      tris2[,(i-1)*nv+j] <- c((i-1L)*nv+jp1,(ip1-1L)*nv+jp1,(ip1-1L)*nv+j)   
    }
  }
  out <- tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    homogeneous = FALSE
  )
  addNormals(out)
}

m <- TubularKnotMesh(p = 3, q = 11, a = 0.15, 10*60, 60)

open3d(windowRect = c(50, 50, 562, 562))
bg3d(rgb(54, 57, 64, maxColorValue = 255))
view3d(0, 0, zoom = 0.8)
shade3d(m, color = "magenta")

# anim
mesh <- cgalMesh$new(m)
rglmesh <- mesh$getMesh()
rglmesh0 <- rglmesh
vs <- mesh$vertices()
# estimated geodesic distances
index <- which.min(vs[, 1L])
geodists <- mesh$geoDists(index)
x <- seq(min(geodists), max(geodists), length.out = 90L)
# coloring function
fcolor <- colorRamp(viridisLite::magma(200L))
cols <- fcolor(geodists / max(geodists))
colors0 <- rgb(cols[, 1L], cols[, 2L], cols[, 3L], maxColorValue = 255)

for(i in seq_along(x)) {
  rglmesh <- rglmesh0
  colors <- colors0
  red <- geodists >= x[i]
  rglmesh[["vb"]][, which(red)] <- NA_real_
  rglmesh[["material"]] <- list("color" = colors)
  # plot
  open3d(windowRect = 50 + c(0, 0, 512, 512))
  view3d(0, 0, zoom = 0.8)
  shade3d(rglmesh0, alpha = 0)
  shade3d(rglmesh)
  snapshot3d(sprintf("zzpic%03d.png", i), webshot = FALSE)
  close3d()
}

for(k in 1:10){
  file.copy("zzpic090.png", sprintf("zzpic%03d.png", 90+k))
}

command <- "convert -delay 1x11 -duplicate 1,-2-1 -layers OptimizePlus zzpic*.png knot-3-11.gif"
system(command)
file.remove(Sys.glob("zzpic*.png"))

rglmesh <- m 
mesh <- cgalMesh$new(rglmesh)
# estimated geodesic distances
geodists <- mesh$geoDists(1L)
# normalization to (0, 1)
geodists <- (geodists - min(geodists)) / (max(geodists) - min(geodists))
# color each vertex according to its geodesic distance from the source
fcolor <- colorRamp(viridisLite::viridis(200L), bias = 2)
colors <- fcolor(geodists)
colors <- rgb(colors[, 1L], colors[, 2L], colors[, 3L], maxColorValue = 255)
rglmesh[["material"]] <- list("color" = colors)
# plot
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, 0, zoom = 0.8)
shade3d(rglmesh)
if(!rgl.useNULL()) {
  play3d(spin3d(axis = c(0, 1, 1), rpm = 5), duration = 20)  
}



# animation ####
movie3d(
  spin3d(axis = c(0, 0, 1), rpm = 20),
  duration = 3, fps = 20,
  movie = "zzpic", dir = ".",
  convert = FALSE,
  startTime = 1/20,
  webshot = FALSE
)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png")[1L:30L], # 120 / n'lobes'
  gif_file = "TwistedTubularKnot_p2_q8.gif",
  width = 512,
  height = 512,
  delay = 1/12
)

pngfiles <- Sys.glob("zzpic*.png")
file.remove(pngfiles)
