setwd("C:/SL/MyPackages/cgalMeshes/inst/trash")

library(cgalMeshes)
library(rgl)

n <- 2
Enneper <- function(phi, r) {
  rbind(
    r*cos(phi) - r^(2*n-1)*cos((2*n-1)*phi)/(2*n-1),
    r*sin(phi) + r^(2*n-1)*sin((2*n-1)*phi)/(2*n-1),
    2*r^n*cos(n*phi)/n
  )
}

rmesh <- parametricMesh(
  Enneper, urange = c(0, 2*pi), vrange = c(0, 1.2),
  periodic = c(TRUE, FALSE), nu = 512L, nv = 128L, clean = TRUE
)

shade3d(rmesh, col = "deeppink4")
bbox3d()

# convert to CGAL mesh ####
mesh <- cgalMesh$new(rmesh)

# take a look at the edge lengths
edges <- mesh$getEdges()
summary(edges[["length"]])

# add vertices in order that the checkeboard has regular lines ####
mesh$isotropicRemeshing(0.008, iterations = 3, relaxSteps = 2)

# compute mesh parameterization ####
UV <- mesh$parameterization(method = "ARAP", lambda = 1000)

png("EnneperARAP_UVspace.png", width = 384, height = 384)
opar <- par(mar = c(4, 4, 0, 0))
plot(UV, asp = 1, pch = 19, xlab = "u", ylab = "v", axes = FALSE)
axis(1, at = seq(-3.5, 1.5, by = 0.5))
axis(2, at = seq(-1.5, 2.5, by = 0.5))
par(opar)
dev.off()


vs <- mesh$getVertices()
vsz <- vs[, 3L]
odec <- order(vsz, decreasing = TRUE)
vz1 <- odec[1L]
vz2 <- odec[2L]

uv1 <- UV[vz1, ]
uv2 <- UV[vz2, ]
alpha <- atan2(uv2[1L]-uv1[1L], uv2[2L]-uv1[2L])

rotation <- function(alpha, uv) {
  t(rbind(
    c(cos(alpha), -sin(alpha)),
    c(sin(alpha),  cos(alpha))
  ) %*% t(uv))
}

UVrot <- rotation(alpha, UV)
Urot <- UVrot[, 1L]
Vrot <- UVrot[, 2L]

png("EnneperARAP_UVspaceWithTheTwoExtremeZPoints.png", width = 384, height = 384)
opar <- par(mar = c(4, 4, 0, 0))
plot(UV, asp = 1, pch = ".", xlab = "u", ylab = "v", axes = FALSE)
axis(1, at = seq(-3.5, 1.5, by = 0.5))
axis(2, at = seq(-1.5, 2.5, by = 0.5))
points(UV[c(v1,v2),], col = c("red", "green"), pch = 19, cex = 2)
par(opar)
dev.off()


# ARAP checkerboard ####
# U <- UV[, 1L]
# V <- UV[, 2L]
# Un <- (U - min(U)) / (max(U) - min(U))
# Vn <- (V - min(V)) / (max(V) - min(V))
# checkerboard <- ifelse(
#   (floor(5 * U) %% 2) == (floor(5 * V) %% 2), 
#   "yellow", "navy"
# )

Un <- (Urot - min(Urot)) / (max(Urot) - min(Urot))
Vn <- (Vrot - min(Vrot)) / (max(Vrot) - min(Vrot))
checkerboard <- ifelse( # c'est ça le bon
  (floor(5 * Un) %% 2) == (floor(5 * Vn) %% 2), 
  "yellow", "navy"
)

png("Enneper_UVcheckerboard_ARAP.png", width = 384, height = 384)
opar <- par(mar = c(4, 4, 0, 0))
plot(
  U, V, asp = 1, pch = ".", xlab = "u", ylab = "v", 
  axes = FALSE, col = checkerboard
)
axis(1, at = seq(-3.5, 1.5, by = 0.5))
axis(2, at = seq(-1.5, 2.5, by = 0.5))
par(opar)
dev.off()

# add normals, convert to 'rgl' mesh, and add colors ####
mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh[["material"]] <- list("color" = checkerboard)
b <- getBoundary3d(rmesh, sorted = TRUE, color = "black")

# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
bg3d("#363940")
par3d(userMatrix = um)
shade3d(rmesh, meshColor = "vertices", polygon_offset = 1)
shade3d(b, lwd = 4)

snapshot3d(sprintf("EnneperWithCheckerboard_ARAP.png"), webshot = FALSE)
saveRDS(um, "userMatrix.rds")

# animation ####
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
  movie = "zzpic",
  convert = FALSE, webshot = FALSE
)

# mount animation
library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "TennisBallWithCheckerboard_FourCorners.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))


