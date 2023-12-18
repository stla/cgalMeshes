setwd("C:/SL/MyPackages/cgalMeshes/inst/trash")

library(cgalMeshes)
library(rgl)
library(jacobi)

g2 <- (beta(1/4, 1/4)/2)^4
omega <- c(0.5, 0.5i)
f <- Vectorize(function(r, theta) {
  z <- r * exp(1i*theta)
  zeta <- wzeta(z, omega = omega)
  pprime <- wp(z, omega = omega, derivative = 1L)
  p <- wp(z, omega = omega)
  c(
    Re(pi*z - zeta - pi/g2*pprime),
    Im(pi*z + zeta - pi/g2*pprime),
    sqrt(6*pi/g2) * Re(p)
  )
})

rmesh0 <- parametricMesh(
  f, c(0.2, 0.8), c(-pi, pi), periodic = c(FALSE, TRUE), 
  nu = 256L, nv = 512L, clean = TRUE
)

# it self-intersects so we reconstruct it from its vertices
mesh <- alphaWrap(t(rmesh0$vb[-4, ]), 100, 60000)
#mesh <- AFSreconstruction(t(rmesh0$vb[-4, ]), jetSmoothing = 2)
rmesh <- mesh$getMesh()


rmesh2 <- Rvcg::vcgClean(rmesh, c(6,7), tol = 0.3)
shade3d(rmesh2, col = "deeppink4", back = "culled")
shade3d(rmesh2, col = "gold", front = "culled")
m <- cgalMesh$new(rmesh2, clean=TRUE)
rmesh <- m$getMesh()

# library(AlphaHull3D)
# ahull <- fullAhull3d(t(rmesh0$vb[-4, ]))
# rmesh <- setAlpha(ahull, alpha = 1)

shade3d(rmesh, col = "deeppink4", back = "culled")
shade3d(rmesh, col = "gold", front = "culled")
bbox3d()

rmesh2 <- Rvcg::vcgClean(rmesh, 0:5, tol = 1e-6)

# convert to CGAL mesh ####
mesh <- cgalMesh$new(rmesh)

# take a look at the edge lengths
edges <- mesh$getEdges()
summary(edges[["length"]])

# add vertices in order that the checkeboard has regular lines ####
mesh$isotropicRemeshing(0.06, iterations = 3, relaxSteps = 2)

# compute mesh parameterization ####
UV <- mesh$parameterization(method = "DCP", UVborder = "circle")

plot(UV, asp = 1, pch = ".", xlab = "u", ylab = "v", axes = TRUE)

# radial checkerboard ####
UV0 <- UV
UV <- 10 * (UV0 - 0.5)
radii <- sqrt(apply(UV, 1L, crossprod))
angles <- 10 * (1 + atan2(UV[, 2L], UV[, 1L])/pi)
clrs <- ifelse(
  floor(radii) %% 2 == 0,
  ifelse(
    floor(angles) %% 2 == 0, "navy", "yellow"
  ),
  ifelse(
    floor(angles) %% 2 == 0, "yellow", "navy"
  )
)

# check the checkerboard is correct ####
plot(
  UV0, type = "p", asp = 1, pch = ".", col = clrs, 
  xlab = "u", ylab = "v", xlim = c(0,1), ylim = c(0,1)
)

# add normals, convert to 'rgl' mesh, and add colors ####
mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh[["material"]] <- list("color" = clrs)
b <- getBoundary3d(rmesh, sorted = TRUE, color = "black")

# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
bg3d("#363940")
shade3d(rmesh, meshColor = "vertices", polygon_offset = 1)

vs <- mesh$getVertices()
bdrs <- mesh$getBorders()
points3d(vs[bdrs[[1L]][, 2L], ], color = "red")
points3d(vs[bdrs[[2L]][, 2L], ], color = "green")


shade3d(b, lwd = 4)

movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 20,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/20)

#snapshot3d(sprintf("EnneperWithCheckerboard_ARAP.png"), webshot = FALSE)

# mount animation
library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "EnneperCheckerboard_ARAP.gif",
  width = 512,
  height = 512,
  delay = 1/10
)
file.remove(Sys.glob("zzpic*.png"))


