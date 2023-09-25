library(cgalMeshes)
setwd("C:/SL/MyPackages/cgalMeshes/inst/trash")
library(rgl)

# Enneper surface parameteritazion ####
n <- 3
Enneper <- function(phi, r) {
  rbind(
    r*cos(phi) - r^(2*n-1)*cos((2*n-1)*phi)/(2*n-1),
    r*sin(phi) + r^(2*n-1)*sin((2*n-1)*phi)/(2*n-1),
    2*r^n*cos(n*phi)/n
  )
}

# do the mesh ####
rmesh <- parametricMesh(
  Enneper, urange = c(0, 2*pi), vrange = c(0, 1.3),
  periodic = c(TRUE, FALSE), nu = 512, nv = 512, clean = TRUE
)

# convert to CGAL mesh
mesh <- cgalMesh$new(rmesh)

# add vertices in order that the checkeboard has regular lines
mesh$isotropicRemeshing(0.01, iterations = 3, relaxSteps = 2)


#mesh$writeMeshFile("torus.off")
#M <- cgalMeshes:::testparam(normalizePath("torus.off"), 1L)

# compute mesh parameterization ####
UV <- mesh$parameterization(method = "DCP", UVborder = "circle")

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

# check the checkboard is correct ####
plot(
  UV0, type = "p", asp = 1, pch = ".", col = clrs, 
  xlab = "u", ylab = "v", xlim = c(0,1), ylim = c(0,1)
)

svg("x.svg")
opar <- par(mar=c(3,3,0,0))
plot(
  UV0, type = "p", asp = 1, pch = ".", col = clrs, 
  xlab = "u", ylab = "v", xlim = c(0,1), ylim = c(0,1)
)
par(opar)
dev.off()
rsvg::rsvg_png("x.svg", "checkerboard_flat.png", width = 512, height = 512)


# compute normals, convert to 'rgl' mesh, and add colors ####
mesh$computeNormals()
rmesh <- mesh$getMesh()
rmesh$material <- list(color = clrs)

# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, -20, zoom = 0.7)
shade3d(rmesh, meshColor = "vertices")
snapshot3d(sprintf("Enneper_checkerboard.png"), webshot = FALSE)

# animation rotating checkboard ####
angles0 <- angles
fclrs <- function(alpha) {
  tests <- floor(angles0 + alpha) %% 2 == 0
  ifelse(
    floor(radii) %% 2 == 0,
    ifelse(
      tests, "navy", "yellow"
    ),
    ifelse(
      tests, "yellow", "navy"
    )
  )
}

# make animation frames ####
alpha_ <- seq(0, 2, length.out = 19L)[-1L]
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, -20, zoom = 0.7)
for(i in seq_along(alpha_)) {
  clrs <- fclrs(alpha_[i])
  rmesh$material <- list(color = clrs)
  shade3d(rmesh, meshColor = "vertices")
  snapshot3d(sprintf("zzpic%03d.png", i), webshot = FALSE)
  clear3d()
}


# mount animation ####
library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "Enneper-radialCheckerboard-rotating.gif",
  width = 512,
  height = 512,
  delay = 1/8
)


file.remove(Sys.glob("zzpic*.png"))

appendGIFs("Enneper-nopsi-DCP.gif", "Enneper-psi-DCP.gif", delay = 8)
