library(cgalMeshes)
setwd("C:/SL/MyPackages/cgalMeshes/inst/trash")
library(rgl)

f <- function(x, y, z) { # Enneper surface: f=0
  64*z^9 - 128*z^7 + 64*z^5 - 702*x^2*y^2*z^3 - 18*x^2*y^2*z + 
    144*(y^2*z^6-x^2*z^6) + 162*(y^4*z^2-x^4*z^2) + 27*(y^6-x^6) +
    9*(x^4*z+y^4*z) + 48*(x^2*z^3+y^2*z^3) - 432*(x^2*z^5+y^2*z^5) +
    81*(x^4*y^2-x^2*y^4) + 240*(y^2*z^4-x^2*z^4) - 135*(x^4*z^3+y^4*z^3)
}

smesh <- sphereMesh(r = 0.5, iterations = 5L)
rmesh1 <- clipMesh3d(smesh, f, greater = TRUE, minVertices = 20000L)
rmesh2 <- clipMesh3d(smesh, f, greater = FALSE, minVertices = 20000L)

# convert to CGAL mesh
mesh1 <- cgalMesh$new(rmesh1)
mesh2 <- cgalMesh$new(rmesh2)

# # add vertices in order that the checkeboard has regular lines
mesh1$isotropicRemeshing(0.005, iterations = 3, relaxSteps = 2)
mesh2$isotropicRemeshing(0.005, iterations = 3, relaxSteps = 2)



#mesh$writeMeshFile("torus.off")
#M <- cgalMeshes:::testparam(normalizePath("torus.off"), 1L)

# compute mesh1 parameterization ####
UV <- mesh1$parameterization(method = "DCP", UVborder = "circle")

# radial checkerboard ####
UV0 <- UV
UV <- 10 * (UV0 - 0.5)
radii <- sqrt(apply(UV, 1L, crossprod))
angles <- 10 * (1 + atan2(UV[, 2L], UV[, 1L])/pi)
clrs1 <- ifelse(
  floor(radii) %% 2 == 0,
  ifelse(
    floor(angles) %% 2 == 0, "navy", "yellow"
  ),
  ifelse(
    floor(angles) %% 2 == 0, "yellow", "navy"
  )
)

# compute mesh2 parameterization ####
UV <- mesh2$parameterization(method = "DCP", UVborder = "circle")

# radial checkerboard ####
UV0 <- UV
UV <- 10 * (UV0 - 0.5)
radii <- sqrt(apply(UV, 1L, crossprod))
angles <- 10 * (1 + atan2(UV[, 2L], UV[, 1L])/pi)
clrs2 <- ifelse(
  floor(radii) %% 2 == 0,
  ifelse(
    floor(angles) %% 2 == 0, "darkviolet", "orangered"
  ),
  ifelse(
    floor(angles) %% 2 == 0, "orangered", "darkviolet"
  )
)


# compute normals, convert to 'rgl' mesh, and add colors ####
mesh1$computeNormals()
rmesh1 <- mesh1$getMesh()
rmesh1$material <- list(color = clrs1)

mesh2$computeNormals()
rmesh2 <- mesh2$getMesh()
rmesh2$material <- list(color = clrs2)

# plot ####
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, -20, zoom = 0.7)
shade3d(rmesh1, meshColor = "vertices")
shade3d(rmesh2, meshColor = "vertices")



snapshot3d(sprintf("MeeksMobius_checkerboard.png"), webshot = FALSE)




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


