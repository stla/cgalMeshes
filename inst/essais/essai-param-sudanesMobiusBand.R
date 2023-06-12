library(rgl)
library(cgalMeshes)
library(onion)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

oooo <- as.quaternion(c(1, 1, 1, 1), single = TRUE)/2

f <- Vectorize(function(theta, phi, tau = 1/2) {
  q <- as.quaternion(c(
    cos(theta)*cos(phi), 
    cos(theta)*sin(phi), 
    sin(theta)*cos(tau*phi), 
    sin(theta)*sin(tau*phi)
  ), single = TRUE)
  x <- as.numeric(oooo * q)
  x[1L:3L] / (1 - x[4L])
})

rmesh <- parametricMesh(f, c(-pi/2, pi/2), c(0, 2*pi), nu = 256, nv = 256, clean = TRUE, periodic = c(FALSE, FALSE))
#rmesh <- Rvcg::vcgClean(rmesh, sel=6, tol = 1e-6)
#Rvcg::checkFaceOrientation(rmesh)
shade3d(addNormals(rmesh), color = "red")
shade3d(
  getBoundary3d(
    rmesh, sorted = TRUE, color = "black", 
    lwd = 3, line_antialias = TRUE
  )
)

u <- v <- seq(0, 1, length.out = 257L)[-257L]
Grid <- expand.grid(U = u, V = v)
checkerboard <- ifelse(
  (floor(10*Grid$U) %% 2) == (floor(10*Grid$V) %% 2), 
  "yellow", "navy"
)
rmesh$material <- list(color = checkerboard)




mesh <- cgalMesh$new(rmesh)

# look at the edge lengths
summary(mesh$getEdges()[["length"]])

# do an isotropic remeshing to get smaller faces
mesh$isotropicRemeshing(1e-2, iterations = 3L, relaxSteps = 2L)
summary(mesh$getEdges()[["length"]])

# compute the discrete conformal parameterization
UV <- mesh$parameterization("DCP", UVborder = "circle")
head(UV)

# make a checkerboard with these points
UVcheckerboard <- ifelse(
  (floor(10*UV[, 1L]) %% 2) == (floor(10*UV[, 2L]) %% 2), 
  "yellow", "navy"
)

# compute mesh normals and convert to rgl mesh
mesh$computeNormals()
rmesh <- mesh$getMesh()

# assign these colors to the mesh and plot
rmesh$material <- list(color = UVcheckerboard)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-10, -35, zoom = 0.9)
shade3d(rmesh, meshColor = "vertices")


# animation ####
movie3d(spin3d(axis = c(0, 0, 1), rpm = 10),
        duration = 6, fps = 10,
        movie = "zzpic", dir = ".",
        convert = FALSE, webshot = FALSE,
        startTime = 1/10)

library(gifski)
gifski(
  png_files = Sys.glob("zzpic*.png"),
  gif_file = "sudanesMobiusBand.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))

