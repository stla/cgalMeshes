library(rgl)
library(cgalMeshes)
setwd("~/Documents/R/MyPackages/cgalMeshes/inst/trash")

# soccer star
# isosurface f=0
a1 <- a2 <- a3 <- a4 <- -100
f <- function(x, y, z){
  u <- x*x + y*y + z*z
  v <- -z*(2*x+z)*(x^4 - x*x*z*z + z^4 + 2*(x^3*z-x*z^3) + 5*(y^4-y*y*z*z) + 10*(x*y*y*z-x*x*y*y))
  w <- (4*x*x + z*z - 6*x*z) * 
    (z^4 - 2*z^3*x - x*x*z*z + 2*z*x^3 + x^4 - 25*y*y*z*z - 30*x*y*y*z - 10*x*x*y*y + 5*y^4) *
    (z^4 + 8*z^3*x + 14*x*x*z*z - 8*z*x^3 + x^4 - 10*y*y*z*z - 10*x*x*y*y + 5*y^4)
  1 + ((128565+115200*sqrt(5))/1295029 * a3 + (49231296000*sqrt(5)-93078919125)/15386239549 * a4 - a1 - 3*a2 - 3) * u +
    ((-230400*sqrt(5) - 257130)/1295029 * a3 + (238926989250-126373248000*sqrt(5))/15386239549 * a4 + 3*a1 + 8*a2 + 3) * u * u + 
    ((115200*sqrt(5)+128565)/1295029 * a3 + (91097280000*sqrt(5)-172232645625)/15386239549 * a4 - 3*a1 - 6*a2 - 1) * u*u*u + 
    (a3 + (121075-51200*sqrt(5))/11881 * a4) * v + ((102400*sqrt(5)-242150)/11881 - 2*a3) * u * v + 
    a1 * u^4 + a2 * u^5 + a3 * u*u*v + a4*w
}
# make the isosurface
nx <- 500L; ny <- 400L; nz <- 200L
x <- seq(-1.5, 1.5, length.out = nx) 
y <- seq(-1.5, 1.5, length.out = ny)
z <- seq(0, 1.5, length.out = nz) 
Grid <- expand.grid(X = x, Y = y, Z = z)
voxel <- array(with(Grid, f(X, Y, Z)), dim = c(nx, ny, nz))
library(rmarchingcubes)
cont <- contour3d(voxel, level = 0, x = x, y = y, z = z)
# plot
library(rgl)
rmesh <- tmesh3d(
  vertices = t(cont[["vertices"]]),
  indices  = t(cont[["triangles"]]),
  normals  = cont[["normals"]],
  homogeneous = FALSE
)

view3d(0, 0)
shade3d(rmesh, color = "hotpink")

#
mesh <- cgalMesh$new(rmesh)
mesh$writeMeshFile("soccer.off")

# look at the edge lengths
summary(mesh$getEdges()[["length"]])
# do an isotropic remeshing to get smaller faces
mesh$isotropicRemeshing(0.005, iterations = 3L, relaxSteps = 2L)
summary(mesh$getEdges()[["length"]])
# parameterize
UV <- mesh$parameterization("DCP")

# make a radial checkerboard
UVnew <- UV - 0.5
radii <- sqrt(apply(10 * UVnew, 1L, crossprod))
angles <- 5 * (1 + atan2(UVnew[, 2L], UVnew[, 1L])/pi)
UVcheckerboard <- ifelse(
  floor(radii) %% 2 == 0,
  ifelse(
    floor(angles) %% 2 == 0, "navy", "yellow"
  ),
  ifelse(
    floor(angles) %% 2 == 0, "yellow", "navy"
  )
)

# compute mesh normals and convert to rgl mesh
mesh$computeNormals()
rmesh <- mesh$getMesh()

# assign these colors to the mesh and plot
rmesh$material <- list(color = UVcheckerboard)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, 0, zoom = 0.65)
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
  gif_file = "halfSoccerStar.gif",
  width = 512,
  height = 512,
  delay = 1/8
)
file.remove(Sys.glob("zzpic*.png"))


