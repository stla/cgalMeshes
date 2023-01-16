#' @title Mesh of an algebraic surface
#' @description Computes a mesh of an algebraic surface (isosurface defined 
#'   by a polynomial).
#'
#' @param polynomial the polynomial defining the isosurface; this can be either 
#'   a \code{\link[spray]{spray}} object or a list with two fields: 
#'   \code{exponents}, an integer matrix with three columns, and \code{coeffs}, 
#'   a numeric vector
#' @param isolevel a number, the value of the polynomial defining the 
#'   isosurface: the isosurface is defined by \code{P(x,y,z)=isolevel}, 
#'   where \code{P} is the polynomial 
#' @param sphereCenter,sphereRadius center and radius of a sphere bounding the 
#'   isosurface; the value of the polynomial at the center of 
#'   this sphere must be less than the isolevel 
#' @param angleBound lower bound in degrees for the angles of the faces of 
#'   the mesh
#' @param radiusBound upper bound for the radii of surface Delaunay balls; 
#'   a surface Delaunay ball is a ball circumscribing a face and centered 
#'  at the surface
#' @param distanceBound upper bound for the distance between the 
#'   circumcenter of a face and the center of the surface Delaunay 
#'   ball of this face
#'
#' @return A \code{cgalMesh} object. The mesh has normals computed 
#'   with the gradient of the polynomial.
#'   
#' @note If either \code{angleBound} is too high, \code{radiusBound} is 
#'   too small, or \code{distanceBound} is too small, the computation could 
#'   never finishes. Therefore, it is recommended to start with a small 
#'   \code{angleBound} and large \code{radiusBound} and \code{distanceBound}, 
#'   and then to adjust these parameters if the resulting mesh is not nice.  
#'   
#' @export
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' # the Barth decic
#' phi <- (1 + sqrt(5)) / 2
#' f <- function(x, y, z) {
#'   (18+30*phi) * x**4 * y**4 +
#'     (-36 - 30*phi + 44*phi**2 - 10*phi**3) * x**2 +
#'     (-18 - 24*phi + 10*phi**2) * x**6 +
#'     (3 + 5*phi) * x**8 +
#'     (36 + 60*phi) * x**2 * y**2 * z**4 +
#'     (12 + 20*phi) * x**2 * y**6 +
#'     phi +
#'     (-16 + 8*phi**4 - 8*phi**8 + 16*phi**12) * x**2 * y**4 * z**4 +
#'     (8 * phi**8) * y**2 * z**8 +
#'     (-18 - 24*phi + 10*phi**2) * z**6 +
#'     (-8*phi**4 - 16*phi**8) * x**6 * z**4 +
#'     (16*phi**4 + 8*phi**8) * x**6 * y**4 +
#'     (-8*phi**4) * y**8 * z**2 +
#'     (-18 - 24*phi + 10*phi**2) * y**6 +
#'     (12 + 20*phi) * x**2 * z**6 +
#'     (36 + 60*phi) * x**4 * y**2 * z**2 +
#'     (36 + 60*phi) * x**2 * y**4 * z**2 +
#'     (8 + 16*phi**4 -16*phi**8 - 8*phi**12) * x**2 * y**2 * z**6 +
#'     (-54 - 72*phi + 30*phi**2) * y**4 * z**2 +
#'     (-8*phi**4) * x**8 * y**2 +
#'     (16*phi**4 + 8*phi**8) * y**6 * z**4 +
#'     (12 + 20*phi) * y**2 * z**6 +
#'     (3 + 5*phi) * z**8 +
#'     (-8*phi**4) * x**2 * z**8 +
#'     (39 + 41*phi - 37*phi**2 + 5*phi**3) * z**4 +
#'     (-54 - 72*phi + 30*phi**2) * x**2 * y**4 +
#'     (8 + 16*phi**4 -16*phi**8 - 8*phi**12) * x**6 * y**2 * z**2 +
#'     (-54 - 72*phi + 30*phi**2) * x**2 * z**4 +
#'     (12 + 20*phi) * x**6 * z**2 +
#'     (-16 + 8*phi**4 - 8*phi**8 + 16*phi**12) * x**4 * y**2 * z**4 +
#'     (16*phi**4 + 8*phi**8) * x**4 * z**6 +
#'     (39 + 41*phi - 37*phi**2 + 5*phi**3) * y**4 +
#'     (-36 - 30*phi + 44*phi**2 - 10*phi**3) * z**2 +
#'     (8*phi**8) * x**2 * y**8 +
#'     (12 + 20*phi) * y**6 * z**2 +
#'     (8*phi**8) * x**8 * z**2 +
#'     (-36 - 30*phi + 44*phi**2 - 10*phi**3) * y**2 +
#'     (12 + 20*phi) * x**6 * y**2 +
#'     (-8*phi**4 - 16*phi**8) * y**4 * z**6 +
#'     (-16 + 8*phi**4 - 8*phi**8 + 16*phi**12) * x**4 * y**4 * z**2 +
#'     (78 + 82*phi - 74*phi**2 + 10*phi**3) * x**2 * z**2 +
#'     (18+30*phi) * x**4 * z**4 +
#'     (-8*phi**4 - 16*phi**8) * x**4 * y**6 +
#'     (-54 - 72*phi + 30*phi**2) * x**4 * y**2 +
#'     (-54 - 72*phi + 30*phi**2) * x**4 * z**2 +
#'     (-54 - 72*phi + 30*phi**2) * y**2 * z**4 +
#'     (78 + 82*phi - 74*phi**2 + 10*phi**3) * x**2 * y**2 +
#'     (-108 - 144*phi + 60*phi**2) * x**2 * y**2 * z**2 +
#'     (18+30*phi) * y**4 * z**4 +
#'     (3 + 5*phi) * y**8 +
#'     (78 + 82*phi - 74*phi**2 + 10*phi**3) * y**2 * z**2 +
#'     (8 + 16*phi**4 -16*phi**8 - 8*phi**12) * x**2 * y**6 * z**2 +
#'     (39 + 41*phi - 37*phi**2 + 5*phi**3) * x**4
#' }
#' # define f as a polynomial
#' library(spray)
#' x <- lone(1, 3)
#' y <- lone(2, 3)
#' z <- lone(3, 3)
#' P <- f(x, y, z)
#' # compute the mesh of the algebraic surface
#' # we have to negate P in order that the value of P at the 
#' # center of the bounding sphere is less than the isovalue
#' \donttest{mesh <- algebraicMesh(
#'   polynomial = -P, isolevel = 0,
#'   sphereCenter = c(0, 0, 0),
#'   sphereRadius = sqrt((5+sqrt(5))/2),
#'   angleBound    = 20, 
#'   radiusBound   = 0.015, 
#'   distanceBound = 0.015
#' )
#' # plot
#' rmesh <- mesh$getMesh()
#' open3d(windowRect = c(50, 50, 562, 562))
#' view3d(20, 20, zoom= 0.75)
#' bg3d(rgb(54, 57, 64, maxColorValue = 255))
#' shade3d(rmesh, color = "orangered")}
algebraicMesh <- function(
  polynomial, isolevel, 
  sphereCenter, sphereRadius,
  angleBound, radiusBound, distanceBound
) {
  stopifnot(isNumber(isolevel))
  stopifnot(isVector3(sphereCenter))
  stopifnot(isPositiveNumber(sphereRadius))
  stopifnot(isPositiveNumber(angleBound))
  stopifnot(isPositiveNumber(radiusBound))
  stopifnot(isPositiveNumber(distanceBound))
  if(inherits(polynomial, "spray")) {
    exponents <- polynomial[["index"]]
    coeffs    <- polynomial[["value"]]
  } else if(is.list(polynomial)) {
    exponents <- polynomial[["exponents"]]
    stopifnot(is.matrix(exponents))
    stopifnot(ncol(exponents) == 3L)
    storage.mode(exponents) <- "integer"
    if(anyNA(exponents)) {
      stop("Found missing values in the exponents of the polynomial.")
    }
    if(any(exponents < 0L)) {
      stop("The exponents of the polynomial must be positive.")
    }
    coeffs <- polynomial[["coeffs"]]
    storage.mode(coeffs) <- "double"
    if(anyNA(coeffs)) {
      stop("Found missing values in the coefficients of the polynomial.")
    }
    stopifnot(nrow(exponents) == length(coeffs))
  } else {
    stop("Invalid `polynomial` argument.")
  }
  xptr <- AlgebraicMesh(
    exponents, coeffs, isolevel, 
    sphereCenter, sphereRadius, 
    angleBound, radiusBound, distanceBound
  )
  cgalMesh$new(clean = xptr)
}

#' @title Extract mesh from a NifTI file.
#' @description Make a triangle mesh from a NifTI file and an isovalue.
#'
#' @param fileName path to the NifTI file (extension \code{nii} 
#'   or \code{nii.gz})
#' @param isolevel the level at which to construct the isosurface
#' @param sphereCenter,sphereRadius center and radius of a bounding sphere; 
#'   if \code{NULL} is given, a default bounding sphere is used, covering 
#'   all the space
#' @param angleBound lower bound in degrees for the angles of the faces of 
#'   the mesh
#' @param radiusBound upper bound for the radii of surface Delaunay balls; 
#'   a surface Delaunay ball is a ball circumscribing a face and centered 
#'  at the surface
#' @param distanceBound upper bound for the distance between the 
#'   circumcenter of a face and the center of the surface Delaunay 
#'   ball of this face
#'
#' @return A \code{cgalMesh} object.
#' 
#' @note If either \code{angleBound} is too high, \code{radiusBound} is 
#'   too small, or \code{distanceBound} is too small, the computation could 
#'   never finishes. Therefore, it is recommended to start with a small 
#'   \code{angleBound} and large \code{radiusBound} and \code{distanceBound}, 
#'   and then to adjust these parameters if the resulting mesh is not nice.
#'   
#' @export
#' @importFrom RNifti readNifti writeAnalyze niftiHeader
#' @examples 
#' library(cgalMeshes)
#' library(rgl)
#' teapot <- system.file("extdata", "teapot.nii.gz", package = "cgalMeshes")
#' mesh <- voxel2mesh(
#'   teapot, isolevel = 6,
#'   sphereCenter = c(128, 128, 89), 
#'   sphereRadius = 200,
#'   angleBound    = 20, 
#'   radiusBound   = 5, 
#'   distanceBound = 5
#' )
#' mesh$computeNormals()
#' rmesh <- rotate3d(mesh$getMesh(), pi, 1, 0, 0)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(0, 0)
#' shade3d(rmesh, color = "red")
#' wire3d(rmesh)
voxel2mesh <- function(
    fileName, isolevel,
    sphereCenter = NULL, sphereRadius = NULL,
    angleBound, radiusBound, distanceBound
) {
  stopifnot(isFilename(fileName))
  stopifnot(isNumber(isolevel))
  stopifnot(is.null(sphereCenter) || isVector3(sphereCenter))
  stopifnot(is.null(sphereRadius) || isPositiveNumber(sphereRadius))
  stopifnot(isPositiveNumber(angleBound))
  stopifnot(isPositiveNumber(radiusBound))
  stopifnot(isPositiveNumber(distanceBound))
  checkFile <- grepl("\\.nii\\.gz$", tolower(fileName)) || 
    grepl("\\.nii$", tolower(fileName))
  if(!checkFile) {
    stop("Not a NifTI file.")
  }
  if(is.null(sphereCenter) || is.null(sphereRadius)) {
    hdr <- niftiHeader(fileName)
    spacing <- attr(hdr, "pixdim")
    imgdim  <- spacing * attr(hdr, "imagedim")
    sphereCenter <- imgdim / 2
    sphereRadius <- sqrt(sum((imgdim - sphereCenter)^2))
  }
  hdrFile <- tempfile(fileext = ".hdr")
  img <- readNifti(fileName)
  . <- writeAnalyze(img, hdrFile)
  xptr <- VoxelToMesh(
    hdrFile, isolevel, 
    sphereCenter, 
    sphereRadius,
    angleBound, radiusBound, distanceBound
  )
  cgalMesh$new(clean = xptr)
}
