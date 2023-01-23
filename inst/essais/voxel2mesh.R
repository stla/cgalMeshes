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
#' @param errorBound a relative error bound used in the computations
#'
#' @return A \code{cgalMesh} object.
#' 
#' @note If either \code{angleBound} is too high, \code{radiusBound} is 
#'   too small, or \code{distanceBound} is too small, the computation could 
#'   never finishes. Therefore, it is recommended to start with a small 
#'   \code{angleBound} and large \code{radiusBound} and \code{distanceBound}, 
#'   and then to adjust these parameters if the resulting mesh is not nice.
#'   
#' @noRd
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
    angleBound, radiusBound, distanceBound,
    errorBound = 1e-3
) {
  stopifnot(isFilename(fileName))
  stopifnot(isNumber(isolevel))
  stopifnot(is.null(sphereCenter) || isVector3(sphereCenter))
  stopifnot(is.null(sphereRadius) || isPositiveNumber(sphereRadius))
  stopifnot(isPositiveNumber(angleBound))
  stopifnot(isPositiveNumber(radiusBound))
  stopifnot(isPositiveNumber(distanceBound))
  stopifnot(isPositiveNumber(errorBound))
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
    sphereCenter, sphereRadius,
    angleBound, radiusBound, distanceBound,
    errorBound
  )
  cgalMesh$new(clean = xptr)
}
