#' @title Scale-space surface reconstruction
#' @description Reconstruction of a surface from a cloud of 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param scaleIterations number of iterations used to increase the scale
#' @param neighbors number of neighbors used to smooth the points cloud
#' @param samples number of samples used to smooth the points cloud
#' @param separateShells Boolean, whether to separate the shells
#' @param forceManifold Boolean, whether to force a manifold output mesh
#' @param borderAngle bound on the angle in degrees used to detect border edges
#'
#' @return A \code{cgalMesh} object.
#'
#' @details See \href{https://doc.cgal.org/latest/Scale_space_reconstruction_3/index.html}{Scale-space Surface Reconstruction}.
#'
#' @export
#'
#' @examples 
#' library(cgalMeshes)
#' mesh <- SSSreconstruction(
#'   SolidMobiusStrip, scaleIterations = 4, 
#'   forceManifold = TRUE, neighbors = 30
#' )
#' mesh$computeNormals()
#' rglMesh <- mesh$getMesh()
#' library(rgl)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(20, -40, zoom = 0.85)
#' shade3d(rglMesh, color = "tomato")
SSSreconstruction <- function(
  points, scaleIterations = 1, 
  neighbors = 12, samples = 300,
  separateShells = FALSE, forceManifold = TRUE, borderAngle = 45
){
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(points) != 3L){
    stop("The `points` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(points) <= 3L){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(anyNA(points)){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  stopifnot(isStrictPositiveInteger(scaleIterations))
  stopifnot(isPositiveInteger(neighbors), neighbors >= 2)
  stopifnot(isBoolean(separateShells))
  stopifnot(isBoolean(forceManifold))
  stopifnot(isNonNegativeNumber(borderAngle))
  xptr <- SSSreconstruction_cpp(
    t(points), as.integer(scaleIterations),
    as.integer(neighbors), as.integer(samples), 
    separateShells, forceManifold, as.double(borderAngle)
  )
  cgalMesh$new(clean = xptr)
}
