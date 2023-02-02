#' @title Advancing front surface reconstruction
#' @description Reconstruction of a surface from a cloud of 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param jetSmoothing if not \code{NULL}, must be an integer higher than two, 
#'   and then the points cloud is smoothed before the reconstruction, using 
#'   this integer as the number of neighbors for the smoothing; note that this 
#'   smoothing preprocessing relocates the points and then should not be used 
#'   if the points have been sampled without noise on the surface
#'
#' @return A \code{cgalMesh} object.
#'
#' @details See \href{https://doc.cgal.org/latest/Advancing_front_surface_reconstruction/index.html#Chapter_Advancing_Front_Surface_Reconstruction}{Advancing Front Surface Reconstruction}.
#'
#' @export
#'
#' @examples 
#' library(cgalMeshes)
#' data(bunny, package = "onion")
#' mesh <- AFSreconstruction(bunny)
#' rglMesh <- mesh$getMesh()
#' \donttest{library(rgl)
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' shade3d(rglMesh, color = "firebrick")}
#' 
#' # jet smoothing example ####
#' library(cgalMeshes)
#' # no smoothing
#' mesh1 <- AFSreconstruction(SolidMobiusStrip)
#' mesh1$computeNormals()
#' rglMesh1 <- mesh1$getMesh()
#' # jet smoothing
#' mesh2 <- AFSreconstruction(SolidMobiusStrip, jetSmoothing = 30)
#' mesh2$computeNormals()
#' rglMesh2 <- mesh2$getMesh()
#' # plot
#' \donttest{library(rgl)
#' open3d(windowRect = 50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(20, -40, zoom = 0.85)
#' shade3d(rglMesh1, color = "gold")
#' next3d()
#' view3d(20, -40, zoom = 0.85)
#' shade3d(rglMesh2, color = "gold")}
AFSreconstruction <- function(points, jetSmoothing = NULL){
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
  if(!is.null(jetSmoothing)) {
    stopifnot(isPositiveInteger(jetSmoothing), jetSmoothing >= 2L)
  } else {
    jetSmoothing <- 0L
  }
  xptr <- AFSreconstruction_cpp(t(points), as.integer(jetSmoothing))
  cgalMesh$new(clean = xptr)
}
