#' @title Advancing front surface reconstruction
#' @description Reconstruction of a surface from a cloud of 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
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
#' library(rgl)
#' shade3d(rglMesh, color = "firebrick")
AFSreconstruction <- function(
    points
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
  xptr <- AFSreconstruction_cpp(t(points))
  cgalMesh$new(clean = xptr)
}
