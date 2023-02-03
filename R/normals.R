#' @title Normals for a points could
#' @description Returns a function which computes some normals for a 3D points
#'   cloud.
#'
#' @param nbNeighbors integer, number of neighbors used to compute the normals
#' @param method one of \code{"pca"} or \code{"jet"}
#'
#' @return A function which takes just one argument: a numeric matrix with
#'   three columns, each row represents a point, and the function returns a
#'   matrix of the same size as the input matrix, whose each row gives one
#'   unit normal for each point.
#'
#' @note The \code{getSomeNormals} function is intended to be used in the
#'   \code{\link{PoissonReconstruction}} function. If you want to use it for
#'   another purpose, be careful because the function it returns does not
#'   check the matrix it takes as argument.
#' @export
#'
#' @examples 
#' library(cgalMeshes)
#' library(rgl)
#' psr <- PoissonReconstruction( XXXXXXXXXXXXXXXXXXXXXXXXXXXX,
#'   SolidMobiusStrip, normals = "jet", neighbors = 6
#' )
#' open3d()
#' shade3d(psr, color = "cyan")
#' wire3d(psr)
getSomeNormals <- function(neighbors, method = "pca") {
  method <- match.arg(method, c("pca", "jet"))
  stopifnot(isPositiveInteger(neighbors))
  neighbors <- as.integer(neighbors)
  if(neighbors <= 2L) {
    stop("There must be at least two neighbors.", call. = TRUE)
  }
  if(method == "pca") {
    out <- function(points) t(pca_normals(t(points), neighbors))
  } else {
    out <- function(points) t(jet_normals(t(points), neighbors))
  }
  out
}
