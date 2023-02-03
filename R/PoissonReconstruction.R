#' @title Poisson surface reconstruction
#' @description Poisson reconstruction of a surface, from a cloud of 3D points.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param normals either a numeric matrix which stores the normals, one normal 
#'   per row (i.e. it must have the same size as the \code{points} matrix), or 
#'   one character string, either \code{"jet"} or \code{"pca"}; in this case, 
#'   the normals are computed with the help of the \code{\link{getSomeNormals}} 
#'   function, and you also have to specify the \code{neighbors} argument, 
#'   which is passed to \code{\link{getSomeNormals}}
#' @param neighbors argument passed to \code{\link{getSomeNormals}}, see above
#' @param spacing size parameter; smaller values increase the precision of the
#'   output mesh at the cost of higher computation time; set to \code{NULL}
#'   (the default) for a reasonable automatic value: an average spacing whose
#'   value will be displayed in a message and that you can also get in the
#'   \code{"spacing"} attribute of the output
#' @param sm_angle bound for the minimum facet angle in degrees
#' @param sm_radius relative bound for the radius of the surface Delaunay balls
#' @param sm_distance relative bound for the center-center distances
#'
#' @return A \code{cgalMesh} object.
#'
#' @details See \href{https://doc.cgal.org/latest/Poisson_surface_reconstruction_3/index.html}{Poisson Surface Reconstruction}.
#'
#' @export
#' 
#' @examples 
#' library(cgalMeshes)
#' library(rgl)
#' mesh <- PoissonReconstruction(SolidMobiusStrip)
#' mesh$computeNormals()
#' rmesh <- mesh$getMesh()
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(20, -40, zoom = 0.85)
#' shade3d(rmesh, color = "darkorange")
#' wire3d(rmesh)
PoissonReconstruction <- function(
    points, normals = "jet", neighbors = 12, spacing = NULL,
    sm_angle = 20, sm_radius = 30, sm_distance = 0.375
){
  if(!is.matrix(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  if(anyNA(points)){
    stop("Missing values in the `points` matrix are not allowed.", call. = TRUE)
  }
  if(ncol(points) != 3L){
    stop("The `points` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(points) <= 3L){
    stop("Insufficient number of points.", call. = TRUE)
  }
  if(!is.matrix(normals)) {
    if(!isString(normals) && !is.element(normals, c("jet", "pca"))) {
      stop("Invalid `normals` argument.")
    }
    normals <- getSomeNormals(neighbors, method = normals)(points)
  } else {
    if(ncol(normals) != 3L){
      stop("The `normals` matrix must have three columns.", call. = TRUE)
    }
    storage.mode(normals) <- "double"
    if(anyNA(normals)){
      stop("Missing values in the `normals` matrix are not allowed.")
    }
    if(nrow(points) != nrow(normals)){
      stop(
        "The `points` matrix and the `normals` matrix must have the same ",
        "number of rows.", call. = TRUE
      )
    }
  }
  if(is.null(spacing)) {
    spacing <- -1
  } else {
    stopifnot(isPositiveNumber(spacing))
  }
  stopifnot(isPositiveNumber(sm_angle))
  stopifnot(isPositiveNumber(sm_radius))
  stopifnot(isPositiveNumber(sm_distance))
  Psr <- Poisson_reconstruction_cpp(
    t(points), t(normals), spacing, sm_angle, sm_radius, sm_distance
  )
  out <- cgalMesh$new(clean = Psr[["xptr"]])
  if(spacing == -1){
    message(sprintf(
      "Poisson reconstruction using average spacing: %s.",
      formatC(Psr[["spacing"]])
    ))
    attr(out, "spacing") <- Psr[["spacing"]]
  }
  out
}
