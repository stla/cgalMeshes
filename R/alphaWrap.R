#' @title 3D alpha wrapping
#' @description Reconstruction of a surface from a cloud of 3D points by 
#'   alpha wrapping.
#'
#' @param points numeric matrix which stores the points, one point per row
#' @param ralpha relative alpha parameter; the alpha parameter (see details) 
#'   is then defined as the length of the diagonal of the bounding box of the 
#'   points cloud divided by the relative alpha parameter
#' @param roffset relative offset; the offset parameter (see details) 
#'   is then defined as the length of the diagonal of the bounding box of the 
#'   points cloud divided by the relative offset parameter
#'
#' @return A \code{cgalMesh} object.
#'
#' @details See \href{https://doc.cgal.org/latest/Alpha_wrap_3/index.html}{3D Alpha Wrapping} 
#'   for details. The smallest alpha parameter, the smallest triangles in the 
#'   output mesh. The offset is the distance from the input points to the output 
#'   mesh.
#'
#' @export
#'
#' @examples 
#' library(cgalMeshes)
#' library(rgl)
#' # take the diplodocus mesh
#' off <- system.file("extdata", "diplodocus.off", package = "cgalMeshes")
#' diplodocusMesh <- cgalMesh$new(off)
#' diplodocusMesh$computeNormals()
#' # reconstruct the mesh from its vertices
#' pts <- diplodocusMesh$getVertices()
#' \donttest{wrapMesh <- alphaWrap(pts, 70, 3000)
#' wrapMesh$computeNormals()
#' # plot
#' diplodocusRglMesh <- diplodocusMesh$getMesh()
#' wrapRglMesh <- wrapMesh$getMesh()
#' open3d(windowRect = 50 + c(0, 0, 800, 400))
#' mfrow3d(1, 2)
#' view3d(20, 0, zoom = 0.85)
#' shade3d(diplodocusRglMesh, color = "forestgreen")
#' next3d()
#' view3d(20, 0, zoom = 0.85)
#' shade3d(wrapRglMesh, color = "forestgreen")}
alphaWrap <- function(points, ralpha, roffset){
  stopifnot(isPositiveNumber(ralpha))
  stopifnot(isPositiveNumber(roffset))
  if(!is.matrix(points) || !is.numeric(points)){
    stop("The `points` argument must be a numeric matrix.", call. = TRUE)
  }
  if(ncol(points) != 3L){
    stop("The `points` matrix must have three columns.", call. = TRUE)
  }
  if(nrow(points) <= 3L){
    stop("Insufficient number of points.", call. = TRUE)
  }
  storage.mode(points) <- "double"
  if(anyNA(points)){
    stop("Points with missing values are not allowed.", call. = TRUE)
  }
  xptr <- alphaWrap_cpp(t(points), ralpha, roffset)
  cgalMesh$new(clean = xptr)
}
