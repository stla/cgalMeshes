convexHull <- function(points) {
  stopifnot(is.matrix(points))
  stopifnot(ncol(points) == 3L)
  storage.mode(points) <- "double"
  if(anyNA(points)) {
    stop("Found missing values in the `points` matrix.")
  }
  xptr <- cxhull(t(points))
  cgalMesh$new(clean = xptr)
}
