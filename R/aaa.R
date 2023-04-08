#' @useDynLib cgalMeshes, .registration=TRUE
#' @importFrom Rcpp evalCpp setRcppClass
#' @importFrom methods new
NULL

CGALmesh <- setRcppClass("CGALmesh")

CGALversion <- function() {
  5.4
}