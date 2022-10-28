#' @useDynLib cgalMeshes, .registration=TRUE
#' @importFrom Rcpp evalCpp setRcppClass
#' @importFrom methods new
NULL

CGALmesh <- Rcpp::setRcppClass("CGALmesh")

utils::globalVariables("exterior")
