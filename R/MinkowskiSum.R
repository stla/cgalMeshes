#' @title Minkowski addition
#' @description Minkowski addition of two meshes.
#'
#' @param mesh1,mesh2 two \code{cgalMesh} objects representing triangle 
#'   meshes
#'
#' @return A \code{cgalMesh} object. The mesh it represents is not triangle 
#'   in general.
#' @export
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' # cube + icosahedron
#' mesh1 <- cgalMesh$new(cube3d())$triangulate()
#' mesh2 <- cgalMesh$new(icosahedron3d())
#' mesh <- MinkowskiSum(mesh1, mesh2)
#' # get the edges before triangulation for plotting
#' edges <- mesh$getEdges()
#' # triangulation
#' mesh$triangulate()
#' # plot
#' rmesh <- mesh$getMesh()
#' open3d(windowRect = 50 + c(0, 0, 512, 512))
#' view3d(20, 40, zoom = 0.85)
#' shade3d(rmesh, color = "turquoise")
#' plotEdges(
#'   mesh$getVertices(), edges[, c("i1", "i2")], 
#'   color = "darkred", tubesRadius = 0.07, spheresRadius = 0.1
#' )
MinkowskiSum <- function(mesh1, mesh2) {
  stopifnot(isCGALmesh(mesh1))
  stopifnot(isCGALmesh(mesh2))
  xptr1 <- getXPtr(mesh1)
  xptr2 <- getXPtr(mesh2)
  xptr  <- MinkowskiSum_cpp(xptr1, xptr2)
  cgalMesh$new(clean = xptr)
}