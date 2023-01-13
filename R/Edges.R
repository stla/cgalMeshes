#' @title Exterior edges of a mesh
#' @description Returns the edges of a mesh whose corresponding dihedral angles 
#'   are not too flat.
#'
#' @param edgesDF the dataframe returned by the \code{edges} method of 
#'   \code{\link{cgalMesh}}
#' @param angleThreshold maximum deviation in degrees from the flat angle; for 
#'   example if \code{angleThreshold=1}, then an edge is considered as exterior 
#'   if its corresponding dihedral angle is lower than 179 or higher than 181
#'
#' @return An integer matrix giving the vertex indices of the exterior edges.
#' 
#' @note Once you get the exterior edges, say in \code{extEdges}, then you can 
#'   get the indices of the exterior vertices with 
#'   \code{which(table(extEdges) != 2)}.
#' 
#' @export
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' mesh <- cgalMesh$new(dodecahedron3d())
#' extEdges <- exteriorEdges(mesh$getEdges())
#' vertices <- mesh$getVertices()
#' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
#' shade3d(dodecahedron3d(), color = "tomato")
#' plotEdges(vertices, extEdges)
exteriorEdges <- function(edgesDF, angleThreshold = 1) {
  stopifnot(is.data.frame(edgesDF))
  stopifnot(isNonNegativeNumber(angleThreshold))
  angles <- edgesDF[["angle"]]
  keep <- abs(angles - 180) >= angleThreshold
  as.matrix(edgesDF[, c("i1", "i2")])[keep, ]
}

#' @title Plot some edges
#' @description Plot the given edges with \strong{rgl}.
#'
#' @param vertices a three-columns matrix giving the coordinates of the vertices
#' @param edges a two-columns integer matrix giving the edges by pairs of
#'   vertex indices
#' @param color a color for the edges
#' @param lwd line width, a positive number, ignored if \code{edgesAsTubes=TRUE}
#' @param edgesAsTubes Boolean, whether to draw the edges as tubes
#' @param tubesRadius the radius of the tubes when \code{edgesAsTubes=TRUE}
#' @param verticesAsSpheres Boolean, whether to draw the vertices as spheres
#' @param spheresRadius the radius of the spheres when
#'   \code{verticesAsSpheres=TRUE}
#' @param spheresColor the color of the spheres when
#'   \code{verticesAsSpheres=TRUE}
#'
#' @return No value, just produces a 3D graphic.
#'
#' @importFrom rgl cylinder3d shade3d segments3d spheres3d
#' @export
#'
#' @examples
#' library(cgalMeshes)
#' library(rgl)
#' 
#' mesh <- cgalMesh$new(pentagrammicPrism, clean = FALSE)
#' vertices <- mesh$getVertices()
#' edges <- mesh$getEdges()
#' extEdges <- exteriorEdges(edges)
#' tmesh <- mesh$triangulate()$getMesh()
#' 
#' \donttest{open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
#' shade3d(tmesh, color = "navy")
#' # we plot the exterior edges only
#' plotEdges(
#'   vertices, extEdges, color = "gold",
#'   tubesRadius = 0.02, spheresRadius = 0.02
#' )}
#' 
#' # or only plot the edges whose corresponding dihedral angle is acute:
#' sharpEdges <- as.matrix(subset(edges, angle <= 91, select = c("i1", "i2")))
#' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
#' shade3d(tmesh, color = "maroon")
#' plotEdges(
#'   vertices, sharpEdges, color = "darkred", 
#'   tubesRadius = 0.02, spheresRadius = 0.02
#' )
plotEdges <- function(
    vertices,
    edges,
    color = "black",
    lwd = 2,
    edgesAsTubes = TRUE,
    tubesRadius = 0.03,
    verticesAsSpheres = TRUE,
    spheresRadius = 0.05,
    spheresColor = color
){
  edges <- as.matrix(edges[, c(1L, 2L)])
  nedges <- nrow(edges)
  if(edgesAsTubes) {
    for(i in 1L:nedges){
      edge <- edges[i, ]
      tube <- cylinder3d(
        vertices[edge, ], radius = tubesRadius, sides = 90
      )
      shade3d(tube, color = color)
    }
  } else {
    pts <- matrix(NA_real_, nrow = 2L*nedges, ncol = 3L)
    is_v1 <- rep(c(TRUE, FALSE), times = nedges)
    pts[is_v1, ]  <- vertices[edges[, 1L], ]
    pts[!is_v1, ] <- vertices[edges[, 2L], ]
    segments3d(pts, color = color, lwd = lwd)
  }
  if(verticesAsSpheres){
    only <- unique(c(edges))
    vertices <- vertices[only, ]
    spheres3d(vertices, radius = spheresRadius, color = spheresColor)
  }
  invisible(NULL)
}
