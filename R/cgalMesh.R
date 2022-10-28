#' @title R6 class to represent a CGAL mesh
#' @description R6 class to represent a CGAL mesh.
#'
#' @importFrom R6 R6Class
#' @export
cgalMesh <- R6Class(
  "cgalMesh",
  
  lock_class = TRUE,
  
  cloneable = FALSE,
  
  private = list(
    .meshXPtr = NULL
  ),
  
  public = list(
    
    #' @description Creates a new \code{cgalMesh} object.
    #' @param mesh either xxx
    #' @param vertices a numeric matrix with three columns
    #' @param faces either a matrix of integers (each row gives the vertex 
    #'   indices of a face) or a list of vectors of integers (each one gives 
    #'   the vertex indices of a face)
    #' @param clean Boolean, whether to clean the mesh (merge duplicated 
    #'   vertices and duplicated faces, remove isolated vertices); set to 
    #'   \code{FALSE} if you know your mesh is already clean
    #' @return A \code{cgalMesh} object.
    #' @examples 
    #' library(rgl)
    #' cgalMesh$new(cube3d())
    "initialize" = function(mesh, vertices, faces, clean = TRUE){
      stopifnot(isBoolean(clean))
      if(!missing(mesh)) {
        if(inherits(mesh, "mesh3d")) {
          mesh <- getVF(mesh)
        }
        VF <- checkMesh(mesh[["vertices"]], mesh[["faces"]], aslist = TRUE)
      } else {
        VF <- checkMesh(vertices, faces, aslist = TRUE)
      }
      private[[".meshXPtr"]] <- 
        CGALmesh$new(VF[["vertices"]], VF[["faces"]], clean)
      invisible(self)
    }
    
  )
  
)