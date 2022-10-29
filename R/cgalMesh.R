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
    #' @importFrom rgl tmesh3d qmesh3d
    #' @examples 
    #' library(rgl)
    #' cgalMesh$new(cube3d())
    "initialize" = function(mesh, vertices, faces, clean = TRUE){
      if(inherits(clean, "externalptr")) {
        private[[".meshXPtr"]] <- CGALmesh$new(clean)
        return(invisible(self))
      }
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
    },
    
    #' @description Print a \code{cgalMesh} object.
    #' @param ... ignored
    #' @return No value returned, just prints some information about the mesh.
    "print" = function(...) {
      private[[".meshXPtr"]]$print()
    },
    
    #' @description Check whether the mesh bounds a volume. The mesh must be 
    #'   triangle.
    #' @return A Boolean value, whether the mesh bounds a volume.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(tetrahedron3d())
    #' mesh$boundsVolume() # TRUE
    #' mesh$reverseOrientation()
    #' mesh$boundsVolume() # TRUE
    "boundsVolume" = function() {
      if(!self$isTriangle()) {
        stop("The mesh is not triangle.")
      }
      private[[".meshXPtr"]]$doesBoundVolume()
    },
    
    #' @description Centroid of the mesh. The mesh must be triangle.
    #' @return The Cartesian coordinates of the centroid of the mesh.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(icosahedron3d())
    #' mesh$centroid()
    "centroid" = function() {
      if(!self$isTriangle()) {
        stop("The mesh is not triangle.")
      }
      private[[".meshXPtr"]]$centroid()
    },

    #' @description Copy the mesh.
    #' @return A new \code{cgalMesh} object.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())
    #' tmesh <- mesh$copy()$triangulate()
    #' tmesh$isTriangle() # TRUE
    #' mesh$isTriangle() # FALSE
    "copy" = function() {
      xptr <- private[[".meshXPtr"]]$clone()
      cgalMesh$new(clean = xptr)
    },

    #' @description Get the edges of the mesh.
    #' @return xxxxxx
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(dodecahedron3d())
    #' mesh$edges()
    "edges" = function() {
      private[[".meshXPtr"]]$edges()
    },
    
    #' @description Get the mesh.
    #' @param normals Boolean, whether to return the per-vertex normals 
    #' @param rgl Boolean, whether to return a \strong{rgl} mesh if possible, 
    #'   i.e. if the mesh only has triangular or quadrilateral faces
    #' @param ... arguments passed to \code{\link[rgl:mesh3d]{mesh3d}} (if 
    #'   a \strong{rgl} mesh is returned)
    #' @return xxxxxx
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())$triangulate()
    #' mesh$getMesh(FALSE)
    "getMesh" = function(normals = TRUE, rgl = TRUE, ...) {
      stopifnot(isBoolean(normals))
      stopifnot(isBoolean(rgl))
      mesh <- private[[".meshXPtr"]]$getRmesh(normals)
      if(rgl) {
        if(is.matrix(mesh[["faces"]])) {
          nsides <- nrow(mesh[["faces"]])
          if(nsides == 3L) {
            mesh <- mesh3d(
              x         = mesh[["vertices"]],
              triangles = mesh[["faces"]],
              normals   = mesh[["normals"]],
              ...
            )
          } else {
            mesh <- mesh3d(
              x       = mesh[["vertices"]],
              quads   = mesh[["faces"]],
              normals = mesh[["normals"]],
              ...
            )
          }
        } else {
          faces <- split(mesh[["faces"]], lengths(mesh[["faces"]]))
          if(all(names(faces) %in% c("3", "4"))) {
            mesh <- mesh3d(
              x         = mesh[["vertices"]],
              normals   = mesh[["normals"]],
              triangles = do.call(cbind, faces[["3"]]),
              quads     = do.call(cbind, faces[["4"]]),
              ...
            )
          } else {
            warning("Cannot make a rgl mesh.")
            mesh[["vertices"]] <- t(mesh[["vertices"]])
            mesh[["faces"]] <- t(mesh[["faces"]])
            if(normals) {
              mesh[["normals"]] <- t(mesh[["normals"]])
            }
          }
        }
      } else {
        mesh[["vertices"]] <- t(mesh[["vertices"]])
        mesh[["faces"]] <- t(mesh[["faces"]])
        if(normals) {
          mesh[["normals"]] <- t(mesh[["normals"]])
        }
      }
      mesh
    },

    #' @description Check whether the mesh is closed.
    #' @return A Boolean value, whether the mesh is closed.
    "isClosed" = function() {
      private[[".meshXPtr"]]$isClosed()
    },
    
    #' @description Check whether the mesh is triangle.
    #' @return A Boolean value, whether the mesh is triangle.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())
    #' mesh$isTriangle()
    "isTriangle" = function() {
      private[[".meshXPtr"]]$isTriangle()
    },

    #' @description Reverse the orientation of the faces of the mesh.
    #' @return The modified \code{cgalMesh} object. \strong{WARNING}: even if 
    #'   you store the result in a new variable, the original mesh is modified. 
    "reverseOrientation" = function() {
      private[[".meshXPtr"]]$reverseFaceOrientations()
      invisible(self)
    },
    
    #' @description Check whether the mesh self-intersects.
    #' @return A Boolean value, whether the mesh self-intersects.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())
    #' mesh$selfIntersects()
    "selfIntersects" = function() {
      private[[".meshXPtr"]]$doesSelfIntersect()
    },
    
    #' @description Triangulate mesh.
    #' @return The modified \code{cgalMesh} object. \strong{WARNING}: even if 
    #'   you store the result in a new variable, the original mesh is modified 
    #'   (see the example). You may want to triangulate a copy of the mesh; 
    #'   see the \code{copy} method.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())
    #' mesh$isTriangle() # FALSE
    #' # warning: triangulating the mesh modifies it
    #' x <- mesh$triangulate()
    #' mesh$isTriangle() # TRUE
    "triangulate" = function() {
      private[[".meshXPtr"]]$triangulate()
      invisible(self)
    }
    
  )
  
)