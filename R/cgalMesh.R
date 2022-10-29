#' @title R6 class to represent a CGAL mesh
#' @description R6 class to represent a CGAL mesh.
#'
#' @importFrom R6 R6Class
#' @importFrom rgl mesh3d
#' @export
cgalMesh <- R6Class(
  "cgalMesh",
  
  lock_class = TRUE,
  
  cloneable = FALSE,
  
  private = list(
    ".CGALmesh" = NULL
  ),
  
  public = list(
    
    #' @description Creates a new \code{cgalMesh} object.
    #' @param mesh there are four possibilities for this argument: it can be 
    #'   missing, in which case the arguments \code{vertices} and \code{faces} 
    #'   must be given, or it can be the path to a mesh file (accepted formats: 
    #'   \code{off}, \code{obj}, \code{stl}, \code{ply}, \code{ts}, \code{vtp}),
    #'   or it can be a \strong{rgl} mesh (i.e. a \code{mesh3d} object), or it 
    #'   can be a list containing (at least) the fields \code{vertices} 
    #'   (numeric matrix with three columns) and \code{faces} (matrix of 
    #'   integers or list of vectors of integers)
    #' @param vertices if \code{mesh} is missing, must be a numeric matrix with 
    #'   three columns
    #' @param faces if \code{mesh} is missing, must be either a matrix of 
    #'   integers (each row gives the vertex indices of a face) or a list of 
    #'   vectors of integers (each one gives the vertex indices of a face)
    #' @param clean Boolean, no effect if the mesh is given by a file, 
    #'   otherwise it indicates whether to clean the mesh (merge duplicated 
    #'   vertices and duplicated faces, remove isolated vertices); set to 
    #'   \code{FALSE} if you know your mesh is already clean
    #' @return A \code{cgalMesh} object.
    #' @importFrom rgl tmesh3d qmesh3d
    #' @examples 
    #' library(cgalMeshes)
    #' meshFile <- system.file(
    #'   "extdata", "bigPolyhedron.off", package = "cgalMeshes"
    #' )
    #' mesh <- cgalMesh$new(meshFile)
    #' rglmesh <- mesh$getMesh(normals = FALSE)
    #' library(rgl)
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglmesh, color = "tomato")
    #' \donttest{plotEdges(
    #'   mesh$vertices(), exteriorEdges(mesh$edges()), color = "darkred"
    #' )}
    "initialize" = function(
      mesh, vertices, faces, clean = TRUE
    ){
      if(inherits(clean, "externalptr")) {
        private[[".CGALmesh"]] <- CGALmesh$new(clean)
        return(invisible(self))
      }
      stopifnot(isBoolean(clean))
      if(!missing(mesh)) {
        if(inherits(mesh, "mesh3d")) {
          mesh <- getVF(mesh)
        }
        if(is.list(mesh)) {
          VF <- checkMesh(mesh[["vertices"]], mesh[["faces"]], aslist = TRUE)
        } else if(isFilename(mesh)) {
          private[[".CGALmesh"]] <- CGALmesh$new(path.expand(mesh), TRUE)
          return(invisible(self))
        } else {
          stop("Invalid `mesh` argument.")
        }
      } else {
        VF <- checkMesh(vertices, faces, aslist = TRUE)
      }
      private[[".CGALmesh"]] <- 
        CGALmesh$new(VF[["vertices"]], VF[["faces"]], clean)
      invisible(self)
    },
    
    #' @description Print a \code{cgalMesh} object.
    #' @param ... ignored
    #' @return No value returned, just prints some information about the mesh.
    "print" = function(...) {
      private[[".CGALmesh"]]$print()
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
      private[[".CGALmesh"]]$doesBoundVolume()
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
      private[[".CGALmesh"]]$centroid()
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
      xptr <- private[[".CGALmesh"]]$clone()
      cgalMesh$new(clean = xptr)
    },
    
    #' @description Get the edges of the mesh.
    #' @return A dataframe with four columns; the first two ones give the 
    #'   vertex indices of each edge (one edge per row), the third one gives 
    #'   the lengths of each edge, and the fourth one gives the dihedral angles 
    #'   in degrees between the two faces adjacent to each edge 
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(dodecahedron3d())
    #' mesh$edges()
    "edges" = function() {
      private[[".CGALmesh"]]$edges()
    },
    
    #' @description Get the mesh.
    #' @param normals Boolean, whether to return the per-vertex normals 
    #' @param rgl Boolean, whether to return a \strong{rgl} mesh if possible, 
    #'   i.e. if the mesh only has triangular or quadrilateral faces
    #' @param ... arguments passed to \code{\link[rgl:mesh3d]{mesh3d}} (if 
    #'   a \strong{rgl} mesh is returned)
    #' @return A \strong{rgl} mesh or a list with two or three fields: 
    #'   \code{vertices}, \code{faces}, and \code{normals} if the argument 
    #'   \code{normals} is set to \code{TRUE}
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())$triangulate()
    #' mesh$getMesh(normals = FALSE)
    "getMesh" = function(normals = TRUE, rgl = TRUE, ...) {
      stopifnot(isBoolean(normals))
      stopifnot(isBoolean(rgl))
      mesh <- private[[".CGALmesh"]]$getRmesh(normals)
      if(rgl) {
        if(is.matrix(mesh[["faces"]])) {
          nsides <- nrow(mesh[["faces"]])
          if(nsides == 3L) {
            mesh <- mesh3d(
              x         = t(mesh[["vertices"]]),
              triangles = mesh[["faces"]],
              normals   = mesh[["normals"]],
              ...
            )
          } else {
            mesh <- mesh3d(
              x       = t(mesh[["vertices"]]),
              quads   = mesh[["faces"]],
              normals = mesh[["normals"]],
              ...
            )
          }
        } else {
          faces <- split(mesh[["faces"]], lengths(mesh[["faces"]]))
          if(all(names(faces) %in% c("3", "4"))) {
            mesh <- mesh3d(
              x         = t(mesh[["vertices"]]),
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
      private[[".CGALmesh"]]$isClosed()
    },
    
    #' @description Check whether the mesh is triangle.
    #' @return A Boolean value, whether the mesh is triangle.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())
    #' mesh$isTriangle()
    "isTriangle" = function() {
      private[[".CGALmesh"]]$isTriangle()
    },
    
    #' @description Reverse the orientation of the faces of the mesh.
    #' @return The modified \code{cgalMesh} object. \strong{WARNING}: even if 
    #'   you store the result in a new variable, the original mesh is modified. 
    "reverseOrientation" = function() {
      private[[".CGALmesh"]]$reverseFaceOrientations()
      invisible(self)
    },
    
    #' @description Check whether the mesh self-intersects.
    #' @return A Boolean value, whether the mesh self-intersects.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())
    #' mesh$selfIntersects()
    "selfIntersects" = function() {
      private[[".CGALmesh"]]$doesSelfIntersect()
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
      private[[".CGALmesh"]]$triangulate()
      invisible(self)
    },

    #' @description Get the vertices of the mesh.
    #' @return The vertices in a matrix.
    "vertices" = function() {
      t(private[[".CGALmesh"]]$vertices())
    },
    
    #' @description Write mesh to a file.
    #' @param filename path to the file to be written, with extension 
    #'   \code{off} or \code{ply}
    #' @param precision a positive integer, the desired number of decimal 
    #'   places
    #' @param binary Boolean, whether to write a binary PLY file if 
    #'   \code{filename} has the \code{ply} extension
    #' @return Nothing, just writes a file.
    "writeMeshFile" = function(filename, precision = 17, binary = FALSE) {
      stopifnot(isString(filename))
      stopifnot(isPositiveInteger(precision))
      stopifnot(isBoolean(binary))
      filename <- path.expand(filename)
      private[[".CGALmesh"]]$writeFile(filename, as.integer(precision))
    }
    
  )
  
)