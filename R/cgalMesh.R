getXPtr <- function(cMesh){
  cMesh[[".__enclos_env__"]][["private"]][[".CGALmesh"]][["xptr"]]
}

#' @title R6 class to represent a CGAL mesh
#' @description R6 class to represent a CGAL mesh.
#'
#' @importFrom R6 R6Class
#' @importFrom rgl mesh3d
#' @importFrom grDevices col2rgb
#' @importFrom tools file_ext
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
    #'   integers or list of vectors of integers), and optionally a field 
    #'   \code{normals} (numeric matrix with three columns); if 
    #'   this argument is a \strong{rgl} mesh containing some colors, these 
    #'   colors will be assigned to the created \code{cgalMesh} object
    #' @param vertices if \code{mesh} is missing, must be a numeric matrix with 
    #'   three columns
    #' @param faces if \code{mesh} is missing, must be either a matrix of 
    #'   integers (each row gives the vertex indices of a face) or a list of 
    #'   vectors of integers (each one gives the vertex indices of a face)
    #' @param normals if \code{mesh} is missing, must be \code{NULL} or a 
    #'   numeric matrix with three columns and as many rows as vertices
    #' @param clean Boolean, no effect if the mesh is given by a file, 
    #'   otherwise it indicates whether to clean the mesh (merge duplicated 
    #'   vertices and duplicated faces, remove isolated vertices); set to 
    #'   \code{FALSE} if you know your mesh is already clean
    #' @return A \code{cgalMesh} object.
    #' @examples 
    #' library(cgalMeshes)
    #' meshFile <- system.file(
    #'   "extdata", "bigPolyhedron.off", package = "cgalMeshes"
    #' )
    #' mesh <- cgalMesh$new(meshFile)
    #' rglmesh <- mesh$getMesh()
    #' \donttest{library(rgl)
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglmesh, color = "tomato")
    #' plotEdges(
    #'   mesh$getVertices(), mesh$getEdges(), color = "darkred"
    #' )}
    #' 
    #' # this one has colors: ####
    #' meshFile <- system.file(
    #'   "extdata", "pentagrammicDipyramid.ply", package = "cgalMeshes"
    #' )
    #' mesh <- cgalMesh$new(meshFile)
    #' rmesh <- mesh$getMesh()
    #' \donttest{library(rgl)
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.85)
    #' shade3d(rmesh, meshColor = "faces")}
    "initialize" = function(
      mesh, vertices, faces, normals = NULL, clean = FALSE
    ){
      # one can also initialize from an external pointer, but 
      # this is hidden to the user
      if(inherits(clean, "externalptr")) {
        private[[".CGALmesh"]] <- CGALmesh$new(clean)
        return(invisible(self))
      }
      stopifnot(isBoolean(clean))
      vcolors <- NULL
      fcolors <- NULL
      if(!missing(mesh)) {
        if(inherits(mesh, "mesh3d")) {
          normals <- mesh[["normals"]][1L:3L, ]
          colors <- mesh[["material"]]$color
          if(!is.null(colors)){
            nv <- ncol(mesh[["vb"]])
            if(length(colors) == nv) {
              vcolors <- colors
            }
            nf <- 0L
            if(!is.null(it <- mesh[["it"]])) {
              nf <- nf + ncol(it)
            }
            if(!is.null(ib <- mesh[["ib"]])) {
              nf <- nf + ncol(ib)
            }
            if(length(colors) == nf) {
              fcolors <- colors
            }
          }
          mesh <- getVF(mesh)
        } else if(is.list(mesh)) {
          normals <- mesh[["normals"]]
          if(!is.null(normals)) {
            normals <- t(normals)
          }
        }
        if(is.list(mesh)) {
          VF <- checkMesh(mesh[["vertices"]], mesh[["faces"]], aslist = TRUE)
        } else if(isFilename(mesh)) {
          binary <- FALSE
          if(tolower(file_ext(mesh)) == "ply") {
            r <- readBin(mesh, "raw", 20L)
            binary <- grepl("binary", rawToChar(r))
          }
          private[[".CGALmesh"]] <- CGALmesh$new(path.expand(mesh), binary)
          return(invisible(self))
        } else {
          stop("Invalid `mesh` argument.")
        }
      } else {
        VF <- checkMesh(vertices, faces, aslist = TRUE)
        if(!is.null(normals)) {
          normals <- t(normals)
        }
      }
      if(!is.null(normals)) {
        if(nrow(normals) != 3L) {
          stop("The normals must be three-dimensional")
        }
        if(ncol(normals) != ncol(VF[["vertices"]])) {
          stop("The number of normals does not match the number of vertices.")
        }
      }
      private[[".CGALmesh"]] <- 
        CGALmesh$new(
          VF[["vertices"]], VF[["faces"]], clean, normals, vcolors, fcolors
        )
      invisible(self)
    },
    
    #' @description Print a \code{cgalMesh} object.
    #' @param ... ignored
    #' @return No value returned, just prints some information about the mesh.
    "print" = function(...) {
      private[[".CGALmesh"]]$print()
    },

    #' @description Compute the area of the mesh. The mesh must be triangle 
    #'   and must not self-intersect.
    #' @return A number, the mesh area.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())$triangulate()
    #' mesh$area()
    "area" = function() {
      private[[".CGALmesh"]]$area()
    },

    #' @description Assign colors (or any character strings) to the faces of 
    #'   the mesh.
    #' @param colors a character vector whose length equals the number 
    #'  of faces, or a single character string to be assigned to each face 
    #' @return The current \code{cgalMesh} object, invisibly.
    "assignFaceColors" = function(colors) {
      stopifnot(isStringVector(colors))
      . <- private[[".CGALmesh"]]$assignFaceColors(colors)
      invisible(self)
    },

    #' @description Assign scalars to the faces of the mesh.
    #' @param scalars a numeric vector whose length equals the number of faces
    #' @return The current \code{cgalMesh} object, invisibly.
    "assignFaceScalars" = function(scalars) {
      stopifnot(is.numeric(scalars))
      . <- private[[".CGALmesh"]]$assignFaceScalars(scalars)
      invisible(self)
    },

    #' @description Assign per-vertex normals to the mesh.
    #' @param normals a numeric matrix with three columns and as many rows as 
    #'   the number of vertices
    #' @return The current \code{cgalMesh} object, invisibly.
    "assignNormals" = function(normals) {
      stopifnot(is.matrix(normals), ncol(normals) == 3L)
      storage.mode(normals) <- "double"
      . <- private[[".CGALmesh"]]$assignNormals(t(normals))
      invisible(self)
    },
    
    #' @description Assign colors (or any character strings) to the vertices of 
    #'   the mesh.
    #' @param colors a character vector whose length equals the number of vertices
    #' @return The current \code{cgalMesh} object, invisibly.
    "assignVertexColors" = function(colors) {
      stopifnot(isStringVector(colors))
      . <- private[[".CGALmesh"]]$assignVertexColors(colors)
      invisible(self)
    },

    #' @description Assign scalars to the vertices of the mesh.
    #' @param scalars a numeric vector whose length equals the number 
    #' of vertices
    #' @return The current \code{cgalMesh} object, invisibly.
    "assignVertexScalars" = function(scalars) {
      stopifnot(is.numeric(scalars))
      . <- private[[".CGALmesh"]]$assignVertexScalars(scalars)
      invisible(self)
    },

    #' @description Bounding box of the mesh.
    #' @return A list containing the smallest corner point and the largest 
    #'   corner point of the bounding box, named \code{lcorner} and 
    #'  \code{ucorner} respectively. Use \code{\link{isoCuboidMesh}} to get a 
    #'  mesh of this bounding box.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' rmesh <- cyclideMesh(a = 97, c = 32, mu = 57)
    #' mesh <- cgalMesh$new(rmesh)
    #' bbox <- mesh$boundingBox()
    #' bxmesh <- isoCuboidMesh(bbox[["lcorner"]], bbox[["ucorner"]])
    #' open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(0, -60)
    #' shade3d(rmesh, color = "gold")
    #' wire3d(bxmesh, color = "black")}
    "boundingBox" = function() {
      private[[".CGALmesh"]]$boundingBox()
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
      private[[".CGALmesh"]]$doesBoundVolume()
    },

    #' @description Performs the Catmull-Clark subdivision and deformation. 
    #'   The mesh must be triangle.
    #' @param iterations number of iterations
    #' @return The modified reference mesh, invisibly.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' hopfMesh <- HopfTorusMesh(nu = 80, nv = 40)
    #' mesh <- cgalMesh$new(hopfMesh)
    #' mesh$CatmullClark(iterations = 2)
    #' mesh$computeNormals()
    #' rmesh <- mesh$getMesh()
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 800, 400))
    #' mfrow3d(1, 2)
    #' view3d(0, 0, zoom = 0.9)
    #' shade3d(hopfMesh, color = "red")
    #' wire3d(hopfMesh, color = "black")
    #' next3d()
    #' view3d(0, 0, zoom = 0.9)
    #' shade3d(rmesh, color = "red")
    #' wire3d(rmesh, color = "black")}
    "CatmullClark" = function(iterations = 1) {
      stopifnot(isStrictPositiveInteger(iterations))
      private[[".CGALmesh"]]$CatmullClark(as.integer(iterations))
      invisible(self)
    },
    
    #' @description Centroid of the mesh. The mesh must be triangle.
    #' @return The Cartesian coordinates of the centroid of the mesh.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' mesh <- cgalMesh$new(icosahedron3d())
    #' mesh$centroid()}
    "centroid" = function() {
      private[[".CGALmesh"]]$centroid()
    },
    
    #' @description Clip mesh to the volume bounded by another mesh. The
    #'   mesh must be triangle. Face properties (colors and scalars) 
    #'   are preserved. 
    #'   \strong{WARNING}: the reference mesh is then replaced by its 
    #'   clipped version.
    #'
    #' @param clipper a \code{cgalMesh} object; it must represent a closed 
    #'   triangle mesh which doesn't self-intersect
    #' @param clipVolume Boolean, whether the clipping has to be done on the 
    #'   volume bounded by the reference mesh rather than on its surface (i.e. 
    #'   the reference mesh will be kept closed if it is closed); if 
    #'   \code{TRUE}, the mesh to be clipped must not self-intersect
    #' @return The reference mesh is always replaced by the result of the 
    #'   clipping. If \code{clipVolume=TRUE}, this function returns two 
    #'   \code{cgalMesh} objects: the two parts of the clipped mesh contained 
    #'   in the reference mesh and the clipping mesh respectively. Otherwise, 
    #'   this function returns the modified reference mesh.
    #' @examples 
    #' # cube clipped to sphere ####
    #' library(cgalMeshes)
    #' library(rgl)
    #' mesh    <- cgalMesh$new(cube3d())$triangulate()
    #' clipper <- cgalMesh$new(sphereMesh(r= sqrt(2)))
    #' mesh$assignFaceColors("blue")
    #' clipper$assignFaceColors("red")
    #' meshes <- mesh$clip(clipper, clipVolume = TRUE)
    #' mesh1 <- meshes[[1]]
    #' mesh2 <- meshes[[2]]
    #' mesh2$computeNormals()
    #' rglmesh1 <- mesh1$getMesh()
    #' rglmesh2 <- mesh2$getMesh()
    #' \donttest{open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(45, 45, zoom = 0.9)
    #' shade3d(rglmesh1, meshColor = "faces")
    #' shade3d(rglmesh2, meshColor = "faces")}
    #' 
    #' # Togliatti surface clipped to a ball ####
    #' \donttest{library(rmarchingcubes)
    #' library(rgl)
    #' library(cgalMeshes)
    #' # Togliatti surface equation: f(x,y,z) = 0
    #' f <- function(x, y, z) {
    #'   64*(x-1) *
    #'     (x^4 - 4*x^3 - 10*x^2*y^2 - 4*x^2 + 16*x - 20*x*y^2 + 5*y^4 + 16 - 20*y^2) - 
    #'     5*sqrt(5-sqrt(5))*(2*z - sqrt(5-sqrt(5))) * 
    #'     (4*(x^2 + y^2 - z^2) + (1 + 3*sqrt(5)))^2
    #' }
    #' # grid
    #' n <- 200L
    #' x <- y <- seq(-5, 5, length.out = n)
    #' z <- seq(-4, 4, length.out = n)
    #' Grid <- expand.grid(X = x, Y = y, Z = z)
    #' # calculate voxel
    #' voxel <- array(with(Grid, f(X, Y, Z)), dim = c(n, n, n))
    #' # calculate isosurface
    #' contour_shape <- contour3d(
    #'   griddata = voxel, level = 0, x = x, y = y, z = z
    #' )
    #' # make rgl mesh (plotted later)
    #' rglMesh <- tmesh3d(
    #'   vertices = t(contour_shape[["vertices"]]),
    #'   indices  = t(contour_shape[["triangles"]]),
    #'   normals  = contour_shape[["normals"]],
    #'   homogeneous = FALSE
    #' )
    #' # make CGAL mesh
    #' mesh <- cgalMesh$new(rglMesh)
    #' # clip to sphere of radius 4.8
    #' sphere <- sphereMesh(r = 4.8)
    #' clipper <- cgalMesh$new(sphere)
    #' mesh$clip(clipper, clipVolume = FALSE)
    #' rglClippedMesh <- mesh$getMesh()
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 900, 450))
    #' mfrow3d(1L, 2L)
    #' view3d(0, -70, zoom = 0.8)
    #' shade3d(rglMesh, color = "firebrick")
    #' next3d()
    #' view3d(0, -70, zoom = 0.8)
    #' shade3d(rglClippedMesh, color = "firebrick")
    #' shade3d(sphere, color = "yellow", alpha = 0.15)}
    "clip" = function(clipper, clipVolume) {
      stopifnot(isCGALmesh(clipper))
      stopifnot(isBoolean(clipVolume))
      clipperXPtr <- getXPtr(clipper)
      if(clipVolume) {
        xptrs <- private[[".CGALmesh"]]$clipMesh(clipperXPtr, TRUE)
        lapply(xptrs, function(xptr) cgalMesh$new(clean = xptr))
      } else {
        . <- private[[".CGALmesh"]]$clipMesh(clipperXPtr, FALSE)
        invisible(self)
      }
    },
    
    #' @description Clip the mesh to a plane. The mesh must be triangle. Face 
    #'   properties (colors, scalars) are preserved.
    #' @param planePoint numeric vector of length three, a point belonging 
    #'  to the plane
    #' @param planeNormal numeric vector of length three, a vector orthogonal
    #'   to the plane 
    #' @param clipVolume Boolean, whether to clip on the volume
    #' @return If \code{clipVolume=FALSE}, the modified reference mesh is 
    #'   invisibly returned. If \code{clipVolume=TRUE}, a list of two 
    #'   \code{cgalMesh} objects is returned: the first one is the part of 
    #'   the clipped mesh corresponding to the original mesh, the second 
    #'   one is the part of the clipped mesh corresponding to the plane.
    #' @examples 
    #' library(cgalMeshes)
    #' library(rgl)
    #' rmesh <- sphereMesh()
    #' mesh <- cgalMesh$new(rmesh)
    #' nfaces <- nrow(mesh$getFaces())
    #' if(require("randomcoloR")) {
    #'   colors <- 
    #'     randomColor(nfaces, hue = "random", luminosity = "dark")
    #' } else {
    #'   colors <- rainbow(nfaces)
    #' }
    #' mesh$assignFaceColors(colors)
    #' meshes <- mesh$clipToPlane(
    #'   planePoint  = c(0, 0, 0), 
    #'   planeNormal = c(0, 0, 1), 
    #'   clipVolume = TRUE
    #' )
    #' mesh1 <- meshes[[1]]
    #' mesh2 <- meshes[[2]]
    #' mesh1$computeNormals()
    #' rClippedMesh1 <- mesh1$getMesh()
    #' rClippedMesh2 <- mesh2$getMesh()
    #' \donttest{open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(70, 0)
    #' shade3d(rClippedMesh1, meshColor = "faces")
    #' shade3d(rClippedMesh2, color = "orange")}
    "clipToPlane" = function(planePoint, planeNormal, clipVolume) {
      check <- is.numeric(planePoint) && length(planePoint) == 3L && 
        !anyNA(planePoint)
      if(!check) {
        stop("Invalid `planePoint` vector.")
      }
      check <- is.numeric(planeNormal) && length(planeNormal) == 3L && 
        !anyNA(planeNormal)
      if(!check) {
        stop("Invalid `planeNormal` vector.")
      }
      if(c(crossprod(planeNormal)) == 0) {
        stop("The `planeNormal` vector cannot be null.")
      }
      stopifnot(isBoolean(clipVolume))
      if(clipVolume) {
        xptrs <- private[[".CGALmesh"]]$clipToPlane(
          planePoint, planeNormal, TRUE
        )
        lapply(xptrs, function(xptr) cgalMesh$new(clean = xptr))
      } else {
        . <- private[[".CGALmesh"]]$clipToPlane(
          planePoint, planeNormal, FALSE
        )
        invisible(self)
      }
    },

    #' @description Clip the mesh to an iso-oriented cuboid. The mesh must be 
    #'   triangle. Face properties (colors, scalars) are preserved.
    #' @param lcorner,ucorner two diagonally opposite vertices of the 
    #'   iso-oriented cuboid, the smallest and the largest (see 
    #'   \code{\link{isoCuboidMesh}})
    #' @param clipVolume Boolean, whether to clip on the volume
    #' @return If \code{clipVolume=FALSE}, the modified reference mesh is 
    #'   invisibly returned. If \code{clipVolume=TRUE}, a list of two 
    #'   \code{cgalMesh} objects is returned: the first one is the part of 
    #'   the clipped mesh corresponding to the original mesh, the second 
    #'   one is the part of the clipped mesh corresponding to the cuboid.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' rmesh <- HopfTorusMesh(nu = 200, nv = 200)
    #' mesh <- cgalMesh$new(rmesh)
    #' mesh$assignFaceColors("orangered")
    #' lcorner <- c(-7, -7, -5)
    #' ucorner <- c(7, 6, 5)
    #' bxmesh <- isoCuboidMesh(lcorner, ucorner)
    #' mesh$clipToIsoCuboid(
    #'   lcorner, ucorner, clipVolume = FALSE
    #' )
    #' mesh$computeNormals()
    #' rClippedMesh <- mesh$getMesh()
    #' open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(-40, 0)
    #' shade3d(rClippedMesh, meshColor = "faces")
    #' shade3d(bxmesh, color = "cyan", alpha = 0.3)}
    "clipToIsoCuboid" = function(lcorner, ucorner, clipVolume) {
      check <- is.numeric(lcorner) && length(lcorner) == 3L && 
        !anyNA(lcorner)
      if(!check) {
        stop("Invalid `lcorner` vector.")
      }
      check <- is.numeric(ucorner) && length(ucorner) == 3L && 
        !anyNA(ucorner)
      if(!check) {
        stop("Invalid `ucorner` vector.")
      }
      stopifnot(all(ucorner > lcorner))
      stopifnot(isBoolean(clipVolume))
      if(clipVolume) {
        xptrs <- private[[".CGALmesh"]]$clipToIsoCuboid(
          lcorner, ucorner, TRUE
        )
        lapply(xptrs, function(xptr) cgalMesh$new(clean = xptr))
      } else {
        . <- private[[".CGALmesh"]]$clipToIsoCuboid(
          lcorner, ucorner, FALSE
        )
        invisible(self)
      }
    },

    #' @description tmp. 
    #' @return tmp.
    "collectGarbage" = function() {
      . <- private[[".CGALmesh"]]$collectGarbage()
      invisible(self)
    },
    
    #' @description Compute per-vertex normals of the mesh. 
    #' @return The current \code{cgalMesh} object, invisibly. 
    #'  To get the normals, use the \code{getNormals} method.
    "computeNormals" = function() {
      . <- private[[".CGALmesh"]]$computeNormals()
      invisible(self)
    },
    
    #' @description Decomposition into connected components. All face 
    #'   properties (colors, scalars) and vertex properties 
    #'   (colors, scalars, normals) are preserved.
    #' @param triangulate Boolean, whether to triangulate the connected 
    #'   components.
    #' @return A list of \code{cgalMesh} objects, one for each connected 
    #'   component.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rmarchingcubes)
    #' # isosurface function (slice of a seven-dimensional toratope)
    #' f <- function(x, y, z, a) {
    #'   (sqrt(
    #'     (sqrt((sqrt((x*sin(a))^2 + (z*cos(a))^2) - 5)^2 + (y*sin(a))^2) - 2.5)^2 + 
    #'       (x*cos(a))^2) - 1.25
    #'   )^2 + (sqrt((sqrt((z*sin(a))^2 + (y*cos(a))^2) - 2.5)^2) - 1.25)^2
    #' }
    #' # make grid
    #' n <- 200L
    #' x <- seq(-10, 10, len = n)
    #' y <- seq(-10, 10, len = n)
    #' z <- seq(-10, 10, len = n)
    #' Grid <- expand.grid(X = x, Y = y, Z = z)
    #' # compute isosurface
    #' voxel <- array(with(Grid, f(X, Y, Z, a = pi/2)), dim = c(n, n, n))
    #' isosurface <- contour3d(voxel, level = 0.25, x = x, y = y, z = z)
    #' # make CGAL mesh
    #' mesh <- cgalMesh$new(
    #'   vertices = isosurface[["vertices"]], 
    #'   faces = isosurface[["triangles"]],
    #'   normals = isosurface[["normals"]]
    #' )
    #' # connected components
    #' components <- mesh$connectedComponents()
    #' ncc <- length(components)
    #' # plot
    #' library(rgl)
    #' colors <- rainbow(ncc)
    #' open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(30, 50)
    #' for(i in 1L:ncc) {
    #'   rglMesh <- components[[i]]$getMesh()
    #'   shade3d(rglMesh, color = colors[i])
    #' }}
    "connectedComponents" = function(triangulate = FALSE) {
      stopifnot(isBoolean(triangulate))
      xptrs <- private[[".CGALmesh"]]$connectedComponents(triangulate)
      lapply(xptrs, function(xptr) cgalMesh$new(clean = xptr))
    },
    
    #' @description Decomposition into convex parts. The mesh must be triangle.
    #' @param triangulate Boolean, whether to triangulate the convex parts
    #' @return A list of \code{cgalMesh} objects, one for each convex part.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' mesh <- cgalMesh$new(pentagrammicPrism)$triangulate()
    #' cxparts <- mesh$convexParts()
    #' ncxparts <- length(cxparts)
    #' colors <- hcl.colors(ncxparts, palette = "plasma")
    #' open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(20, -20, zoom = 0.8)
    #' for(i in 1L:ncxparts) {
    #'   cxmesh <- cxparts[[i]]$getMesh()
    #'   shade3d(cxmesh, color = colors[i])
    #' }}
    "convexParts" = function(triangulate = TRUE) {
      stopifnot(isBoolean(triangulate))
      xptrs <- private[[".CGALmesh"]]$convexParts(triangulate)
      lapply(xptrs, function(xptr) cgalMesh$new(clean = xptr))
    },

    #' @description Copy the mesh. This can change the order of the vertices.
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
    
    #' @description Distance from one or more points to the mesh. The mesh 
    #'   must be triangle.
    #' @param points either one point given as a numeric vector or several 
    #'   points given as a numeric matrix with three columns
    #' @return A numeric vector providing the distances between the given 
    #'   point(s) to the mesh.
    #' @examples
    #' # cube example ####
    #' library(cgalMeshes)
    #' mesh <- cgalMesh$new(rgl::cube3d())$triangulate()
    #' points <- rbind(
    #'   c(0, 0, 0),
    #'   c(1, 1, 1)
    #' )
    #' mesh$distance(points) # should be 1 and 0
    #'
    #' # cyclide example ####
    #' \donttest{library(cgalMeshes)
    #' a <- 100; c <- 30; mu <- 80
    #' mesh <- cgalMesh$new(cyclideMesh(a, c, mu, nu = 100L, nv = 100L))
    #' O2 <- c(c, 0, 0)
    #' # should be a - mu = 20 (see ?cyclideMesh):
    #' mesh$distance(O2)}    
    "distance" = function(points){
      if(!is.matrix(points)){
        points <- rbind(points)
      }
      if(ncol(points) != 3L){
        stop(
          "The `points` argument must be a vector of length three or ",
          "a matrix with three columns."
        )
      }
      stopifnot(is.numeric(points))
      storage.mode(points) <- "double"
      if(anyNA(points)){
        stop("Found missing values in `points`.")
      }
      private[[".CGALmesh"]]$distance(t(points))
    },

    #' @description Performs the Doo-Sabin subdivision and deformation.
    #' @param iterations number of iterations
    #' @return The modified reference mesh, invisibly.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' hopfMesh <- HopfTorusMesh(nu = 80, nv = 40)
    #' mesh <- cgalMesh$new(hopfMesh)
    #' mesh$DooSabin(iterations = 2)
    #' mesh$triangulate()
    #' mesh$computeNormals()
    #' rmesh <- mesh$getMesh()
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 800, 400))
    #' mfrow3d(1, 2)
    #' view3d(0, 0, zoom = 0.9)
    #' shade3d(hopfMesh, color = "brown")
    #' wire3d(hopfMesh, color = "black")
    #' next3d()
    #' view3d(0, 0, zoom = 0.9)
    #' shade3d(rmesh, color = "brown")
    #' wire3d(rmesh, color = "black")}
    "DooSabin" = function(iterations = 1) {
      stopifnot(isStrictPositiveInteger(iterations))
      private[[".CGALmesh"]]$DooSabin(as.integer(iterations))
      invisible(self)
    },
    
    #' @description Faces containing a given vertex.
    #' @param v a vertex index
    #' @return An integer vector, the indices of the faces containing \code{v}.
    "facesAroundVertex" = function(v) {
      stopifnot(isStrictPositiveInteger(v))
      private[[".CGALmesh"]]$facesAroundVertex(as.integer(v) - 1L)
    },
    
    #' @description Fair a region of the mesh, i.e. make it smooth. The mesh 
    #'   must be triangle. This modifies the reference mesh. All 
    #'   face properties and vertex properties except the normals 
    #'   are preserved.
    #' @param indices the indices of the vertices in the region to be faired
    #' @return The modified \code{cgalMesh} object.
    #' @examples 
    #' library(cgalMeshes)
    #' rglHopf <- HopfTorusMesh(nu = 100, nv = 100)
    #' hopf <- cgalMesh$new(rglHopf)
    #' # squared norms of the vertices
    #' normsq <- apply(hopf$getVertices(), 1L, crossprod)
    #' # fair the region where the squared norm is > 19
    #' indices <- which(normsq > 19)
    #' hopf$fair(indices)
    #' rglHopf_faired <- hopf$getMesh()
    #' # plot
    #' \donttest{library(rgl)
    #' open3d(windowRect = 50 + c(0, 0, 900, 450))
    #' mfrow3d(1L, 2L)
    #' view3d(0, 0, zoom = 0.8)
    #' shade3d(rglHopf, color = "orangered")
    #' next3d()
    #' view3d(0, 0, zoom = 0.8)
    #' shade3d(rglHopf_faired, color = "orangered")}
    "fair" = function(indices) {
      stopifnot(isAtomicVector(indices))
      stopifnot(is.numeric(indices))
      integers <- isTRUE(all.equal(indices, floor(indices)))
      if(!integers) {
        stop("The indices must be positive integers.")
      }
      positive <- isTRUE(all(indices >= 1))
      if(!positive) {
        stop("The indices must be positive integers.")
      }
      private[[".CGALmesh"]]$fair(as.integer(indices) - 1L)
      invisible(self)
    },

    #' @description Fill a hole in the mesh. The face properties and the 
    #'   vertex properties are preserved.
    #' @param border index of the boundary cycle forming the hole to be 
    #'   filled; the boundary cycles can be identified with 
    #'   \code{$getBorders()} 
    #' @param fair Boolean, whether to fair (i.e. smooth) the filled hole
    #' @return The filled hole as a \code{cgalMesh} object. The reference 
    #'   mesh is updated.
    #' @examples 
    #' library(cgalMeshes)
    #' library(rgl)
    #' # make a sphere
    #' sphere <- sphereMesh()
    #' mesh <- cgalMesh$new(sphere)
    #' # make a hole in this sphere
    #' mesh$clipToPlane(
    #'   planePoint  = c(0.5, 0, 0),
    #'   planeNormal = c(1, 0, 0),
    #'   clipVolume  = FALSE
    #' )
    #' mesh$computeNormals()
    #' rmesh <- mesh$getMesh()
    #' # fill the hole
    #' hole <- mesh$fillBoundaryHole(1, fair = TRUE)
    #' hole$computeNormals()
    #' rhole <- hole$getMesh()
    #' # plot
    #' \donttest{open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(30, 30)
    #' shade3d(rmesh, color = "red")
    #' shade3d(rhole, color = "blue")}
    "fillBoundaryHole" = function(border, fair = TRUE) {
      stopifnot(isStrictPositiveInteger(border))
      stopifnot(isBoolean(fair))
      xptr <- private[[".CGALmesh"]]$fillBoundaryHole(
        as.integer(border) - 1L, fair
      )
      cgalMesh$new(clean = xptr)
    },
    
    #' @description Replace missing face colors with a color.
    #' @param color the color to replace the missing face colors
    #' @return The reference mesh, invisibly.
    "fillFaceColors" = function(color) {
      stopifnot(isString(color))
      colors <- self$getFaceColors()
      if(is.null(colors)) {
        message("The mesh has no face colors.")
      } else {
        colors[colors == ""] <- color
        self$assignFaceColors(colors)
      }
      invisible(self)
    },
    
    #' @description Split the mesh into two meshes according to a 
    #'   given set of selected faces. Face properties are 
    #'   preserved.
    #' @param faces a vector of face indices
    #' @return Two \code{cgalMesh} objects. The first one is the mesh consisting 
    #'  of the faces of the reference mesh given in the \code{faces} 
    #'  argument. The second one is the complementary mesh.
    #' @examples 
    #' \donttest{library(rgl)
    #' library(cgalMeshes)
    #' rmesh <- HopfTorusMesh(nu = 80, nv = 60)
    #' mesh <- cgalMesh$new(rmesh)
    #' areas <- mesh$getFacesInfo()[, "area"]
    #' bigFaces <- which(areas > 1)
    #' meshes <- mesh$filterMesh(bigFaces)
    #' rmesh1 <- meshes[[1]]$getMesh()
    #' rmesh2 <- meshes[[2]]$getMesh()
    #' open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(0, 0)
    #' shade3d(rmesh1, color = "red")
    #' shade3d(rmesh2, color = "blue")
    #' wire3d(rmesh)}
    "filterMesh" = function(faces) {
      stopifnot(isAtomicVector(faces))
      stopifnot(is.numeric(faces))
      integers <- isTRUE(all.equal(faces, floor(faces)))
      if(!integers) {
        stop("The face indices must be positive integers.")
      }
      positive <- isTRUE(all(faces >= 1))
      if(!positive) {
        stop("The face indices must be positive integers.")
      }
      xptrs <- private[[".CGALmesh"]]$filterMesh(as.integer(faces) - 1L)
      list(
        "mesh1" = cgalMesh$new(clean = xptrs[["fmesh1"]]),
        "mesh2" = cgalMesh$new(clean = xptrs[["fmesh2"]])
      )
    },

    #' @description Duplicate non-manifold vertices.
    #' @return The possibly modified reference mesh, invisibly.
    "fixManifoldness" = function() {
      private[[".CGALmesh"]]$fixManifoldness()
      invisible(self)
    },
    
    #' @description Estimated geodesic distances between vertices. The mesh 
    #'   must be triangle.
    #' @param index index of the source vertex
    #' @return The estimated geodesic distances from the source vertex to each
    #'   vertex.
    #' @examples 
    #' \donttest{# torus ####
    #' library(cgalMeshes)
    #' library(rgl)
    #' rglmesh <- torusMesh(R = 3, r = 2, nu = 90, nv = 60)
    #' mesh <- cgalMesh$new(rglmesh)
    #' # estimated geodesic distances
    #' geodists <- mesh$geoDists(1L)
    #' # normalization to (0, 1)
    #' geodists <- geodists / max(geodists)
    #' # color each vertex according to its geodesic distance from the source
    #' fcolor <- colorRamp(viridisLite::turbo(200L))
    #' colors <- fcolor(geodists)
    #' colors <- rgb(colors[, 1L], colors[, 2L], colors[, 3L], maxColorValue = 255)
    #' rglmesh[["material"]] <- list("color" = colors)
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
    #' shade3d(rglmesh)
    #' wire3d(rglmesh, color = "black")
    #' if(!rgl.useNULL()) {
    #'   play3d(spin3d(axis = c(1, 1, 1), rpm = 5), duration = 20)  
    #' }}
    #' 
    #' # a trefoil knot (taken from `?rgl::cylinder3d`) ####
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' theta <- seq(0, 2*pi, length.out = 50L)
    #' knot <- cylinder3d(
    #'   center = cbind(
    #'     sin(theta) + 2*sin(2*theta), 
    #'     2*sin(3*theta), 
    #'     cos(theta) - 2*cos(2*theta)),
    #'   e1 = cbind(
    #'     cos(theta) + 4*cos(2*theta), 
    #'     6*cos(3*theta), 
    #'     sin(theta) + 4*sin(2*theta)),
    #'   radius = 0.8, 
    #'   closed = TRUE)
    #' knot <- subdivision3d(knot, depth = 2)
    #' mesh <- cgalMesh$new(knot)$triangulate()
    #' rglmesh <- mesh$getMesh()
    #' # estimated geodesic distances
    #' geodists <- mesh$geoDists(1L)
    #' # normalization to (0, 1)
    #' geodists <- geodists / max(geodists)
    #' # color each vertex according to its geodesic distance from the source
    #' fcolor <- colorRamp(viridisLite::inferno(200L))
    #' colors <- fcolor(geodists)
    #' colors <- rgb(colors[, 1L], colors[, 2L], colors[, 3L], maxColorValue = 255)
    #' rglmesh[["material"]] <- list("color" = colors)
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
    #' shade3d(rglmesh)
    #' if(!rgl.useNULL()) {
    #'   play3d(spin3d(axis = c(1, 1, 0), rpm = 5), duration = 20)  
    #' }}
    "geoDists" = function(index) {
      stopifnot(isStrictPositiveInteger(index))
      private[[".CGALmesh"]]$geoDists(as.integer(index) - 1L)
    },

    #' @description Get the borders of the mesh.
    #' @return A list of matrices representing the boundary cycles. Each matrix 
    #'  has three columns: \code{"edge"}, an edge index, and 
    #'  \code{"v1"} and \code{"v2"}, the vertex indices of this edge.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' # isosurface f=0
    #' f <- function(x, y, z) {
    #'   sin_x <- sin(x)
    #'   sin_y <- sin(y)
    #'   sin_z <- sin(z)
    #'   cos_x <- cos(x)
    #'   cos_y <- cos(y)
    #'   cos_z <- cos(z)
    #'   d <- sqrt(
    #'     (-sin_x * sin_y + cos_x * cos_z) ** 2
    #'     + (-sin_y * sin_z + cos_y * cos_x) ** 2
    #'     + (-sin_z * sin_x + cos_z * cos_y) ** 2
    #'   )
    #'   (
    #'     cos(
    #'       x - (-sin_x * sin_y + cos_x * cos_z) / d
    #'     )
    #'     * sin(
    #'       y - (-sin_y * sin_z + cos_y * cos_x) / d
    #'     )
    #'     + cos(
    #'       y - (-sin_y * sin_z + cos_y * cos_x) / d
    #'     )
    #'     * sin(
    #'       z - (-sin_z * sin_x + cos_z * cos_y)/ d
    #'     )
    #'     + cos(
    #'       z - (-sin_z * sin_x + cos_z * cos_y) / d
    #'     )
    #'     * sin(
    #'       x - (-sin_x * sin_y + cos_x * cos_z) / d
    #'     )
    #'   ) * (
    #'     (
    #'       cos(
    #'         x + (-sin_x * sin_y + cos_x * cos_z) / d
    #'       )
    #'       * sin(
    #'         y + (-sin_y * sin_z + cos_y * cos_x) / d
    #'       )
    #'       + cos(
    #'         y + (-sin_y * sin_z + cos_y * cos_x) / d
    #'       )
    #'       * sin(
    #'         z + (-sin_z * sin_x + cos_z * cos_y) / d
    #'       )
    #'       + cos(
    #'         z + (-sin_z * sin_x + cos_z * cos_y) / d
    #'       )
    #'       * sin(
    #'         x + (-sin_x * sin_y + cos_x * cos_z) / d
    #'       )
    #'     )
    #'   )
    #' }
    #' # construct the isosurface f=0
    #' ngrid <- 200L
    #' x <- y <- z <- seq(-8.1, 8.1, length.out = ngrid)
    #' Grid <- expand.grid(X = x, Y = y, Z = z)
    #' voxel <- array(
    #'   with(Grid, f(X, Y, Z)), dim = c(ngrid, ngrid, ngrid)
    #' )
    #' library(rmarchingcubes)
    #' contour_shape <- contour3d(
    #'   griddata = voxel, level = 0,
    #'   x = x, y = y, z = z
    #' )
    #' # make mesh
    #' mesh <- cgalMesh$new(
    #'   list(
    #'     "vertices" = contour_shape[["vertices"]],
    #'     "faces"    = contour_shape[["triangles"]]
    #'   )
    #' )
    #' # clip the mesh to the ball of radius 8
    #' spheremesh <- cgalMesh$new(sphereMesh(r = 8))
    #' mesh$clip(spheremesh, clipVolume = FALSE)
    #' # compute normals
    #' mesh$computeNormals()
    #' # we will plot the borders
    #' borders <- mesh$getBorders()
    #' # plot
    #' rmesh <- mesh$getMesh()
    #' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.7)
    #' shade3d(rmesh, color = "darkred")
    #' vertices <- mesh$getVertices()
    #' for(border in borders){
    #'   plotEdges(
    #'     vertices, border[, c("v1", "v2")], color = "gold",
    #'     lwd = 3, edgesAsTubes = FALSE, verticesAsSpheres = FALSE
    #'   )
    #' }}
    "getBorders" = function() {
      private[[".CGALmesh"]]$getBorders()
    },
    
    #' @description Get the edges of the mesh.
    #' @return A dataframe with five columns; the first two ones give the 
    #'   vertex indices of each edge (one edge per row), the third one gives 
    #'   the lengths of each edge, the fourth one indicates whether the edges 
    #'   is a border edge, and the fifth one gives the dihedral angles 
    #'   in degrees between the two faces adjacent to each edge 
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(dodecahedron3d())
    #' head(mesh$getEdges())
    "getEdges" = function() {
      private[[".CGALmesh"]]$edges()
    },
    
    #' @description Get the faces of the mesh.
    #' @return The faces in a matrix if the mesh is triangle or quad, 
    #'   otherwise in a list.
    "getFaces" = function() {
      if(
        private[[".CGALmesh"]]$isTriangle() ||
        private[[".CGALmesh"]]$isQuad()
      ){
        private[[".CGALmesh"]]$getFacesMatrix()  
      } else {
        private[[".CGALmesh"]]$getFacesList()
      }
    },

    #' @description Get the centroids and the areas of the faces, for a 
    #'   triangle mesh only.
    #' @return A matrix with four columns: the first three ones provide the 
    #'   Cartesian coordinates of the centroids, the fourth one provides the 
    #'   areas.
    "getFacesInfo" = function() {
      private[[".CGALmesh"]]$getFacesInfo()
    },
    
    #' @description Get the face colors (if there are).
    #' @return The vector of colors (or any character vector) attached to 
    #'   the faces of the mesh, or \code{NULL} if nothing is assigned to 
    #'   the faces.
    "getFaceColors" = function() {
      private[[".CGALmesh"]]$getFcolors()
    },

    #' @description Get the face scalars (if there are).
    #' @return The vector of scalars attached to 
    #'   the faces of the mesh, or \code{NULL} if nothing is assigned to 
    #'   the faces.
    "getFaceScalars" = function() {
      private[[".CGALmesh"]]$getFscalars()
    },

    #' @description Get the vertex colors (if there are).
    #' @return The vector of colors (or any character vector) attached to 
    #'   the vertices of the mesh, or \code{NULL} if nothing is assigned to 
    #'   the vertices.
    "getVertexColors" = function() {
      private[[".CGALmesh"]]$getVcolors()
    },

    #' @description Get the vertex scalars (if there are).
    #' @return The vector of scalars attached to 
    #'   the vertices of the mesh, or \code{NULL} if nothing is assigned to 
    #'   the vertices.
    "getVertexScalars" = function() {
      private[[".CGALmesh"]]$getVscalars()
    },
    
    #' @description Get the per-vertex normals (if there are).
    #' @return The matrix of per-vertex normals if they have been given or 
    #'   computed (see \code{computeNormals}, or \code{NULL} otherwise.
    "getNormals" = function() {
      private[[".CGALmesh"]]$getNormals()
    },
    
    #' @description Get the mesh.
    #' @param rgl Boolean, whether to return a \strong{rgl} mesh if possible, 
    #'   i.e. if the mesh only has triangular or quadrilateral faces
    #' @param ... arguments passed to \code{\link[rgl:mesh3d]{mesh3d}} (if 
    #'   a \strong{rgl} mesh is returned)
    #' @return A \strong{rgl} mesh or a list with two or three fields: 
    #'   \code{vertices}, \code{faces}, and \code{normals} if XXXXXXXXXXXXXXXXXXXXXXXXX the argument 
    #'   \code{normals} is set to \code{TRUE}
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())$triangulate()
    #' mesh$getMesh()
    "getMesh" = function(rgl = TRUE, ...) {
      stopifnot(isBoolean(rgl))
      mesh <- private[[".CGALmesh"]]$getRmesh()
      if(!is.null(mesh[["normals"]])) {
        mesh[["normals"]] <- t(mesh[["normals"]])
      }
      if(rgl) {
        if(is.matrix(mesh[["faces"]])) {
          nsides <- nrow(mesh[["faces"]])
          if(nsides == 3L) {
            mesh <- mesh3d(
              x         = t(mesh[["vertices"]]),
              triangles = mesh[["faces"]],
              normals   = mesh[["normals"]],
              material  = list("color" = mesh[["colors"]]),
              ...
            )
          } else {
            mesh <- mesh3d(
              x        = t(mesh[["vertices"]]),
              quads    = mesh[["faces"]],
              normals  = mesh[["normals"]],
              material = list("color" = mesh[["colors"]]),
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
              material  = list("color" = mesh[["colors"]]),
              ...
            )
          } else {
            warning("Cannot make a rgl mesh.")
            mesh[["vertices"]] <- t(mesh[["vertices"]])
          }
        }
      } else {
        mesh[["vertices"]] <- t(mesh[["vertices"]])
      }
      mesh
    },
    
    #' @description Get the vertices of the mesh.
    #' @return The vertices in a matrix.
    "getVertices" = function() {
      t(private[[".CGALmesh"]]$getVertices())
    },
    
    #' @description Intersection with another mesh.
    #' @param mesh2 a \code{cgalMesh} object
    #' @return A \code{cgalMesh} object.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' # take two cubes
    #' rglmesh1 <- cube3d()
    #' rglmesh2 <- translate3d(cube3d(), 1, 1, 1)
    #' mesh1 <- cgalMesh$new(rglmesh1)
    #' mesh2 <- cgalMesh$new(rglmesh2)
    #' # the two meshes must be triangle
    #' mesh1$triangulate()
    #' mesh2$triangulate()
    #' # intersection
    #' imesh <- mesh1$intersection(mesh2)
    #' rglimesh <- imesh$getMesh()
    #' # extract edges for plotting
    #' extEdges <- exteriorEdges(imesh$getEdges())
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglimesh, color = "red")
    #' plotEdges(imesh$getVertices(), extEdges)
    #' shade3d(rglmesh1, color = "yellow", alpha = 0.2)
    #' shade3d(rglmesh2, color = "cyan", alpha = 0.2)}
    "intersection" = function(mesh2) {
      stopifnot(isCGALmesh(mesh2))
      xptr2 <- getXPtr(mesh2)
      ixptr <- private[[".CGALmesh"]]$intersection(xptr2)
      cgalMesh$new(clean = ixptr)
    },
    
    #' @description Check whether the mesh is closed.
    #' @return A Boolean value, whether the mesh is closed.
    "isClosed" = function() {
      private[[".CGALmesh"]]$isClosed()
    },

    #' @description Isotropic remeshing.
    #' @param targetEdgeLength positive number, the target edge length of the
    #'   remeshed mesh
    #' @param iterations number of iterations, a positive integer
    #' @param relaxSteps number of relaxation steps, a positive integer
    #' @return The modified \code{cgalMesh} object, invisibly.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' mesh <- cgalMesh$new(HopfTorusMesh(nu = 80, nv = 50))
    #' mesh$isotropicRemeshing(targetEdgeLength = 0.7)
    #' # squared norms of the vertices
    #' normsq <- apply(mesh$getVertices(), 1L, crossprod)
    #' # fair the region where the squared norm is > 19
    #' mesh$fair(which(normsq > 19))
    #' # plot
    #' mesh$computeNormals()
    #' rmesh <- mesh$getMesh()
    #' open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(0, 0)
    #' shade3d(rmesh, color = "maroon")
    #' wire3d(rmesh)}
    "isotropicRemeshing" = function(
      targetEdgeLength, iterations = 1, relaxSteps = 1
    ) {
      stopifnot(isPositiveNumber(targetEdgeLength))
      stopifnot(isStrictPositiveInteger(iterations))
      stopifnot(isStrictPositiveInteger(relaxSteps))
      private[[".CGALmesh"]]$isotropicRemeshing(
        targetEdgeLength, as.integer(iterations), as.integer(relaxSteps)        
      )
      invisible(self)
    },
    
    #' @description Check whether the mesh is outward oriented. The mesh must 
    #'   be triangle.
    #' @return A Boolean value, whether the mesh is outward oriented.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(tetrahedron3d())
    #' mesh$isOutwardOriented() # TRUE
    #' mesh$reverseOrientation()
    #' mesh$isOutwardOriented() # FALSE
    "isOutwardOriented" = function() {
      private[[".CGALmesh"]]$isOutwardOriented()
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

    #' @description Check whether the mesh is valid.
    #' @return A Boolean value, whether the mesh is valid.
    "isValid" = function() {
      private[[".CGALmesh"]]$isValid()
    },

    #' @description Check whether the mesh is valid.
    #' @return A Boolean value, whether the mesh is valid.
    "isValidFaceGraph" = function() {
      private[[".CGALmesh"]]$isValidFaceGraph()
    },
    
    #' @description Check whether the mesh is valid.
    #' @return A Boolean value, whether the mesh is valid.
    "isValidHalfedgeGraph" = function() {
      private[[".CGALmesh"]]$isValidHalfedgeGraph()
    },

    #' @description Check whether the mesh is valid.
    #' @return A Boolean value, whether the mesh is valid.
    "isValidPolygonMesh" = function() {
      private[[".CGALmesh"]]$isValidPolygonMesh()
    },

    #' @description Performs the Loop subdivision and deformation.
    #' @param iterations number of iterations
    #' @return The modified reference mesh, invisibly.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' hopfMesh <- HopfTorusMesh(nu = 80, nv = 40)
    #' mesh <- cgalMesh$new(hopfMesh)
    #' mesh$LoopSubdivision(iterations = 2)
    #' mesh$computeNormals()
    #' rmesh <- mesh$getMesh()
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 800, 400))
    #' mfrow3d(1, 2)
    #' view3d(0, 0, zoom = 0.9)
    #' shade3d(hopfMesh, color = "gold")
    #' wire3d(hopfMesh, color = "black")
    #' next3d()
    #' view3d(0, 0, zoom = 0.9)
    #' shade3d(rmesh, color = "gold")
    #' wire3d(rmesh, color = "black")}
    "LoopSubdivision" = function(iterations = 1) {
      stopifnot(isStrictPositiveInteger(iterations))
      private[[".CGALmesh"]]$LoopSubdivision(as.integer(iterations))
      invisible(self)
    },
    
    #' @description Merge the mesh and another mesh.
    #' @param mesh2 a \code{cgalMesh} object
    #' @return The updated reference mesh, invisibly.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' mesh1 <- cgalMesh$new(sphereMesh())
    #' mesh1$assignFaceColors("red")
    #' mesh2 <- cgalMesh$new(sphereMesh(x = 3))
    #' mesh2$assignFaceColors("blue")
    #' mesh1$merge(mesh2)
    #' rmesh <- mesh1$getMesh()
    #' open3d(windowRect = c(50, 50, 562, 562))
    #' shade3d(rmesh, meshColor = "faces")}
    "merge" = function(mesh2) {
      stopifnot(isCGALmesh(mesh2))
      mesh2XPtr <- getXPtr(mesh2)
      . <- private[[".CGALmesh"]]$merge(mesh2XPtr)
      invisible(self)
    },
    
    #' @description Reorient the connected components of the mesh in order that 
    #' it bounds a volume. The mesh must be triangle.
    #' @return The modified \code{cgalMesh} object, invisibly. \strong{WARNING}: 
    #'   even if you store the result in a new variable, the original mesh is 
    #'   modified. 
    #' @examples 
    #' # two disjoint tetrahedra ####
    #' vertices <- rbind(
    #'   c(0, 0, 0),
    #'   c(2, 2, 0),
    #'   c(2, 0, 2),
    #'   c(0, 2, 2),
    #'   c(3, 3, 3),
    #'   c(5, 5, 3),
    #'   c(5, 3, 5),
    #'   c(3, 5, 5)
    #' )
    #' faces <- rbind(
    #'   c(3, 2, 1),
    #'   c(3, 4, 2),
    #'   c(1, 2, 4),
    #'   c(4, 3, 1),
    #'   c(5, 6, 7),
    #'   c(6, 8, 7),
    #'   c(8, 6, 5),
    #'   c(5, 7, 8)
    #' )
    #' mesh <- cgalMesh$new(vertices = vertices, faces = faces)
    #' mesh$boundsVolume() # FALSE
    #' mesh$orientToBoundVolume()
    #' mesh$boundsVolume() # TRUE
    "orientToBoundVolume" = function() {
      private[[".CGALmesh"]]$orientToBoundVolume()
      invisible(self)
    },
    
    #' @description Remove self-intersections (experimental). The mesh must 
    #'   be triangle.
    #' @return The modified \code{cgalMesh} object, invisibly.
    "removeSelfIntersections" = function() {
      private[[".CGALmesh"]]$removeSelfIntersections()
      invisible(self)
    },
    
    #' @description Reverse the orientation of the faces of the mesh.
    #' @return The modified \code{cgalMesh} object, invisibly. \strong{WARNING}: 
    #'   even if you store the result in a new variable, the original mesh is 
    #'   modified. 
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(tetrahedron3d())
    #' mesh$isOutwardOriented() # TRUE
    #' mesh$reverseOrientation()
    #' mesh$isOutwardOriented() # FALSE
    "reverseOrientation" = function() {
      private[[".CGALmesh"]]$reverseFaceOrientations()
      invisible(self)
    },
    
    #' @description Check whether the mesh self-intersects. The mesh must be 
    #'   triangle.
    #' @return A Boolean value, whether the mesh self-intersects.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(dodecahedron3d())
    #' mesh$selfIntersects()
    "selfIntersects" = function() {
      private[[".CGALmesh"]]$doesSelfIntersect()
    },

    #' @description Returns edges considered to be sharp according to the given 
    #'   angle bound.
    #' @param angleBound angle bound in degrees; an edge whose corresponding 
    #'   dihedral angle is smaller than this bound is considered as sharp 
    #' @return An integer matrix with three columns: \code{"edge"}, an edge 
    #'   index, and \code{"v1"} and \code{"v2"}, the vertex indices of this 
    #'   edge.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' # astroidal ellipsoid
    #' f <- function(u, v) {
    #'   rbind(
    #'     cos(u)^3 * cos(v)^3,
    #'     sin(u)^3 * cos(v)^3,
    #'     sin(v)^3
    #'   )
    #' }
    #' rmesh <- parametricMesh(
    #'   f, urange = c(0, 2*pi), vrange = c(0, 2*pi), 
    #'   periodic = c(TRUE, TRUE), nu = 120, nv = 110
    #' )
    #' mesh <- cgalMesh$new(rmesh)
    #' sharpEdges <- mesh$sharpEdges(30)
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
    #' shade3d(addNormals(rmesh), color = "chartreuse")
    #' plotEdges(
    #'   mesh$getVertices(), sharpEdges[, c("v1", "v2")], 
    #'   edgesAsTubes = FALSE, lwd = 5, verticesAsSpheres = FALSE
    #' )}
    "sharpEdges" = function(angleBound) {
      stopifnot(isNumber(angleBound))
      private[[".CGALmesh"]]$sharpEdges(angleBound - 180)
    },
    
    #' @description Performs the 'Sqrt3' subdivision and deformation. The mesh 
    #'   must be triangle.
    #' @param iterations number of iterations
    #' @return The modified reference mesh, invisibly.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' hopfMesh <- HopfTorusMesh(nu = 80, nv = 40)
    #' mesh <- cgalMesh$new(hopfMesh)
    #' mesh$Sqrt3Subdivision(iterations = 2)
    #' mesh$computeNormals()
    #' rmesh <- mesh$getMesh()
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 800, 400))
    #' mfrow3d(1, 2)
    #' view3d(0, 0, zoom = 0.9)
    #' shade3d(hopfMesh, color = "cyan")
    #' wire3d(hopfMesh, color = "black")
    #' next3d()
    #' view3d(0, 0, zoom = 0.9)
    #' shade3d(rmesh, color = "cyan")
    #' wire3d(rmesh, color = "black")}
    "Sqrt3Subdivision" = function(iterations = 1) {
      stopifnot(isStrictPositiveInteger(iterations))
      private[[".CGALmesh"]]$Sqrt3Subdivision(as.integer(iterations))
      invisible(self)
    },
    
    #' @description Subtract a mesh. Both meshes must be triangle. Face 
    #'   properties of the two meshes are copied to the new mesh. 
    #'  \strong{WARNING}: this modifies the reference mesh and \code{mesh2}. 
    #' @param mesh2 a \code{cgalMesh} object
    #' @return A \code{cgalMesh} object, the difference between the reference 
    #'   mesh and \code{mesh2}. Both the reference mesh and \code{mesh2} are 
    #'   modified: they are corefined.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' # take two cubes
    #' rglmesh1 <- cube3d()
    #' rglmesh2 <- translate3d(cube3d(), 1, 1, 1)
    #' mesh1 <- cgalMesh$new(rglmesh1)
    #' mesh2 <- cgalMesh$new(rglmesh2)
    #' # the two meshes must be triangle
    #' mesh1$triangulate()
    #' mesh2$triangulate()
    #' # assign colors
    #' mesh1$assignFaceColors("red")
    #' mesh2$assignFaceColors("navy")
    #' # difference
    #' mesh <- mesh1$subtract(mesh2)
    #' rglmesh <- mesh$getMesh()
    #' # extract edges for plotting
    #' extEdges <- exteriorEdges(mesh$getEdges())
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglmesh, meshColor = "faces")
    #' plotEdges(mesh$getVertices(), extEdges)
    #' shade3d(rglmesh2, color = "cyan", alpha = 0.2)}
    "subtract" = function(mesh2) {
      stopifnot(isCGALmesh(mesh2))
      xptr2 <- getXPtr(mesh2)
      xptrs <- private[[".CGALmesh"]]$subtract(xptr2)
      if(length(xptrs) == 3L) {
        private[[".CGALmesh"]] <- 
          CGALmesh$new(xptrs[["mesh1"]])
        mesh2[[".__enclos_env__"]][["private"]][[".CGALmesh"]] <- 
          CGALmesh$new(xptrs[["mesh2"]])
        cgalMesh$new(clean = xptrs[["dmesh"]])
      } else {
        invisible(self)
      }
    },
    
    #' @description Triangulate mesh.
    #' @return The modified \code{cgalMesh} object, invisibly. \strong{WARNING}: 
    #'   even if you store the result in a new variable, the original mesh is 
    #'   modified (see the example). You may want to triangulate a copy of the 
    #'   mesh; see the \code{copy} method.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())
    #' mesh$isTriangle() # FALSE
    #' # warning: triangulating the mesh modifies it
    #' mesh$triangulate()
    #' mesh$isTriangle() # TRUE
    "triangulate" = function() {
      private[[".CGALmesh"]]$triangulate()
      invisible(self)
    },

    #' @description Union with another mesh. Both meshes must be triangle. Face 
    #'   properties of the two united meshes are copied to the union mesh. 
    #'  \strong{WARNING}: this modifies the reference mesh and \code{mesh2}. 
    #' @param mesh2 a \code{cgalMesh} object
    #' @return A \code{cgalMesh} object, the union of the reference mesh with 
    #'   \code{mesh2}. Both the reference mesh and \code{mesh2} are modified: 
    #'   they are corefined.
    #' @examples 
    #' \donttest{library(cgalMeshes)
    #' library(rgl)
    #' # take two cubes
    #' rglmesh1 <- cube3d()
    #' rglmesh2 <- translate3d(cube3d(), 1, 1, 1)
    #' mesh1 <- cgalMesh$new(rglmesh1)
    #' mesh2 <- cgalMesh$new(rglmesh2)
    #' # the two meshes must be triangle
    #' mesh1$triangulate()
    #' mesh2$triangulate()
    #' # assign a color to the faces; they will be retrieved in the union
    #' mesh1$assignFaceColors("yellow")
    #' mesh2$assignFaceColors("navy")
    #' # union
    #' umesh <- mesh1$union(mesh2)
    #' rglumesh <- umesh$getMesh()
    #' # extract edges for plotting
    #' extEdges <- exteriorEdges(umesh$getEdges())
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglumesh, meshColor = "faces")
    #' plotEdges(umesh$getVertices(), extEdges)}
    "union" = function(mesh2) {
      stopifnot(isCGALmesh(mesh2))
      xptr2 <- getXPtr(mesh2)
      xptrs <- private[[".CGALmesh"]]$Union(xptr2)
      if(length(xptrs) == 3L) {
        private[[".CGALmesh"]] <- 
          CGALmesh$new(xptrs[["mesh1"]])
        mesh2[[".__enclos_env__"]][["private"]][[".CGALmesh"]] <- 
          CGALmesh$new(xptrs[["mesh2"]])
      }
      cgalMesh$new(clean = xptrs[["umesh"]])
    },
    
    #' @description Compute the volume of the mesh. The mesh must be closed,
    #'   triangle, and must not self-intersect.
    #' @return A number, the mesh volume.
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(cube3d())$triangulate()
    #' mesh$volume()
    "volume" = function() {
      private[[".CGALmesh"]]$volume()
    },

    #' @description Locate points with respect to a closed triangle mesh.
    #' @param points a numeric matrix with three columns
    #' @return An integer vector taking values \code{-1} for outside, \code{1} 
    #'   for inside, and \code{0} if the point is on the boundary.
    #' @examples 
    #' library(cgalMeshes)
    #' mesh <- cgalMesh$new(sphereMesh())
    #' pt1 <- c(0, 0, 0) # inside
    #' pt2 <- c(2, 0, 0) # outside
    #' mesh$whereIs(rbind(pt1, pt2))
    "whereIs" = function(points) {
      if(!is.matrix(points)) {
        points <- rbind(points)
      }
      stopifnot(ncol(points) == 3L)
      storage.mode(points) <- "double"
      if(anyNA(points)) {
        stop("Found missing values.")
      }
      private[[".CGALmesh"]]$whereIs(t(points))
    },
    
    #' @description Write mesh to a file.
    #' @param filename path to the file to be written, with extension 
    #'   \code{off} or \code{ply}
    #' @param precision a positive integer, the desired number of decimal 
    #'   places
    #' @param comments for \code{ply} extension only, a string to be included 
    #'   in the header of the PLY file
    #' @return Nothing, just writes a file.
    "writeMeshFile" = function(filename, precision = 17, comments = "") {
      stopifnot(isString(filename))
      stopifnot(isPositiveInteger(precision))
      stopifnot(isString(comments))
      #stopifnot(isBoolean(binary))
      filename <- path.expand(filename)
      normals <- self$getNormals()
      if(!is.null(normals)) {
        normals <- t(normals)
      }
      fcolors <- self$getFaceColors()
      if(!is.null(fcolors)) {
        fcolors <- tryCatch({
          col2rgb(fcolors)
        }, error = function(e) {
          warning(
            "Invalid face colors found - skipping."
          )
          NULL
        })
      }
      vcolors <- self$getVertexColors()
      if(!is.null(vcolors)) {
        vcolors <- tryCatch({
          col2rgb(vcolors)
        }, error = function(e) {
          warning(
            "Invalid vertex colors found - skipping."
          )
          NULL
        })
      }
      private[[".CGALmesh"]]$writeFile(
        filename, as.integer(precision), FALSE, comments, 
        normals, fcolors, vcolors
      )
    }
    
  )
  
)