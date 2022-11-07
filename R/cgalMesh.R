getXPtr <- function(cMesh){
  cMesh[[".__enclos_env__"]][["private"]][[".CGALmesh"]][["xptr"]]
}

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
    #' @examples 
    #' library(cgalMeshes)
    #' meshFile <- system.file(
    #'   "extdata", "bigPolyhedron.off", package = "cgalMeshes"
    #' )
    #' mesh <- cgalMesh$new(meshFile)
    #' rglmesh <- mesh$getMesh(normals = FALSE)
    #' \donttest{library(rgl)
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglmesh, color = "tomato")
    #' plotEdges(
    #'   mesh$vertices(), mesh$edges(), color = "darkred"
    #' )}
    "initialize" = function(
      mesh, vertices, faces, clean = TRUE
    ){
      # one can also initialize from an external pointer, but 
      # this is hidden to the user
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
    
    #' @description Clip mesh to the volume bounded by another mesh. 
    #'   \strong{WARNING}: the reference mesh is then replaced by its 
    #'   clipped version.
    #'
    #' @param clipper a \code{cgalMesh} object; it must represent a closed 
    #'   triangle mesh which doesn't self-intersect
    #' @param clipVolume Boolean, whether the clipping has to be done on the 
    #'   volume bounded by the reference mesh rather than on its surface (i.e. 
    #'   the reference mesh will be kept closed if it is closed); if 
    #'   \code{TRUE}, the mesh to be clipped must not self-intersect
    #' @return The modified \code{cgalObject}.
    #' @examples 
    #' # cube clipped to sphere ####
    #' library(cgalMeshes)
    #' library(rgl)
    #' mesh    <- cgalMesh$new(cube3d())$triangulate()
    #' clipper <- cgalMesh$new(sphereMesh(r= sqrt(2)))
    #' mesh$clip(clipper, clipVolume = TRUE)
    #' rglmesh <- mesh$getMesh(normals = FALSE)
    #' \donttest{open3d(windowRect = 50 + c(0, 0, 512, 512))
    #' view3d(45, 45, zoom = 0.9)
    #' shade3d(rglmesh, col = "darkorange")}
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
      private[[".CGALmesh"]]$clipMesh(clipperXPtr, clipVolume)
      self
    },

    #' @description Decomposition into connected components.
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
    #'   vertices = isosurface[["vertices"]], faces = isosurface[["triangles"]]
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
    "connectedComponents" = function(triangulate = TRUE) {
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
    #'   cxmesh <- cxparts[[i]]$getMesh(normals = FALSE)
    #'   shade3d(cxmesh, color = colors[i])
    #' }}
    "convexParts" = function(triangulate = TRUE) {
      stopifnot(isBoolean(triangulate))
      xptrs <- private[[".CGALmesh"]]$convexParts(triangulate)
      lapply(xptrs, function(xptr) cgalMesh$new(clean = xptr))
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

    #' @description Get the edges of the mesh.
    #' @return A dataframe with four columns; the first two ones give the 
    #'   vertex indices of each edge (one edge per row), the third one gives 
    #'   the lengths of each edge, and the fourth one gives the dihedral angles 
    #'   in degrees between the two faces adjacent to each edge 
    #' @examples 
    #' library(rgl)
    #' mesh <- cgalMesh$new(dodecahedron3d())
    #' head(mesh$edges())
    "edges" = function() {
      private[[".CGALmesh"]]$edges()
    },

    #' @description Fair a region of the mesh, i.e. make it smooth. The mesh 
    #'   must be triangle. This modifies the reference mesh.
    #' @param indices the indices of the vertices in the region to be faired
    #' @return The modified \code{cgalMesh} object.
    #' @examples 
    #' library(cgalMeshes)
    #' rglHopf <- HopfTorusMesh(nu = 100, nv = 100)
    #' hopf <- cgalMesh$new(rglHopf)
    #' # squared norms of the vertices
    #' normsq <- apply(hopf$vertices(), 1L, crossprod)
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
      self
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
      if(normals) {
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
          }
        }
      } else {
        mesh[["vertices"]] <- t(mesh[["vertices"]])
      }
      mesh
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
    #' rglimesh <- imesh$getMesh(normals = FALSE)
    #' # extract edges for plotting
    #' extEdges <- exteriorEdges(imesh$edges())
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglimesh, color = "red")
    #' plotEdges(imesh$vertices(), extEdges)
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
    
    #' @description Subtract another mesh. Both meshes must be triangle.
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
    #' # difference
    #' mesh <- mesh1$subtract(mesh2)
    #' rglmesh <- mesh$getMesh(normals = FALSE)
    #' # extract edges for plotting
    #' extEdges <- exteriorEdges(mesh$edges())
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglmesh, color = "red")
    #' plotEdges(mesh$vertices(), extEdges)
    #' shade3d(rglmesh2, color = "cyan", alpha = 0.2)}
    "subtract" = function(mesh2) {
      stopifnot(isCGALmesh(mesh2))
      xptr2 <- getXPtr(mesh2)
      dxptr <- private[[".CGALmesh"]]$subtract(xptr2)
      cgalMesh$new(clean = dxptr)
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

    #' @description Union with another mesh. Both meshes must be triangle.
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
    #' # union
    #' umesh <- mesh1$union(mesh2)
    #' rglumesh <- umesh$getMesh(normals = FALSE)
    #' # extract edges for plotting
    #' extEdges <- exteriorEdges(umesh$edges())
    #' # plot
    #' open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.9)
    #' shade3d(rglumesh, color = "red")
    #' plotEdges(umesh$vertices(), extEdges)}
    "union" = function(mesh2) {
      stopifnot(isCGALmesh(mesh2))
      xptr2 <- getXPtr(mesh2)
      uxptr <- private[[".CGALmesh"]]$Union(xptr2)
      cgalMesh$new(clean = uxptr)
    },
    
    #' @description Get the vertices of the mesh.
    #' @return The vertices in a matrix.
    "vertices" = function() {
      t(private[[".CGALmesh"]]$vertices())
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
      private[[".CGALmesh"]]$writeFile(
        filename, as.integer(precision), binary
      )
    }
    
  )
  
)