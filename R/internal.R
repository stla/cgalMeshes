isAtomicVector <- function(x) {
  is.atomic(x) && is.vector(x)
}

isBoolean <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

#' @importFrom gmp is.bigq is.matrixZQ
#' @importFrom data.table uniqueN
#' @noRd
checkMesh <- function(vertices, faces, gmp, aslist) {
  if(gmp) {
    if(!is.matrixZQ(vertices) || ncol(vertices) != 3L || !is.bigq(vertices)) {
      stop("The `vertices` argument must be a `bigq` matrix with three columns.")
    }
    vertices <- as.character(vertices)
  } else {
    if(!is.matrix(vertices) || ncol(vertices) != 3L) {
      stop("The `vertices` argument must be a matrix with three columns.")
    }
    stopifnot(is.numeric(vertices))
    storage.mode(vertices) <- "double"
  }
  if(anyNA(vertices)) {
    stop("Found missing values in `vertices`.")
  }
  homogeneousFaces <- FALSE
  isTriangle       <- FALSE
  toRGL            <- FALSE
  if(is.matrix(faces)) {
    if(ncol(faces) < 3L) {
      stop("Faces must be given by at least three indices.")
    }
    storage.mode(faces) <- "integer"
    if(anyNA(faces)) {
      stop("Found missing values in `faces`.")
    }
    if(any(faces < 1L)) {
      stop("Faces cannot contain indices lower than 1.")
    }
    if(any(faces > nrow(vertices))) {
      stop("Faces cannot contain indices higher than the number of vertices.")
    }
    homogeneousFaces <- ncol(faces)
    if(homogeneousFaces %in% c(3L, 4L)) {
      isTriangle <- homogeneousFaces == 3L
      toRGL <- homogeneousFaces
    }
    if(aslist) {
      faces <- lapply(1L:nrow(faces), function(i) faces[i, ] - 1L)
    } else {
      faces <- t(faces - 1L)
    }
  }else if(is.list(faces)) {
    check <- all(vapply(faces, isAtomicVector, logical(1L)))
    if(!check) {
      stop("The `faces` argument must be a list of integer vectors.")
    }
    check <- any(vapply(faces, anyNA, logical(1L)))
    if(check) {
      stop("Found missing values in `faces`.")
    }
    faces <- lapply(faces, function(x) as.integer(x) - 1L)
    sizes <- lengths(faces)
    if(any(sizes < 3L)) {
      stop("Faces must be given by at least three indices.")
    }
    check <- any(vapply(faces, function(f) {
      any(f < 0L) || any(f >= nrow(vertices))
    }, logical(1L)))
    if(check) {
      stop(
        "Faces cannot contain indices lower than 1 or higher than the ",
        "number of vertices."
      )
    }
    usizes <- uniqueN(sizes)
    if(usizes == 1L) {
      homogeneousFaces <- sizes[1L]
      isTriangle <- homogeneousFaces == 3L
      if(homogeneousFaces %in% c(3L, 4L)) {
        toRGL <- homogeneousFaces
      }
    }else if(usizes == 2L && all(sizes %in% c(3L, 4L))) {
      toRGL <- 34L
    }
  } else {
    stop("The `faces` argument must be a list or a matrix.")
  }
  list(
    "vertices"         = t(vertices),
    "faces"            = faces,
    "homogeneousFaces" = homogeneousFaces,
    "isTriangle"       = isTriangle,
    "toRGL"            = toRGL
  )
}

getVFT <- function(mesh, beforeCheck = FALSE) {
  transposed <- !beforeCheck
  i0 <- as.integer(transposed)
  if(inherits(mesh, "mesh3d")) {
    triangles <- mesh[["it"]]
    if(!is.null(triangles)) {
      triangles <- lapply(1L:ncol(triangles), function(i) triangles[, i] - i0)
    }
    quads <- mesh[["ib"]]
    isTriangle <- is.null(quads)
    if(!isTriangle) {
      quads <- lapply(1L:ncol(quads), function(i) quads[, i] - i0)
    }
    faces <- c(triangles, quads)
    vertices <- mesh[["vb"]][-4L, ]
    if(!transposed) {
      vertices <- t(vertices)
    }
    rmesh <- list("vertices" = vertices, "faces" = faces)
  } else if(inherits(mesh, "cgalMesh")) {
    isTriangle <- attr(mesh, "toRGL") == 3L
    vertices <- mesh[["vertices"]]
    if(transposed) {
      vertices <- t(vertices)
    }
    faces <- mesh[["faces"]]
    if(is.matrix(faces)) {
      faces <- lapply(1L:nrow(faces), function(i) faces[i, ] - i0)
    }else if(!beforeCheck) {
      faces <- lapply(faces, function(face) face - 1L)
    }
    rmesh <- list("vertices" = vertices, "faces" = faces)
  } else if(is.list(mesh)) {
    rmesh <-
      checkMesh(mesh[["vertices"]], mesh[["faces"]], gmp = FALSE, aslist = TRUE)
    isTriangle <- rmesh[["isTriangle"]]
    if(beforeCheck) {
      rmesh[["vertices"]] <- t(rmesh[["vertices"]])
      rmesh[["faces"]] <- lapply(rmesh[["faces"]], function(face) face + 1L)
    }
  } else {
    stop("Invalid `mesh` argument.", call. = FALSE)
  }
  list("rmesh" = rmesh, "isTriangle" = isTriangle)
}
