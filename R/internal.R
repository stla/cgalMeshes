isAtomicVector <- function(x) {
  is.atomic(x) && is.vector(x)
}

isBoolean <- function(x) {
  is.logical(x) && length(x) == 1L && !is.na(x)
}

isString <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x)
}

isStringVector <- function(x) {
  is.character(x) && !anyNA(x)
}

isFilename <- function(x) {
  isString(x) && file.exists(x)
}

isNumber <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x)
}

isPositiveNumber <- function(x) {
  isNumber(x) && x > 0
}

isNonNegativeNumber <- function(x) {
  isNumber(x) && x >= 0
}

isPositiveInteger <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && floor(x) == x
}

isStrictPositiveInteger <- function(x) {
  isPositiveInteger(x) && x > 0
}

isVector3 <- function(x) {
  is.numeric(x) && length(x) == 3L && !anyNA(x)
}

#' @importFrom R6 is.R6
#' @noRd
isCGALmesh <- function(x) {
  is.R6(x) && inherits(x, "cgalMesh")
}

#' @importFrom data.table uniqueN
#' @noRd
checkMesh <- function(vertices, faces, aslist) {
  if(!is.matrix(vertices) || ncol(vertices) != 3L) {
    stop(
      "The vertices must be given as a matrix with three columns.",
      call. = FALSE
    )
  }
  stopifnot(is.numeric(vertices))
  storage.mode(vertices) <- "double"
  if(anyNA(vertices)) {
    stop("Found missing values in the vertices.", call. = FALSE)
  }
  homogeneousFaces <- FALSE
  isTriangle       <- FALSE
  toRGL            <- FALSE
  if(is.matrix(faces)) {
    if(ncol(faces) < 3L) {
      stop("Faces must be given by at least three indices.", call. = FALSE)
    }
    storage.mode(faces) <- "integer"
    if(anyNA(faces)) {
      stop("Found missing values in `faces`.", call. = FALSE)
    }
    if(any(faces < 1L)) {
      stop("Faces cannot contain indices lower than 1.", call. = FALSE)
    }
    if(any(faces > nrow(vertices))) {
      stop(
        "Faces cannot contain indices higher than the number of vertices.",
        call. = FALSE
      )
    }
    dups <- vapply(1L:nrow(faces), function(i) {
      anyDuplicated(faces[i, ])
    }, integer(1L))
    if(any(dups)) {
      stop("A face cannot contain duplicated indices.", call. = FALSE)
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
  } else if(is.list(faces)) {
    check <- all(vapply(faces, isAtomicVector, logical(1L)))
    if(!check) {
      stop("A face must be given as an integer vector.", call. = FALSE)
    }
    faces <- lapply(faces, function(x) as.integer(x) - 1L)
    someNA <- any(vapply(faces, anyNA, logical(1L)))
    if(someNA) {
      stop("A face cannot contain missing values.", call. = FALSE)
    }
    dups <- vapply(faces, function(x) anyDuplicated(x), integer(1L))
    if(any(dups)) {
      stop("A face cannot contain duplicated indices.", call. = FALSE)
    }
    sizes <- lengths(faces)
    if(any(sizes < 3L)) {
      stop("Faces must be given by at least three indices.", call. = FALSE)
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
    } else if(usizes == 2L && all(sizes %in% c(3L, 4L))) {
      toRGL <- 34L
    }
  } else {
    stop("Faces must be given as a list or a matrix.", call. = FALSE)
  }
  list(
    "vertices"         = t(vertices),
    "faces"            = faces,
    "homogeneousFaces" = homogeneousFaces,
    "isTriangle"       = isTriangle,
    "toRGL"            = toRGL
  )
}

getVF <- function(mesh) {
  triangles <- mesh[["it"]]
  if(!is.null(triangles)) {
    triangles <- lapply(1L:ncol(triangles), function(i) triangles[, i])
  }
  quads <- mesh[["ib"]]
  isTriangle <- is.null(quads)
  if(!isTriangle) {
    quads <- lapply(1L:ncol(quads), function(i) quads[, i])
  }
  faces <- c(triangles, quads)
  h <- mesh[["vb"]][4L, ]
  zeros <- h == 0
  h[zeros] <- 1 # these vertices should not be referenced
  vertices <- t(mesh[["vb"]][-4L, ])
  vertices[, 1L] <- vertices[, 1L] / h 
  vertices[, 2L] <- vertices[, 2L] / h 
  vertices[, 3L] <- vertices[, 3L] / h 
  list("vertices" = vertices, "faces" = faces)
}

# Cross product of two 3D vectors.
crossProduct <- function(v, w){
  c(
    v[2L] * w[3L] - v[3L] * w[2L],
    v[3L] * w[1L] - v[1L] * w[3L],
    v[1L] * w[2L] - v[2L] * w[1L]
  )
}

#' @description Get a unit quaternion whose corresponding rotation 
#'   sends `u` to `v`; the vectors `u` and `v` must be normalized.
#' @importFrom onion as.quaternion
#' @noRd   
quaternionFromTo <- function(u, v) {
  re <- sqrt((1 + sum(u*v))/2)
  w <- crossProduct(u, v) / 2 / re
  as.quaternion(c(re, w), single = TRUE)
}

#' @description Get a rotation matrix sending `u` to `v`; 
#'   the vectors `u` and `v` must be normalized.
#' @importFrom onion as.orthogonal
#' @noRd   
rotationFromTo <- function(u, v) {
  as.orthogonal(quaternionFromTo(u, v))
}

circumcircle <- function(p1, p2, p3) {
  # normal of triangle (p1, p2, p3)
  normal <- crossProduct(p2 - p1, p3 - p1)
  squaredNorm <- c(crossprod(normal))
  if(squaredNorm < sqrt(.Machine$double.eps)) {
    stop("The three points are aligned.", call. = FALSE)
  }
  normal <- normal / sqrt(squaredNorm)
  # circumcenter of this triangle
  v12 <- p2 - p1
  v13 <- p3 - p1
  p12 <- p1 + v12/2
  p13 <- p1 + v13/2
  offset <- c(crossprod(p1, normal))
  A <- rbind(normal, v12, v13)
  b <- c(offset, crossprod(p12, v12), crossprod(p13, v13))
  center <- solve(A, b)
  # circumradius
  radius <- sqrt(c(crossprod(p1 - center)))
  #
  list("center" = center, "radius" = radius, "normal" = normal)
}

sph2cart <- function(rho, theta, phi) {
  c(
    rho * cos(theta) * sin(phi),
    rho * sin(theta) * sin(phi),
    rho * cos(phi)
  )
}
