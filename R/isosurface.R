#' Title
#'
#' @param voxel a three dimensional array from which the contour has to 
#'   be calculated
#' @param x,y,z locations of the grid points at which values in \code{voxel} 
#'   are measured
#' @param isolevel the level at which to construct the isosurface
#' @param sphereCenter,sphereRadius center and radius of a bounding sphere
#' @param angleBound x
#' @param radiusBound x
#' @param distanceBound x  
#'
#' @return A \code{cgalMesh} object.
#' @export
#' @importFrom RNifti readNifti writeAnalyze
#' @importFrom neuroim writeVolume BrainVolume BrainSpace
isocontour <- function(
  voxel, x, y, z, isolevel,
  sphereCenter, sphereRadius,
  angleBound, radiusBound, distanceBound
) {
  nx <- length(x)
  ny <- length(y)
  nz <- length(z)
  ax <- x[nx] - x[1L]
  ay <- y[ny] - y[1L]
  az <- z[nz] - z[1L]
  bvol <- BrainVolume(
    voxel, 
    BrainSpace(
      Dim = c(nx, ny, nz), 
      spacing = c(ax / (nx-1L), ay / (ny-1L), az / (nz-1L)),
      origin = c(0, 0, 0)
    )
  )
  niiFile <- tempfile(fileext = ".nii.gz")
  hdrFile <- tempfile(fileext = ".hdr")
  writeVolume(bvol, niiFile)
  img <- readNifti(niiFile)
  . <- writeAnalyze(img, hdrFile)
  xptr <- Isomesh(
    hdrFile, isolevel, 
    sphereCenter + c(ax, ay, az)/2, 
    sphereRadius,
    angleBound, radiusBound, distanceBound
  )
  cgalMesh$new(clean = xptr)
}