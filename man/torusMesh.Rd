% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/someMeshes.R
\name{torusMesh}
\alias{torusMesh}
\title{Torus mesh}
\usage{
torusMesh(
  R,
  r,
  p1 = NULL,
  p2 = NULL,
  p3 = NULL,
  nu = 50,
  nv = 30,
  normals = TRUE,
  conformal = FALSE
)
}
\arguments{
\item{R, r}{major and minor radii, positive numbers; \code{R} is 
ignored if \code{p1}, \code{p2} and \code{p3} are given}

\item{p1, p2, p3}{three points or \code{NULL}; if not \code{NULL}, 
the function returns a mesh of the torus whose centerline passes 
through these three points and with minor radius \code{r}; if 
\code{NULL}, the torus has equatorial plane z=0 and the 
z-axis as revolution axis}

\item{nu, nv}{numbers of subdivisions, integers (at least 3)}

\item{normals}{a Boolean value, whether to compute the normals of the mesh}

\item{conformal}{a Boolean value, whether to use a conformal 
parameterization of the torus (with aspect ratio 
\code{sqrt(R^2/r^2-1)})}
}
\value{
A triangle \strong{rgl} mesh (class \code{mesh3d}).
}
\description{
Triangle mesh of a torus.
}
\examples{
library(cgalMeshes)
library(rgl)
mesh <- torusMesh(R = 3, r = 1)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(0, 0, zoom = 0.75)
shade3d(mesh, color = "green")
wire3d(mesh)

# Villarceau circles ####
Villarceau <- function(beta, theta0, phi) {
  c(
    cos(theta0 + beta) * cos(phi),
    sin(theta0 + beta) * cos(phi),
    cos(beta) * sin(phi)
  ) / (1 - sin(beta) * sin(phi))
}
ncircles <- 30
if(require("randomcoloR")) {
  colors <- 
    randomColor(ncircles, hue = "random", luminosity = "dark")
} else {
  colors <- rainbow(ncircles)
}
theta0_ <- seq(0, 2*pi, length.out = ncircles+1)[-1L]
phi <- 0.7
\donttest{open3d(windowRect = 50 + c(0, 0, 512, 512), zoom = 0.8)
for(i in seq_along(theta0_)) {
  theta0 <- theta0_[i]
  p1 <- Villarceau(0, theta0, phi)
  p2 <- Villarceau(2, theta0, phi)
  p3 <- Villarceau(4, theta0, phi)
  rmesh <- torusMesh(r = 0.05, p1 = p1, p2 = p2, p3 = p3)
  shade3d(rmesh, color = colors[i])
}}
}
\references{
See J.M. Sullivan, 
\href{https://static1.bridgesmathart.org/2011/cdrom/proceedings/134/paper_134.pdf}{Conformal Tiling on a Torus}
(2011) for information about the conformal parameterization.
}
