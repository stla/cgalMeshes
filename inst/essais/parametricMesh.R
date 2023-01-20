f <- Vectorize(function(u, v) {
  c(0, u, v) # or rbind(0, u, v)
})

u_ <- 1:4
v_ <- 1:5
Grid <- expand.grid(U = u_, V = v_)
varray <- 
  with(Grid, array(f(U, V), dim = c(3, length(u_), length(v_))))

varray2 <- aperm(varray, c(1L, 3L, 2L))
matrix(varray2, nrow = 3L, ncol = 20L)

# f, umin, umax, vmin, vmax, uperiodic, vperiodic

#n = 3; ParametricPlot3D[{, r Sin[phi] + r^(2n - 1) Sin[(2n - 1) phi]/(2n - 1), 2 r^n Cos[n phi]/n, EdgeForm[]}, {phi, 0, 2Pi}, {r, 0, 1.3}, PlotPoints -> {181, 20}, ViewPoint -> {0, 0, 1}, PlotRange -> All]

library(rgl)
n <- 3
f <- Vectorize(function(phi, r) {
  c(
    r*cos(phi) - r^(2*n-1)*cos((2*n-1)*phi)/(2*n-1),
    r*sin(phi) + r^(2*n-1)*sin((2*n-1)*phi)/(2*n-1),
    2*r^n*cos(n*phi)/n
  )
})

nu <- 140
nv <- 150
u_ <- seq(0, 2*pi, length.out = nu+1)[-1]
v_ <- seq(0.1, 1.3, length.out = nv)
Grid <- expand.grid(U = u_, V = v_)
varray <- 
  with(Grid, array(f(U, V), dim = c(3, nu, nv)))
varray2 <- aperm(varray, c(1L, 3L, 2L))
vs <- matrix(varray2, nrow = 3L, ncol = nu*nv)

tris <- cgalMeshes:::meshTopology(nu, nv)

rmesh <- tmesh3d(
  vertices    = vs,
  indices     = tris,
  homogeneous = FALSE
)

shade3d(addNormals(rmesh), color = "yellow")

parametricMesh <- function(
  f, urange, vrange, periodic = c(FALSE, FALSE), nu, nv
) {
 uperiodic <- periodic[1L]
 vperiodic <- periodic[2L]
 if(uperiodic) {
   u_ <- seq(urange[1L], urange[2L], length.out = nu+1)[-1L]
 } else {
   u_ <- seq(urange[1L], urange[2L], length.out = nu)
 }
 if(vperiodic) {
   v_ <- seq(vrange[1L], vrange[2L], length.out = nv+1)[-1L]
 } else {
   v_ <- seq(vrange[1L], vrange[2L], length.out = nv)
 }
 Grid <- expand.grid(U = u_, V = v_)
 varray <- with(Grid, array(f(U, V), dim = c(3, nu, nv)))
 varray2 <- aperm(varray, c(1L, 3L, 2L))
 vs <- matrix(varray2, nrow = 3L, ncol = nu*nv)
 tris <- meshTopology(nu, nv, uperiodic, vperiodic)
 tmesh3d(
   vertices    = vs,
   indices     = tris,
   homogeneous = FALSE
 )
}

