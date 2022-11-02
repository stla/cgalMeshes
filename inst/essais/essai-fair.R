library(cgalMeshes)
rglHopf <- HopfTorusMesh(nu = 200, nv = 200)
hopf <- cgalMesh$new(rglHopf)
# squared norms of the vertices
normsq <- apply(hopf$vertices(), 1L, crossprod)
# fair the region where the squared norm is > 19
indices <- which(normsq > 10)
hopf$fair(indices)
rglHopf_faired <- hopf$getMesh()
# plot
library(rgl)
open3d(windowRect = 50 + c(0, 0, 900, 450))
mfrow3d(1L, 2L)
view3d(0, 0, zoom = 0.7)
shade3d(rglHopf, color = "orangered")
next3d()
view3d(0, 0, zoom = 0.7)
shade3d(rglHopf_faired, color = "orangered")

# animation
rglHopf <- HopfTorusMesh(nu = 200, nv = 200)
x_ <- seq(70, 10, by = -1)
for(i in seq_along(x_)) {
  hopf <- cgalMesh$new(rglHopf)
  # squared norms of the vertices
  normsq <- apply(hopf$vertices(), 1L, crossprod)
  # fair the region where the squared norm is > 19
  indices <- which(normsq > x_[i])
  hopf$fair(indices)
  rglHopf_faired <- hopf$getMesh()
  # plot
  library(rgl)
  open3d(windowRect = 50 + c(0, 0, 512, 512))
  view3d(0, 0, zoom = 0.7)
  shade3d(rglHopf_faired, color = "orangered")
  snapshot3d(sprintf("zzpic%03d.png", i), webshot = FALSE)
  close3d()  
}

command <- "convert -delay 1x11 -duplicate 1,-2-1 zzpic*.png HopfTorus.gif"
system(command)

