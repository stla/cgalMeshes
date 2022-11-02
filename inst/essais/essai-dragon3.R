library(SurfaceReconstruction)
library(cgalMeshes)

rglDragon <- AFSreconstruction(StanfordDragon)
dragon <- cgalMesh$new(rglDragon)

rglDragon <- dragon$getMesh()

vs <- dragon$vertices()

index <- which.max(vs[, 1L])
geoDists <- dragon$geoDists(index)
x <- seq(min(geoDists), max(geoDists), length.out = 75L)

# anim
for(i in seq_along(x)){
  colors <- rep("forestgreen", nrow(vs))
  red <- geoDists <= x[i]
  colors[red] <- "darkslategrey"
  rglDragon[["material"]] <- list("color" = colors)
  # plot
  open3d(windowRect = 50 + c(0, 0, 512, 512))
  view3d(-25, 0, zoom = 0.75)
  shade3d(rglDragon)
  snapshot3d(sprintf("zzpic%03d.png", i), webshot = FALSE)
  close3d()
}

command <- "convert -delay 1x15 -duplicate 1,-2-1 zzpic*.png dragon-venom.gif"
system(command)
file.remove(Sys.glob("zzpic*.png"))
