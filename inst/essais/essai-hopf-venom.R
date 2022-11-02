library(cgalMeshes)

rglDragon <- HopfTorusMesh(nu = 300, nv = 300)
dragon <- cgalMesh$new(rglDragon)

rglDragon <- dragon$getMesh()

vs <- dragon$vertices()

index <- which.max(vs[, 2L])
geoDists <- dragon$geoDists(index)
x <- seq(min(geoDists), max(geoDists), length.out = 75L)

# anim
for(i in seq_along(x)){
  colors <- rep("orange", nrow(vs))
  red <- geoDists <= x[i]
  colors[red] <- "darkslategrey"
  rglDragon[["material"]] <- list("color" = colors)
  # plot
  open3d(windowRect = 50 + c(0, 0, 512, 512))
  view3d(0, 0, zoom = 0.75)
  shade3d(rglDragon)
  snapshot3d(sprintf("zzpic%03d.png", i), webshot = FALSE)
  close3d()
}

command <- "convert -delay 1x15 -duplicate 1,-2-1 -layers OptimizePlus zzpic*.png hopf-venom.gif"
system(command)
file.remove(Sys.glob("zzpic*.png"))
