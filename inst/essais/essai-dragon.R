library(SurfaceReconstruction)
library(cgalMeshes)

rglDragon <- AFSreconstruction(StanfordDragon)
dragon <- cgalMesh$new(rglDragon)
vs <- dragon$vertices()

index <- which.max(vs[, 1L])
geoDists <- dragon$geoDists(index)
geoDists <- geoDists / max(geoDists)

rglDragon <- dragon$getMesh()

fcolor <- colorRamp(viridisLite::plasma(200L))
colors <- fcolor(geoDists)
colors <- rgb(colors[, 1L], colors[, 2L], colors[, 3L], maxColorValue = 255)
rglDragon[["material"]] <- list("color" = colors)
# plot
library(rgl)
open3d(windowRect = 50 + c(0, 0, 512, 512))
view3d(-30, 0, zoom = 0.7)
shade3d(rglDragon)
snapshot3d("dragon.png", webshot = FALSE)
