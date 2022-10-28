library(cgalMeshes)

cube <- rgl::cube3d()
vs <- cube$vb[-4L, ]
fs <- lapply(1:ncol(cube$ib), function(i) cube$ib[, i]-1L)

cm <- cgalMeshes:::CGALmesh$new(vs, fs, TRUE)
