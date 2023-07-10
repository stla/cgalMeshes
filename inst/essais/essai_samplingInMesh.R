library(cgalMeshes)
cyclide <- cyclideMesh(a = 97, c = 32, mu = 57)
mesh <- cgalMesh$new(cyclide)
sims <- mesh$sampleInMesh(100)
mean(sims[, 3] > 0)
summary(sims)

library(rgl)
shade3d(cyclide, color = "green", alpha = 0.2)
points3d(sims)

