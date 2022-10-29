rho <- sqrt((5 - sqrt(5))/10)
vs1 <- t(vapply(0:4, function(i){
  c(rho*cos(2*i*pi/5), rho*sin(2*i*pi/5), 0.1)
}, numeric(3L)))
vs2 <- t(vapply(0:4, function(i){
  c(rho*cos(2*i*pi/5), rho*sin(2*i*pi/5), -0.1)
}, numeric(3L)))
R <- sqrt((25 - 11*sqrt(5))/10)
vs3 <- t(vapply(0:4, function(i){
  c(R*cos(2*i*pi/5 + pi/5), R*sin(2*i*pi/5 + pi/5), 0.1)
}, numeric(3L)))
vs4 <- t(vapply(0:4, function(i){
  c(R*cos(2*i*pi/5 + pi/5), R*sin(2*i*pi/5 + pi/5), -0.1)
}, numeric(3L)))
vertices <- rbind(vs1, vs2, vs3, vs4)

triangles <- list(
  c(15L, 1L, 11L),
  c(11L, 2L, 12L),
  c(12L, 3L, 13L),
  c(13L, 4L, 14L),
  c(14L, 5L, 15L)
)
triangles <- c(
  triangles,
  lapply(triangles, function(x) x + 5L)
)
pentagons <- list(
  11L:15L,
  16L:20L
)
rectangles <- list(
  c(1L, 6L, 16L, 11L),
  c(16L, 11L, 2L, 7L),
  c(2L, 7L, 17L, 12L),
  c(17L, 12L, 3L, 8L),
  c(3L, 8L, 18L, 13L),
  c(18L, 13L, 4L, 9L),
  c(4L, 9L, 19L, 14L),
  c(19L, 14L, 5L, 10L),
  c(5L, 10L, 20L, 15L),
  c(20L, 15L, 1L, 6L)
)
faces <- c(triangles, rectangles, pentagons)

pentagrammicPrism <- list("vertices" = vertices, "faces" = faces)
