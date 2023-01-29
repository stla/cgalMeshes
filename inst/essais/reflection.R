library(cgalMeshes)
library(rgl)

A <- c(1, 0, 0) #/ sqrt(2)
B <- c(0, 1, 0) #/ sqrt(2)
C <- c(0, 0, 1)#/ sqrt(2)

n <- rgl:::xprod(C-A, B-A)
n <- n / sqrt(c(crossprod(n)))
R <- diag(3) - 2 * tcrossprod(n)#/c(crossprod(n))

#a <- n[1]; b <- n[2]; c <- n[3]
offset <- c(crossprod(n, B)) 
M <- rbind(
  cbind(R, -2*offset*n), 
  c(0, 0, 0, 1)
)
M <- t(M)
mesh1 <- sphericalTriangle(A, B, C, iterations = 4)
rmesh1 <- mesh1$getMesh()
#rmesh1$vb[-4,] <- rmesh1$vb[-4,] / sqrt(1.8)
# xx <- -2*offset*n
# rmesh2 <- translate3d(transform3d(rmesh1, R), xx[1], xx[2], xx[3])
# rmesh2 <- transform3d(translate3d(rmesh1, xx[1], xx[2], xx[3]), R)
rmesh2 <- transform3d(rmesh1, M)
rmesh2$normals <- -rmesh2$normals


shade3d(rmesh1, color = "red")
shade3d(rmesh2, color = "green")
