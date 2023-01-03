library(rgl)
library(cgalMeshes)

t_ <- seq(0, 1, length.out = 50)

A1 <- c(-2, 0, 0)
B1 <- c(2, 0, 0)
v1 <- B1 - A1 
path1 <- t(sapply(t_, function(t) A1 + t*v1))
cyl1 <- as.tmesh3d(cylinder3d(path1, radius = 1, sides = 60, closed = -2))
cyl1[["material"]] <- list("color" = rep(c("black", "blue"), each = ncol(cyl1[["it"]])/2))
cyl1[["meshColor"]] <- "faces"
A2 <- c(0, -2, 0)
B2 <- c(0, 2, 0)
v2 <- B2 - A2 
path2 <- t(sapply(t_, function(t) A2 + t*v2))
cyl2 <- as.tmesh3d(cylinder3d(path2, radius = 1, sides = 30, closed = -2))
cyl2[["material"]] <- list("color" = rep(c("yellow", "red"), each = ncol(cyl2[["it"]])/2))
cyl2[["meshColor"]] <- "faces"
A3 <- c(0, 0, -3)
B3 <- c(0, 0, 3)
v3 <- B3 - A3 
path3 <- t(sapply(t_, function(t) A3 + t*v3))
cyl3 <- as.tmesh3d(cylinder3d(path3, radius = 1, sides = 60, closed = -2))
cyl3[["material"]] <- list("color" = rep("green", ncol(cyl3[["it"]])))
cyl3[["meshColor"]] <- "faces"

mesh1 <- cgalMesh$new(cyl1)
mesh2 <- cgalMesh$new(cyl2)
mesh3 <- cgalMesh$new(cyl3)

# x= mesh1$doubleclip(mesh2)
# rmesh1 <- x$getMesh()
# cols <- rmesh1$material$color
# cols[cols==""] <- "yellow"
# rmesh1$material$color <- cols
# stop("ok")

# x = mesh2$clip(mesh1,TRUE)
# stop("ooooooooooo")
x = mesh1$clip(mesh2,TRUE)
print(head(x$fimap))
clipper <- cgalMesh$new(clean = x$clipper)
rcols2 <- cols2
for(j in 1:1148) {
  rcols2[x$fmap2[j,1]+1] <- cols2[x$fmap2[j,2]+1]
}
#m12 <- mesh1$doubleclip(mesh2)
#m13 <- m12$doubleclip(mesh3)
rmesh1 <- mesh1$getMesh()
colors <- character(ncol(rmesh1$it))
vb <- t(rmesh1$vb[-4,])
uu = vb[,1]^2 + vb[,3]^2
fff = integer(0)
ggg = integer(0)
ffff = integer(0)
cols1 <- cyl1[["material"]]$color
cols2 <- cyl2[["material"]]$color
kk = 1
x$fimap[x$fimap==0] <- sort(x$xxx[,2], decreasing = TRUE)
for(fi in 1:ncol(rmesh1$it)) {
  fd = x$fimap[fi]
  if((fi-1) %in% x$fmap1[,1]) {
    # ffff <- c(ffff, x$fmap1[which(x$fmap1[,1] == fd), 2])
    # w = which(x$out[,1] == fd)
    colors[fi] = "black"
  } else if((fi-1) %in% x$fmap2[,1]) {
    colors[fi] = cols1[fi]#cols2[1+x$fmap2[which(x$fmap2[,1] == fd), 2]]
  } else if(fi <= 3000) {
    fff <- c(fff, fd)
    ggg <- c(ggg, fi)
    # w = which(x$xxx[,2] == fd)
    # ww = x$xxx[w,1]
    # tr <- rmesh1$it[,fi]
    # test <- abs(uu[tr[2]] - 1) < 0.01 #&& abs(uu[tr[2]] - 1) < 0.01 && abs(uu[tr[3]] - 1) < 0.01
    if(fd == 0) {
      w = x$fimap[fi+1] + 1
      print(fi)
      print(w)
      kk = which(x$xxx[,2] == w)
      oo = x$xxx[kk,1]
      if(oo >= 3000) {
        oo <- x$fmap2[which(x$fmap2[,1] == oo), 2]
      }
      col2 <- ifelse(oo==0, "purple", cols2[oo+1])
      #oo = which(x$fimap == fd)[kk]
      #kk <- kk+1
    }
    #fd = x$zz[fi]
    colors[fi] = ifelse(fd %in% x$xxx[,2], rcols2[1+x$xxx[which(x$xxx[,2]==fd), 1]], ifelse(fd < 6000, cols1[fd+1], cols1[1+x$fmap1[which(x$fmap1[,1] == fd), 2]]))
  } else {
    colors[fi] = cols1[fi]
  }
}
rmesh1$material$color <- colors