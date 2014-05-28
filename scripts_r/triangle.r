# generate points on the triangle (junction configurations)

rm(list = ls()); # clear 

pmax <- 0.6;
pmin <- (1-pmax)/2;
N <- 3;

p1 <- c(pmax, pmin, pmin);
p2 <- c(pmin, pmax, pmin);
p3 <- c(pmin, pmin, pmax);
p0 <- c(1/3, 1/3, 1/3);
v1 <- p1 - p3;
v2 <- p2 - p1;
a  <- sqrt(sum(v1^2));
v1 <- v1 / a;
v2 <- v2 / a;
L <- a / 2^N;

cnt <- 0;
for (i in 0:2^N) {
  for (j in 0:i) cnt <- cnt + 1;
}

pts <- matrix(0, nrow=cnt, ncol=3);
cnt <- 0;
for (i in 0:2^N) {
  for (j in 0:i) {
    cnt <- cnt + 1;
    pp <- p3 + i * L * v1 + j * L * v2;
    print(pp);
    pts[cnt, 1:3] <- pp;
  }
}

library(plot3D)
scatter3D(pts[,1], pts[,2], pts[,3],
          colvar = NULL,
          type="p",
          phi = 10, theta = 80,
          col = NULL, NAcol = "white", 
          main="triangle points", 
          xlim=c(0,1), ylim=c(0,1), zlim=c(0,1),  
          pch=20
          );

print(paste(nrow(pts), " points", sep=""));

readline(prompt = "press enter to continue...")

# plot distances towards symmetry
d <- vector(mode="numeric", length=nrow(pts));
for (i in 1:nrow(pts)) {
  d[i] <- sqrt(sum((pts[i,]-p0)^2));
}

d <- sort(d);
plot(d, pch=20);

readline(prompt = "press enter to continue...")

clust_d <- cutree(hclust(dist(d)), h=0.06); ## possible to cluster it
plot(clust_d, pch=20)

readline(prompt = "press enter to continue...")

# plot ratios max/min
ratios <- vector(mode="numeric", length=nrow(pts));
for (i in 1:nrow(pts)) {
  ratios[i] <- max(pts[i,])/min(pts[i,]); #sqrt(sum((pts[i,]-p0)^2));
}

plot(sort(ratios), pch=20);

readline(prompt = "press enter to continue...")

vy <- p1 - (p2+p3)/2;
vy <- vy / sqrt(sum(vy^2));
vx <- p3 - (p2+p3)/2;
vx <- vx / sqrt(sum(vx^2));
x <- (pts-p0) %*% vx;
y <- (pts-p0) %*% vy;

plot(x,y, pch=20, asp=1)

readline(prompt = "press enter to continue...")
library(plot3D)
scatter2D(x, y, colvar=ratios, pch=20, xlab="", ylab="");