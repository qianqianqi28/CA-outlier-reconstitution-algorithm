rm(list=ls()) 
library(anacor)
library(readxl)
library(xtable)
library(FactoMineR)
library(cellWise)
library(Matrix)
source("R/reconca.R")

set.seed(1234)

###row

# Step 1: Create random weights and normalize
row_weights <- runif(5)
# row_weights <- c(0.2, 0.2, 0.2, 0.2, 0.2)
row_weights <- row_weights / sum(row_weights)
round(row_weights, 5)

###column

# Step 1: Create random weights and normalize
col_weights <- runif(6)
#col_weights <- c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6)
col_weights <- col_weights / sum(col_weights)
round(col_weights, 5)

# Step 2: Create Dr and Dc (diagonal matrices)
Dr <- diag(row_weights)   # square roots used because Dr appears on both sides
Dc <- diag(col_weights)
# Step 5: All ones matrix
ones_mat <- matrix(1, nrow = 5, ncol = 6)


###row

# Step 2: Create first vector x1 (mean zero, norm 1)
x1_raw <- c(1, 2, 3, 4, 5)
x1_centered <- x1_raw - sum(row_weights * x1_raw)
x1 <- x1_centered / sqrt(sum(row_weights * x1_centered^2))

###column

# Step 2: Create first vector x1 (mean zero, norm 1)
y1_raw <- c(1, 2, 3, 4, 5, 6)
y1_centered <- y1_raw - sum(col_weights * y1_raw)
y1 <- y1_centered / sqrt(sum(col_weights * y1_centered^2))

###row

# Step 3: Create second vector y_raw
x2_raw <- x1_raw^2

# Step 4: Remove weighted mean from y_raw
x2_centered <- x2_raw - sum(row_weights * x2_raw)

# Step 5: Remove projection of y_centered onto x1
x1_x2_dot_product <- sum(row_weights * x1 * x2_centered)
x2_proj <- x2_centered - x1_x2_dot_product * x1

# Step 6: Normalize y to have weighted norm = 1
x2 <- x2_proj / sqrt(sum(row_weights * x2_proj^2))

###column

# Step 3: Create second vector y_raw
y2_raw <- y1_raw^2

# Step 4: Remove weighted mean from y_raw
y2_centered <- y2_raw - sum(col_weights * y2_raw)

# Step 5: Remove projection of y_centered onto y1
y1_y2_dot_product <- sum(col_weights * y1 * y2_centered)
y2_proj <- y2_centered - y1_y2_dot_product * y1

# Step 6: Normalize y to have weighted norm = 1
y2 <- y2_proj / sqrt(sum(col_weights * y2_proj^2))

# Step 7: Verify constraints
list(
  weights = round(row_weights, 3),
  x2 = round(x2, 3),
  weighted_dot = round(sum(row_weights * x1 * x2), 6),
  x2_weighted_mean = round(sum(row_weights * x2), 6),
  x2_norm = round(sum(row_weights * x2^2), 6),
  
  weights = round(col_weights, 3),
  y2 = round(y2, 3),
  weighted_dot = round(sum(col_weights * y1 * y2), 6),
  y2_weighted_mean = round(sum(col_weights * y2), 6),
  y2_norm = round(sum(col_weights * y2^2), 6)
)
# Step 4: Choose singular value
sigma1 <- 0.42
sigma2 <- 0.18
# Step 6: Reconstruct M
M3 <- Dr %*% (ones_mat + sigma1 * outer(x1, y1) + sigma2 * outer(x2, y2)) %*% Dc

# Optional: Round for display
round(M3, 3)

sum(M3)

rankMatrix(M3)

rownames(M3) <- c(1, 2, 3, 4, 5)
colnames(M3) <- c("a", "b", "c", "d", "e", "f")
print(xtable(M3,digits=rep(5, (ncol(M3)+1))), include.rownames=T, include.colnames=T)


M3 <- as.data.frame(M3)
colnames(M3)
rownames(M3)

dt <- M3

#dt[1,1] = 2*max(dt)
rankMatrix(dt)
rownames(dt);colnames(dt)

dt[3,1] <- 4*max(M3)
round(dt[3,1], 5)
print(xtable(dt,digits=rep(5, (ncol(dt)+1))), include.rownames=T, include.colnames=T)

#dt[1,1] = 2*max(dt)
rankMatrix(dt)
rownames(dt);colnames(dt)
X <- dt
#start CA
X.P    <- X/sum(X)
X.r    <- apply(X.P, 1, sum)
X.c    <- apply(X.P, 2, sum)

X.Dr   <- diag(X.r)
X.Dc   <- diag(X.c)
X.Drmh <- diag(1/sqrt(X.r))
X.Dcmh <- diag(1/sqrt(X.c))

X.P   <- as.matrix(X.P)
X.S   <- X.Drmh%*%(X.P-X.r%o%X.c)%*%X.Dcmh
X.svd <- svd(X.S)

#obtain singular values
round((X.svd$d), 5)


colproj    <- X.Dcmh %*% X.svd$v %*% diag(X.svd$d)
rownames(colproj) <- colnames(X)

rowproj    <- X.Drmh %*% X.svd$u %*% diag(X.svd$d)
rownames(rowproj) <- rownames(X)

#setting parameter for plot


#setting parameter for plot
Flip <- "flip y"

xaxis <- 1
yaxis <- 2

if(Flip == "flip x"){
  rowproj[,xaxis] <- -rowproj[,xaxis]
  colproj[,xaxis] <- -colproj[,xaxis]
} else if (Flip == "flip y") {
  rowproj[,yaxis] <- -rowproj[,yaxis]
  colproj[,yaxis] <- -colproj[,yaxis]
} else if (Flip == "flip x and y") {
  rowproj[,c(xaxis, yaxis)] <- -rowproj[,c(xaxis, yaxis)]
  colproj[,c(xaxis, yaxis)] <- -colproj[,c(xaxis, yaxis)]
} else{
  rowproj <- rowproj
  colproj <- colproj
}

axisminx <- min(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))
axismaxx <- max(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))

axisminy <- min(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))
axismaxy <- max(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))


axisminx <- -1
axismaxx <- 1.2

axisminy <- -1
axismaxy <- 1.2

savefigure <- "simulationrank3out31casym"
pdf(file = paste0("add plots/f2/",savefigure,".pdf"), width = 12,
    height = 12)
par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right
par(mgp = c(4, 1.8, 0))
plot(rowproj[,c(xaxis)],rowproj[,c(yaxis)], cex.axis = 2.5, pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
     ylim=c(axisminy,axismaxy),# xaxs="i", yaxs="i",
     xlab=paste0("Dim", xaxis, ": ",format(round((X.svd$d[xaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[xaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), ylab=paste0("Dim", yaxis, ": ", format(round((X.svd$d[yaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[yaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), fg = "gray", cex.lab=2.5, cex = 4.5)

points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 4.5)

abline(v=0, h=0, col="gray")
tempforlab <- rbind(rowproj[,c(xaxis, yaxis)], colproj[,c(xaxis, yaxis)])
autoLab(tempforlab[,c(xaxis)], tempforlab[,c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(rowproj)), rep("red", nrow(colproj))), cex.lab=2.5, cex = 4.5)

dev.off() 

# Perform Chi-square test of independence
chi_square_test <- chisq.test(dt)

cat("\nChi-square Test of Independence:\n")
print(chi_square_test)

round((X.svd$d), 5)
round((X.svd$d^2), 5)
round(100*(X.svd$d^2)/sum(X.svd$d^2), 1)

# check the property for outlying columns and rows

outlycol <- c("a")
outlyrow <- c("3")
outliercell <- c("3 a")

for (i in c(1:length(outlyrow))){
  print(apply(X.P, 1, sum)[which(rownames(dt) == outlyrow[i])])
  print(X.svd$u[which(rownames(dt) == outlyrow[i]),1]^2)
}
apply(X.P, 1, sum)
for (i in c(1:length(outlycol))){
  print(apply(X.P, 2, sum)[which(colnames(dt) == outlycol[i])])
  print(X.svd$v[which(colnames(dt) == outlycol[i]),1]^2)
}

apply(X.P, 2, sum)

vector_2x3 <- do.call(paste, expand.grid(rownames(X), colnames(X), stringsAsFactors = FALSE))
matrix_2x3 <- matrix(vector_2x3, nrow = nrow(X))
for (i in c(1:length(outliercell))){
  print(which(matrix_2x3 == outliercell[i], arr.ind = TRUE))
  print((X.S^2)[which(matrix_2x3 == outliercell[i], arr.ind = TRUE)]/sum(X.S^2))
}

sum(X.svd$v[,2]^2)
sum(X.svd$u[,2]^2)


sum(M3)

rankMatrix(M3)
rankMatrix(dt)

dt <- as.matrix(dt)
#Reconstitution of order 0

outliercell <- c("3 a")

recval <- reconca(dt, ncp = 0, initials = 1, threshold = 1e-15, outliercell = outliercell, maxiter = 100000, verbose=TRUE)

X <- recval[[1]]
class(X)

#get the value for outlying cells; the value is -0.0006491913; therefore, we move to reconstitution of order 0
X[which(rownames(dt) == "3"), which(colnames(dt) == "a")]

rankMatrix(X)

#start CA
X.P    <- X/sum(X)
X.r    <- apply(X.P, 1, sum)
X.c    <- apply(X.P, 2, sum)

X.Dr   <- diag(X.r)
X.Dc   <- diag(X.c)
X.Drmh <- diag(1/sqrt(X.r))
X.Dcmh <- diag(1/sqrt(X.c))

X.P   <- as.matrix(X.P)
X.S   <- X.Drmh%*%(X.P-X.r%o%X.c)%*%X.Dcmh
X.svd <- svd(X.S)


colproj    <- X.Dcmh %*% X.svd$v %*% diag(X.svd$d)
rownames(colproj) <- colnames(X)

rowproj    <- X.Drmh %*% X.svd$u %*% diag(X.svd$d)
rownames(rowproj) <- rownames(X)

xaxis <- 1
yaxis <- 2

axisminx <- min(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))
axismaxx <- max(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))

axisminy <- min(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))
axismaxy <- max(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))


axisminx <- -1
axismaxx <- 1.2

axisminy <- -1
axismaxy <- 1.2

savefigure <- "simulationrank3out31rec0casym"
pdf(file = paste0("add plots/f2/",savefigure,".pdf"), width = 12,
    height = 12)
par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right
par(mgp = c(4, 1.8, 0))

plot(rowproj[,c(xaxis)],rowproj[,c(yaxis)], cex.axis = 2.5, pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
     ylim=c(axisminy,axismaxy),# xaxs="i", yaxs="i",
     xlab=paste0("Dim", xaxis, ": ",format(round((X.svd$d[xaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[xaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), ylab=paste0("Dim", yaxis, ": ", format(round((X.svd$d[yaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[yaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), fg = "gray", cex.lab=2.5, cex = 4.5)

points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 4.5)

abline(v=0, h=0, col="gray")
tempforlab <- rbind(rowproj[,c(xaxis, yaxis)], colproj[,c(xaxis, yaxis)])
autoLab(tempforlab[,c(xaxis)], tempforlab[,c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(rowproj)), rep("red", nrow(colproj))), cex.lab=2.5, cex = 4.5)

dev.off() 

# Perform Chi-square test of independence
chi_square_test <- chisq.test(dt)

cat("\nChi-square Test of Independence:\n")
print(chi_square_test)


round((X.svd$d), 5)
round((X.svd$d^2), 5)
round(100*(X.svd$d^2)/sum(X.svd$d^2), 5)

# check the property for outlying columns and rows

outlycol <- c("a")
outlyrow <- c("3")
outliercell <- c("3 a")

for (i in c(1:length(outlyrow))){
  print(apply(X.P, 1, sum)[which(rownames(dt) == outlyrow[i])])
  print(X.svd$u[which(rownames(dt) == outlyrow[i]),1]^2)
}

for (i in c(1:length(outlycol))){
  print(apply(X.P, 2, sum)[which(colnames(dt) == outlycol[i])])
  print(X.svd$v[which(colnames(dt) == outlycol[i]),1]^2)
}



vector_2x3 <- do.call(paste, expand.grid(rownames(X), colnames(X), stringsAsFactors = FALSE))
matrix_2x3 <- matrix(vector_2x3, nrow = nrow(X))
for (i in c(1:length(outliercell))){
  print(which(matrix_2x3 == outliercell[i], arr.ind = TRUE))
  print((X.S^2)[which(matrix_2x3 == outliercell[i], arr.ind = TRUE)]/sum(X.S^2))
}

sum(X.svd$v[,2]^2)
sum(X.svd$u[,2]^2)


rankMatrix(M3)
rankMatrix(X)


recval <- reconca(dt, ncp = 1, initials = 1, threshold = 1e-15, outliercell = outliercell, maxiter = 100000, verbose=TRUE)

X <- recval[[1]]
class(X)

#get the value for outlying cells; the value is -0.0006491913; therefore, we move to reconstitution of order 0
X[which(rownames(dt) == "3"), which(colnames(dt) == "a")]

rankMatrix(X)

#start CA
X.P    <- X/sum(X)
X.r    <- apply(X.P, 1, sum)
X.c    <- apply(X.P, 2, sum)

X.Dr   <- diag(X.r)
X.Dc   <- diag(X.c)
X.Drmh <- diag(1/sqrt(X.r))
X.Dcmh <- diag(1/sqrt(X.c))

X.P   <- as.matrix(X.P)
X.S   <- X.Drmh%*%(X.P-X.r%o%X.c)%*%X.Dcmh
X.svd <- svd(X.S)

#obtain singular values
round((X.svd$d), 5)


colproj    <- X.Dcmh %*% X.svd$v %*% diag(X.svd$d)
rownames(colproj) <- colnames(X)

rowproj    <- X.Drmh %*% X.svd$u %*% diag(X.svd$d)
rownames(rowproj) <- rownames(X)

xaxis <- 1
yaxis <- 2

axisminx <- min(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))
axismaxx <- max(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))

axisminy <- min(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))
axismaxy <- max(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))


axisminx <- -1
axismaxx <- 1.2

axisminy <- -1
axismaxy <- 1.2

savefigure <- "simulationrank3out31rec1casym"
pdf(file = paste0("add plots/f2/",savefigure,".pdf"), width = 12,
    height = 12)
par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right
par(mgp = c(4, 1.8, 0))

plot(rowproj[,c(xaxis)],rowproj[,c(yaxis)], cex.axis = 2.5, pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
     ylim=c(axisminy,axismaxy),# xaxs="i", yaxs="i",
     xlab=paste0("Dim", xaxis, ": ",format(round((X.svd$d[xaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[xaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), ylab=paste0("Dim", yaxis, ": ", format(round((X.svd$d[yaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[yaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), fg = "gray", cex.lab=2.5, cex = 4.5)

points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 4.5)

abline(v=0, h=0, col="gray")
tempforlab <- rbind(rowproj[,c(xaxis, yaxis)], colproj[,c(xaxis, yaxis)])
autoLab(tempforlab[,c(xaxis)], tempforlab[,c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(rowproj)), rep("red", nrow(colproj))), cex.lab=2.5, cex = 4.5)

dev.off() 

# Perform Chi-square test of independence
chi_square_test <- chisq.test(dt)

cat("\nChi-square Test of Independence:\n")
print(chi_square_test)


round((X.svd$d), 5)
round((X.svd$d^2), 5)
round(100*(X.svd$d^2)/sum(X.svd$d^2), 5)

# check the property for outlying columns and rows

outlycol <- c("a")
outlyrow <- c("3")
outliercell <- c("3 a")

for (i in c(1:length(outlyrow))){
  print(apply(X.P, 1, sum)[which(rownames(dt) == outlyrow[i])])
  print(X.svd$u[which(rownames(dt) == outlyrow[i]),1]^2)
}

for (i in c(1:length(outlycol))){
  print(apply(X.P, 2, sum)[which(colnames(dt) == outlycol[i])])
  print(X.svd$v[which(colnames(dt) == outlycol[i]),1]^2)
}



vector_2x3 <- do.call(paste, expand.grid(rownames(X), colnames(X), stringsAsFactors = FALSE))
matrix_2x3 <- matrix(vector_2x3, nrow = nrow(X))
for (i in c(1:length(outliercell))){
  print(which(matrix_2x3 == outliercell[i], arr.ind = TRUE))
  print((X.S^2)[which(matrix_2x3 == outliercell[i], arr.ind = TRUE)]/sum(X.S^2))
}

sum(X.svd$v[,2]^2)
sum(X.svd$u[,2]^2)


rankMatrix(M3)
rankMatrix(X)


recval <- reconca(dt, ncp = 2, initials = 1, threshold = 1e-15, outliercell = outliercell, maxiter = 100000, verbose=TRUE)

X <- recval[[1]]
class(X)

#get the value for outlying cells; the value is -0.0006491913; therefore, we move to reconstitution of order 0
X[which(rownames(dt) == "3"), which(colnames(dt) == "a")]
M3[1,1]
rankMatrix(X)



recval <- reconca(dt, ncp = 3, initials = 2, threshold = 1e-15, outliercell = outliercell, maxiter = 100000, verbose=TRUE)

X <- recval[[1]]
class(X)

#get the value for outlying cells; the value is -0.0006491913; therefore, we move to reconstitution of order 0
X[which(rownames(dt) == "3"), which(colnames(dt) == "a")]
M3[1,1]

rankMatrix(X)

