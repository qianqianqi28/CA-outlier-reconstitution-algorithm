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

dt[1,1] <- 4*max(M3)

dt <- as.matrix(dt)

X <- dt
dim(X)
dim(X)[1]*dim(X)[2]
sum(X)
sum(X == 0)


##########################
##  check multivariate normal distribution
##########################

X.P    <- X/sum(X)
X.r    <- apply(X.P, 1, sum)
X.c    <- apply(X.P, 2, sum)
X.Dr   <- diag(X.r)
X.Dc   <- diag(X.c)
X.Drmh <- diag(1/sqrt(X.r))
X.Dcmh <- diag(1/sqrt(X.c))

X.P   <- as.matrix(X.P)
X.S   <- X.Drmh%*%(X.P-X.r%o%X.c)%*%X.Dcmh
S <- X.S
rownames(S) <- rownames(X)
colnames(S) <- colnames(X)

set.seed(123)

ggdensity(c(S), 
          main = "Density plot of S",
          xlab = "S")

ggqqplot(c(S))
shapiro.test(c(S))

##########################
##  MacroPCA
##########################
set.seed(123)

d <- ncol(X)

MacroPCApar0 <- list(scale = FALSE, center = rep(0, d), alpha = 0.80)
class(S)
rownames(S)
colnames(S)
MacroPCA.out <- MacroPCA(S, k = 1, MacroPCApars = MacroPCApar0)

# DDC

# outlying rows in the DDC
indrows <- MacroPCA.out$DDC$indrows
rownames(X)[indrows]

# outlying cells in the DDC
length(MacroPCA.out$DDC$indcells)

savefigureddc <- "simulationrank3out11mac1outDDCcelloutlier"

pdf(paste0("add plots/f2/", savefigureddc, ".pdf"), width = 12,
    height = 12)

par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right

ggpmacropca     <- cellMap(R = MacroPCA.out$DDC$stdResid,                        
                           indcells = MacroPCA.out$DDC$indcells,
                           indrows = MacroPCA.out$DDC$indrows,
                           mTitle = "",
                           rowtitle = "",
                           columntitle = "",
                           sizetitles = 2.4)

My_Theme = theme(
  axis.title.x = element_text(size = 40),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 40),
  axis.title.y = element_text(size = 30))
ggpmacropca+My_Theme
dev.off()

# CA plot
V <- MacroPCA.out$loadings
scores <- MacroPCA.out$scores
sVals      <- sqrt(nrow(S)*MacroPCA.out$eigenvalues)
U          <- 1/sVals * scores
rowproj    <-  sVals*(X.Drmh %*% U)
rowproj <- cbind(rowproj, c(0,0,0,0,0))

rownames(rowproj) <- rownames(S)


colproj    <- sVals*(X.Dcmh %*% V)
colproj <- cbind(colproj, c(0,0,0,0,0,0))


rownames(colproj) <- colnames(S)

sort(rowproj[,1])
sort(colproj[,1])

Flip <- "flip x"

xaxis <- 1
yaxis <- 2

savefigure <- "mac1simulationrank3out11casym"


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

pdf(paste0("add plots/f2/", savefigure, ".pdf"), width = 12,
    height = 12)
par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right
par(mgp = c(4, 1.8, 0))
plot(rowproj[,c(xaxis)], rowproj[,c(yaxis)], cex.axis = 2.5, pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
     ylim=c(axisminy,axismaxy), xlab=paste0("Dim", xaxis), ylab=paste0("Dim", yaxis), fg = "gray", cex.lab=2.5, cex = 4.5)

points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 4.5)

abline(v=0, h=0, col="gray")

tempforlab <- rbind(rowproj, colproj)

autoLab(tempforlab[, c(xaxis)], tempforlab[, c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(X)), rep("red", ncol(X))), cex.lab=2.5, cex = 4.5)

dev.off()





##########################
##  MacroPCA
##########################
set.seed(123)

d <- ncol(X)

MacroPCApar0 <- list(scale = FALSE, center = rep(0, d), alpha = 0.80)
class(S)
rownames(S)
colnames(S)
MacroPCA.out <- MacroPCA(S, k = 2, MacroPCApars = MacroPCApar0)

# DDC

# outlying rows in the DDC
indrows <- MacroPCA.out$DDC$indrows
rownames(X)[indrows]

# outlying cells in the DDC
length(MacroPCA.out$DDC$indcells)

savefigureddc <- "simulationrank3out11mac2outDDCcelloutlier"

pdf(paste0("add plots/f2/", savefigureddc, ".pdf"), width = 12,
    height = 12)
par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right

ggpmacropca     <- cellMap(R = MacroPCA.out$DDC$stdResid,                        
                           indcells = MacroPCA.out$DDC$indcells,
                           indrows = MacroPCA.out$DDC$indrows,
                           mTitle = "",
                           rowtitle = "",
                           columntitle = "",
                           sizetitles = 2.4)

My_Theme = theme(
  axis.title.x = element_text(size = 40),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 40),
  axis.title.y = element_text(size = 30))
ggpmacropca+My_Theme
dev.off()

# CA plot
V <- MacroPCA.out$loadings
scores <- MacroPCA.out$scores
sVals      <- sqrt(nrow(S)*MacroPCA.out$eigenvalues)
U          <- scores %*% diag(1/sVals)
rowproj    <- X.Drmh %*% U %*% diag(sVals)

rownames(rowproj) <- rownames(S)


colproj    <- X.Dcmh %*% V %*% diag(sVals)
rownames(colproj) <- colnames(S)

sort(rowproj[,1])
sort(colproj[,1])

Flip <- "flip y"

xaxis <- 1
yaxis <- 2

savefigure <- "mac2simulationrank3out11casym"


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

pdf(paste0("add plots/f2/", savefigure, ".pdf"), width = 12,
    height = 12)
par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right
par(mgp = c(4, 1.8, 0))

plot(rowproj[,c(xaxis)], rowproj[,c(yaxis)], cex.axis = 2.5, pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
     ylim=c(axisminy,axismaxy), xlab=paste0("Dim", xaxis), ylab=paste0("Dim", yaxis), fg = "gray", cex.lab=2.5, cex = 4.5)

points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 4.5)

abline(v=0, h=0, col="gray")

tempforlab <- rbind(rowproj, colproj)

autoLab(tempforlab[, c(xaxis)], tempforlab[, c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(X)), rep("red", ncol(X))), cex.lab=2.5, cex = 4.5)

dev.off()


# ##########################
# ##  MacroPCA
# ##########################
# set.seed(123)
# 
# d <- ncol(X)
# 
# MacroPCApar0 <- list(scale = FALSE, center = rep(0, d), alpha = 0.80)
# class(S)
# rownames(S)
# colnames(S)
# MacroPCA.out <- MacroPCA(S, k = 3, MacroPCApars = MacroPCApar0)
# 
# # DDC
# 
# # outlying rows in the DDC
# indrows <- MacroPCA.out$DDC$indrows
# rownames(X)[indrows]
# 
# # outlying cells in the DDC
# length(MacroPCA.out$DDC$indcells)
# 
# savefigureddc <- "simulationrank3out11mac3outDDCcelloutlier"
# 
# pdf(paste0("add plots/f2/", savefigureddc, ".pdf"), width = 12,
#     height = 12)
# ggpmacropca     <- cellMap(R = MacroPCA.out$DDC$stdResid,                        
#                            indcells = MacroPCA.out$DDC$indcells,
#                            indrows = MacroPCA.out$DDC$indrows,
#                            mTitle = "",
#                            rowtitle = "",
#                            columntitle = "",
#                            sizetitles = 2.4)
# 
# My_Theme = theme(
#   axis.title.x = element_text(size = 40),
#   axis.text.x = element_text(size = 30),
#   axis.text.y = element_text(size = 40),
#   axis.title.y = element_text(size = 30))
# ggpmacropca+My_Theme
# dev.off()
# 
# # CA plot
# V <- MacroPCA.out$loadings
# scores <- MacroPCA.out$scores
# sVals      <- sqrt(nrow(S)*MacroPCA.out$eigenvalues)
# U          <- scores %*% diag(1/sVals)
# rowproj    <- X.Drmh %*% U %*% diag(sVals)
# 
# rownames(rowproj) <- rownames(S)
# 
# 
# colproj    <- X.Dcmh %*% V %*% diag(sVals)
# rownames(colproj) <- colnames(S)
# 
# sort(rowproj[,1])
# sort(colproj[,1])
# 
# Flip <- "no flip"
# Checkoverlap <- FALSE
# 
# xaxis <- 1
# yaxis <- 2
# 
# savefigure <- "mac3simulationrank3out11casym"
# 
# 
# if(Flip == "flip x"){
#   rowproj[,xaxis] <- -rowproj[,xaxis]
#   colproj[,xaxis] <- -colproj[,xaxis]
# } else if (Flip == "flip y") {
#   rowproj[,yaxis] <- -rowproj[,yaxis]
#   colproj[,yaxis] <- -colproj[,yaxis]
# } else if (Flip == "flip x and y") {
#   rowproj[,c(xaxis, yaxis)] <- -rowproj[,c(xaxis, yaxis)]
#   colproj[,c(xaxis, yaxis)] <- -colproj[,c(xaxis, yaxis)]
# } else{
#   rowproj <- rowproj
#   colproj <- colproj
# }
# 
# axisminx <- min(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))
# axismaxx <- max(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))
# 
# axisminy <- min(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))
# axismaxy <- max(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))
# 
# axisminx <- -1
# axismaxx <- 1.2
# 
# axisminy <- -1
# axismaxy <- 1.2
# 
# pdf(paste0("add plots/f2/", savefigure, ".pdf"), width = 12,
#     height = 12)
# 
# plot(rowproj[,c(xaxis)], rowproj[,c(yaxis)], pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
#      ylim=c(axisminy,axismaxy), xlab=paste0("Dim", xaxis), ylab=paste0("Dim", yaxis), fg = "gray", cex.lab=2.5, cex = 4.5)
# 
# points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 4.5)
# 
# abline(v=0, h=0, col="gray")
# 
# tempforlab <- rbind(rowproj, colproj)
# 
# autoLab(tempforlab[, c(xaxis)], tempforlab[, c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(X)), rep("red", ncol(X))), cex.lab=2.5, cex = 4.5)
# 
# dev.off()