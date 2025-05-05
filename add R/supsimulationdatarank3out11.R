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

#Delete outlying rows and columns
roww <- c("1")
col <- c("a")
X <- X[-c(which(rownames(X) %in% roww)),-c(which(colnames(X) %in% col))]
X <- as.matrix(X)
dim(X)
rankMatrix(X)
##########################
##  CA
##########################

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

round((X.svd$d), 3)
round((X.svd$d^2), 3)
round(100*(X.svd$d^2)/sum(X.svd$d^2), 1)

colproj    <- X.Dcmh %*% X.svd$v %*% diag(X.svd$d)
rownames(colproj) <- colnames(X)

rowproj    <- X.Drmh %*% X.svd$u %*% diag(X.svd$d)
rownames(rowproj) <- rownames(X)

#setting parameter for plot
Flip <- "flip y"

xaxis <- 1
yaxis <- 2

savefigure <- "supsimulationrank3out11casym"

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

sort(colproj[,1])
sort(rowproj[,1])


sort(colproj[,2])
sort(rowproj[,2])


##########################
## obtain coordinates for supplementary points
##########################
row_sup <- matrix(data=NA,nrow=length(roww),ncol=ncol(X.svd$v)) 
col_sup <- matrix(data=NA,nrow=length(col),ncol=ncol(X.svd$u)) 

if (class(dt)[1] == c("matrix")){
  print("matrix")
  for (i in 1:length(roww)){
    if (sum(dt[c(which(rownames(dt) %in% roww[i])),-c(which(colnames(dt) %in% col))]) == 0){
      row_sup[i, ] <- 0
    } else{
      row_sup[i, ] <- t(as.matrix(dt[c(which(rownames(dt) %in% roww[i])),-c(which(colnames(dt) %in% col))]/sum(dt[c(which(rownames(dt) %in% roww[i])),-c(which(colnames(dt) %in% col))]))) %*% X.Dcmh %*% X.svd$v
    }
  }
  
  
  for (j in 1:length(col)){
    if (sum(dt[-c(which(rownames(dt) %in% roww)),c(which(colnames(dt) %in% col[j]))]) == 0){
      col_sup[j, ] <- 0
    } else{
      col_sup[j, ] <- t(as.matrix(dt[-c(which(rownames(dt) %in% roww)),c(which(colnames(dt) %in% col[j]))]/sum(dt[-c(which(rownames(dt) %in% roww)),c(which(colnames(dt) %in% col[j]))]))) %*% X.Drmh %*% X.svd$u
    }
  }
  
} else if (class(dt) == "data.frame"){
  print("data.frame'")
  for (i in 1:length(roww)){
    if (sum(dt[c(which(rownames(dt) %in% roww[i])),-c(which(colnames(dt) %in% col))]) == 0){
      row_sup[i, ] <- 0
    } else{
      row_sup[i, ] <- as.matrix(dt[c(which(rownames(dt) %in% roww[i])),-c(which(colnames(dt) %in% col))]/sum(dt[c(which(rownames(dt) %in% roww[i])),-c(which(colnames(dt) %in% col))])) %*% X.Dcmh %*% X.svd$v
    }
  }
  
  for (j in 1:length(col)){
    if (sum(dt[-c(which(rownames(dt) %in% roww)),c(which(colnames(dt) %in% col[j]))]) == 0){
      col_sup[j, ] <- 0
    } else{
      col_sup[j, ] <- t(as.matrix(dt[-c(which(rownames(dt) %in% roww)),c(which(colnames(dt) %in% col[j]))]/sum(dt[-c(which(rownames(dt) %in% roww)),c(which(colnames(dt) %in% col[j]))]))) %*% X.Drmh %*% X.svd$u
    }
  }
  
}


# row_sup[1,2] <- 0
# col_sup[1,2] <- 0

if(Flip == "flip x"){
  row_sup[,xaxis] <- -row_sup[,xaxis]
  col_sup[,xaxis] <- -col_sup[,xaxis]
} else if (Flip == "flip y") {
  row_sup[,yaxis] <- -row_sup[,yaxis]
  col_sup[,yaxis] <- -col_sup[,yaxis]
} else if (Flip == "flip x and y") {
  row_sup[,c(xaxis, yaxis)] <- -row_sup[,c(xaxis, yaxis)]
  col_sup[,c(xaxis, yaxis)] <- -col_sup[,c(xaxis, yaxis)]
} else{
  row_sup <- row_sup
  col_sup <- col_sup
}

axisminx <- min(c(rowproj[,c(xaxis)],row_sup[,c(xaxis)],col_sup[,c(xaxis)], colproj[,c(xaxis)]))
axismaxx <- max(c(rowproj[,c(xaxis)],row_sup[,c(xaxis)],col_sup[,c(xaxis)], colproj[,c(xaxis)]))

axisminy <- min(c(rowproj[,c(yaxis)],row_sup[,c(yaxis)],col_sup[,c(yaxis)], colproj[,c(yaxis)]))
axismaxy <- max(c(rowproj[,c(yaxis)],row_sup[,c(yaxis)],col_sup[,c(yaxis)], colproj[,c(yaxis)]))

axisminx <- -1
axismaxx <- 1.2

axisminy <- -1
axismaxy <- 1.2

pdf(file = paste0("add plots/f2/",savefigure,".pdf"), width = 12,
    height = 12)

par(mar = c(5, 6, 4, 2) + 0.1)  # Bottom, left, top, right
par(mgp = c(4, 1.8, 0))

plot(rowproj[,c(xaxis)], rowproj[,c(yaxis)], cex.axis = 2.5, pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
     ylim=c(axisminy,axismaxy),
     xlab=paste0("Dim", xaxis, ": ",format(round((X.svd$d[xaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[xaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), ylab=paste0("Dim", yaxis, ": ", format(round((X.svd$d[yaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[yaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), fg = "gray", cex.lab=2.5, cex = 4.5)

points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 4.5)
points(row_sup[,c(xaxis)],row_sup[,c(yaxis)], col = "blue", pch = 23, cex = 4.5)
points(col_sup[,c(xaxis)],col_sup[,c(yaxis)], col = "red", pch = 23, cex = 4.5)
abline(v=0, h=0, col="gray")
tempforlab <- rbind(rowproj[,c(xaxis, yaxis)], colproj[,c(xaxis, yaxis)])
autoLab(tempforlab[, c(xaxis)], tempforlab[, c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(rowproj)), rep("red", nrow(colproj))), cex.lab=2.5, cex = 4.5)
text(row_sup[,c(xaxis)],row_sup[,c(yaxis)]+0.01*(axismaxy-axisminy), labels = c(roww), col = c(rep("blue", length(roww))), cex.lab=2.5, cex = 4.5)
text(col_sup[,c(xaxis)],col_sup[,c(yaxis)]+0.01*(axismaxy-axisminy), labels = c(col), col = c(rep("red", length(col))), cex.lab=2.5, cex = 4.5)

dev.off() 

# Perform Chi-square test of independence
chi_square_test <- chisq.test(X)

cat("\nChi-square Test of Independence:\n")
print(chi_square_test)

# check the property for outlying columns and rows

outlycol <- c("a")
outlyrow <- c("1")
outliercell <- c("1 a")

for (i in c(1:length(outlycol))){
  print(apply(X.P, 2, sum)[which(colnames(dt) == outlycol[i])])
  print(X.svd$v[which(colnames(dt) == outlycol[i]),1]^2)
}

for (i in c(1:length(outlyrow))){
  print(apply(X.P, 1, sum)[which(rownames(dt) == outlyrow[i])])
  print(X.svd$u[which(rownames(dt) == outlyrow[i]),1]^2)
}

vector_2x3 <- do.call(paste, expand.grid(rownames(X), colnames(X), stringsAsFactors = FALSE))
matrix_2x3 <- matrix(vector_2x3, nrow = nrow(X))
for (i in c(1:length(outliercell))){
  print(which(matrix_2x3 == outliercell[i], arr.ind = TRUE))
  print((X.S^2)[which(matrix_2x3 == outliercell[i], arr.ind = TRUE)]/sum(X.S^2))
}

sum(X.svd$v[,2]^2)
sum(X.svd$u[,2]^2)

round((X.svd$d), 5)
round((X.svd$d^2), 5)
round(100*(X.svd$d^2)/sum(X.svd$d^2), 2)

