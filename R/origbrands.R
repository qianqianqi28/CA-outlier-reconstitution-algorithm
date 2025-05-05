# This is for CA: Figure 1 (a) in the Section Empirical studies/Results
rm(list=ls()) 
library(readxl)
library(xtable)
library(FactoMineR)
library(cellWise)
# The brands dataset can also be loaded by data("data_brands"), same as load except for the row names and column names.
# data("data_brands")

##########################
##  check data
##########################

#load dataset in data/processed
dt <- read.csv("data/processed/brands.csv", header = TRUE)

dt <- as.data.frame(dt)
rownames(dt) <- dt[,1]
dt <- dt[,-1]
dim(dt)

#delete the 0 vectors
if (sum(apply(dt, 1, sum) == 0) != 0){
  dt <- dt[-c(which(apply(dt, 1, sum) == 0)),]
} else {
  dt <- dt}

if (sum(apply(dt, 2, sum) == 0) != 0){
  dt <- dt[,-c(which(apply(dt, 2, sum) == 0))]
} else {
  dt <- dt}
dt <- as.matrix(dt)

X <- dt
dim(X)
dim(X)[1]*dim(X)[2]
sum(X)
sum(X == 0)

#check whether there are identical columns
v1 <-  do.call(paste, as.data.frame(t(X)))
out <- cbind(seq_len(nrow(t(X))), match(v1, unique(v1)))
colrep <- split(out[,1], out[,2])
for (i in 1: length(colrep)){
  if (length(colrep[[i]]) > 1){
    print(colnames(X)[colrep[[i]]])
  }
}

#check whether there are identical rows
v1 <-  do.call(paste, as.data.frame(X))
out <- cbind(seq_len(nrow(X)), match(v1, unique(v1)))
rowrep <- split(out[,1], out[,2])
for (i in 1: length(rowrep)){
  if (length(rowrep[[i]]) > 1){
    print(rownames(X)[rowrep[[i]]])
  }
}


##########################
##  CA
##########################

#start CA
X.P    <- X/sum(X)
X.r    <- apply(X.P, 1, sum)
X.c    <- apply(X.P, 2, sum)
X.Rr    <- apply(X, 1, sum)
X.Cc    <- apply(X, 2, sum)

#convenience for user using latex for original dataset
print(xtable(cbind(dt, X.Rr, X.r),digits=c(rep(0, (ncol(cbind(dt, X.Rr, X.r)))), 3)), include.rownames=T, include.colnames=T)
print(xtable(t(as.matrix(X.Cc)),digits=c(rep(0, (ncol(t(as.matrix(X.Cc))))), 0)))#, include.rownames=T, include.colnames=T)
print(xtable(t(as.matrix(X.c)),digits=c(rep(3, (ncol(t(as.matrix(X.c))))), 3)))#, include.rownames=T, include.colnames=T)


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


# check the property for outlying columns and rows

outlycol <- c("Safety")
outlyrow <- c("Volvo")
outliercell <- c("Volvo Safety")

for (i in c(1:length(outlycol))){
  print(apply(X.P, 2, sum)[which(colnames(dt) == outlycol[i])])
  print(X.svd$v[which(colnames(dt) == outlycol[i]),2]^2)
}

for (i in c(1:length(outlyrow))){
  print(apply(X.P, 1, sum)[which(rownames(dt) == outlyrow[i])])
  print(X.svd$u[which(rownames(dt) == outlyrow[i]),2]^2)
}

vector_2x3 <- do.call(paste, expand.grid(rownames(X), colnames(X), stringsAsFactors = FALSE))
matrix_2x3 <- matrix(vector_2x3, nrow = nrow(X))
for (i in c(1:length(outliercell))){
  print(which(matrix_2x3 == outliercell[i], arr.ind = TRUE))
  print((X.S^2)[which(matrix_2x3 == outliercell[i], arr.ind = TRUE)]/sum(X.S^2))
}

sum(X.svd$v[,2]^2)
sum(X.svd$u[,2]^2)

#setting parameter for plot
Flip <- "flip y"
Checkoverlap <- TRUE

xaxis <- 1
yaxis <- 2

savefigure <- "brandscasym"

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

axisminx <- min(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))
axismaxx <- max(c(rowproj[,c(xaxis)], colproj[,c(xaxis)]))

axisminy <- min(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))
axismaxy <- max(c(rowproj[,c(yaxis)], colproj[,c(yaxis)]))

pdf(file = paste0("plots/",savefigure,".pdf"), width = 12,
    height = 12)



if (Checkoverlap == TRUE){
  plot(rowproj[!duplicated(X),c(xaxis)],rowproj[!duplicated(X),c(yaxis)], pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
       ylim=c(axisminy,axismaxy),# xaxs="i", yaxs="i",
       xlab=paste0("Dim", xaxis, ": ",format(round((X.svd$d[xaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[xaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), ylab=paste0("Dim", yaxis, ": ", format(round((X.svd$d[yaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[yaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), fg = "gray", cex.lab=1.5, cex = 1.5)
  
  points(colproj[!duplicated(t(X)),c(xaxis)], colproj[!duplicated(t(X)),c(yaxis)], col = "red", pch = 17, cex = 1.5)
  
  abline(v=0, h=0, col="gray")
  tempforlab <- rbind(rowproj[!duplicated(X),c(xaxis, yaxis)], colproj[!duplicated(t(X)),c(xaxis, yaxis)])
  autoLab(tempforlab[,c(xaxis)], tempforlab[,c(yaxis)], labels = c(rownames(rowproj[!duplicated(X), ]),rownames(colproj[!duplicated(t(X)), ])), col = c(rep("blue", nrow(rowproj[!duplicated(X), ])), rep("red", nrow(colproj[!duplicated(t(X)), ]))), cex.lab=1.5, cex = 1.5)
  
} else{
  
  plot(rowproj[,c(xaxis)], rowproj[,c(yaxis)], pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
       ylim=c(axisminy,axismaxy),
       xlab=paste0("Dim", xaxis, ": ",format(round((X.svd$d[xaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[xaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), ylab=paste0("Dim", yaxis, ": ", format(round((X.svd$d[yaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[yaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), fg = "gray", cex.lab=1.5, cex = 1.5)
  
  points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 1.5)
  
  abline(v=0, h=0, col="gray")
  tempforlab <- rbind(rowproj[,c(xaxis, yaxis)], colproj[,c(xaxis, yaxis)])
  autoLab(tempforlab[, c(xaxis)], tempforlab[, c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(rowproj)), rep("red", nrow(colproj))), cex.lab=1.5, cex = 1.5)
}
dev.off() 

