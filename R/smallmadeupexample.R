#This code is for madeup data in Section 3 How outliers originate
rm(list=ls()) 
library(anacor)
library(readxl)
library(xtable)

#Please load dataset in data/madeup

#dataset <- "madeupblock"
#dataset <- "madeupapproxiblock"
dataset <- "madeuplargevalue"



if (dataset == "madeupblock"){
  dt <- read_excel("data/madeup/madeupblock.xlsx")
} else if (dataset == "madeupapproxiblock"){
  dt <- read_excel("data/madeup/madeupapproxiblock.xlsx")
} else if (dataset == "madeuplargevalue"){
  dt <- read_excel("data/madeup/madeuplargevalue.xlsx")
}

dt <- as.data.frame(dt)
rownames(dt) <- dt[,1]
dt <- dt[,-1]
dim(dt)

#delete 0 vector
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
round((X.svd$d), 2)


colproj    <- X.Dcmh %*% X.svd$v
rownames(colproj) <- colnames(X)

rowproj    <- X.Drmh %*% X.svd$u
rownames(rowproj) <- rownames(X)

#obtain standard coordinates for rows
round(rowproj, 2)
#obtain standard coordinates for columns
round(colproj, 2)

#reorder matrix based on the first dimension
X[order(rowproj[,1]), order(colproj[,1])]


#convenience for user using latex

#obtain standard coordinates for rows
print(xtable(rowproj,digits=rep(2, (ncol(X) + 1))), include.rownames=T, include.colnames=T)

#obtain standard coordinates for columns
print(xtable(colproj,digits=rep(2, (ncol(X) + 1))), include.rownames=T, include.colnames=T)

#reorder matrix based on the first dimension
print(xtable(X[order(rowproj[,1]), order(colproj[,1])],digits=rep(0, (ncol(X) + 1))), include.rownames=T, include.colnames=T)

