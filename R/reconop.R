# This is for CA: Figure 3(b) in the Section Empirical studies/Results
rm(list=ls()) 
library(readxl)
library(xtable)
library(FactoMineR)

source("R/reconca.R")


##########################
##  check data
##########################

#Please load dataset in data/processed
dt <- read_excel("data/processed/OP for Rno3273.xlsx")
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
##  reconstitution of order 2
##########################

#outlying cells
outliercell <- c("17 Resp.C.I", "59 Resp.C.I")

#reconstitution of order 2
recval <- reconca(dt, ncp = 2, initials = 1, threshold = 1e-15, outliercell = outliercell, maxiter = 100000, verbose=TRUE)

#give the reconstituted matrix to X
X <- recval[[1]]
class(X)

#get the value for outlying cells; the value is -0.0006491913; therefore, we move to reconstitution of order 0
X[which(rownames(dt) == "17"), which(colnames(dt) == "Resp.C.I")]
X[which(rownames(dt) == "59"), which(colnames(dt) == "Resp.C.I")]

##The comments are to check whether the convergent values equals to the values for the ca decomposed; equation 16 in the paper.
# #start CA
# X.P    <- X/sum(X)
# X.r    <- apply(X.P, 1, sum)
# X.c    <- apply(X.P, 2, sum)
# X.Dr   <- diag(X.r)
# X.Dc   <- diag(X.c)
# X.Drmh <- diag(1/sqrt(X.r))
# X.Dcmh <- diag(1/sqrt(X.c))
# 
# X.P   <- as.matrix(X.P)
# X.S   <- X.Drmh%*%(X.P-X.r%o%X.c)%*%X.Dcmh
# X.svd <- svd(X.S)
# 
# #check whether the convergent values equals to the values for the ca decomposed; equation 16 in the paper.
# ncp <- 2
# recon <- X.svd$u[, 1:ncp] %*% diag(X.svd$d[1:ncp]) %*% (t(X.svd$v[, 1:ncp]))
# recon <- sum(X) * (t(t(recon * sqrt(X.r)) * sqrt(X.c)) + X.r %*% t(X.c))
# recon[which(rownames(dt) == "17"), which(colnames(dt) == "Resp.C.I")]

##########################
##  reconstitution of order 0
##########################

#Reconstitution of order 0
recval <- reconca(dt, ncp = 0, initials = 1, threshold = 1e-15, outliercell = outliercell, maxiter = 100000, verbose=TRUE)

X <- recval[[1]]
class(X)

#get the value for outlying cells; the value is 0.006535948
X[which(rownames(dt) == "17"), which(colnames(dt) == "Resp.C.I")]
X[which(rownames(dt) == "59"), which(colnames(dt) == "Resp.C.I")]

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

# #check whether the convergent values equals to the values for the ca decomposed; equation 15 in the paper.
recon <- sum(X) * (X.r %*% t(X.c))
recon[which(rownames(dt) == "17"), which(colnames(dt) == "Resp.C.I")]
recon[which(rownames(dt) == "59"), which(colnames(dt) == "Resp.C.I")]


round((X.svd$d), 3)
round((X.svd$d^2), 3)
round(100*(X.svd$d^2)/sum(X.svd$d^2), 1)

colproj    <- X.Dcmh %*% X.svd$v %*% diag(X.svd$d)
rownames(colproj) <- colnames(X)

rowproj    <- X.Drmh %*% X.svd$u %*% diag(X.svd$d)
rownames(rowproj) <- rownames(X)

# check the property for outlying columns and rows

outlycol <- c("Resp.C.I")
outlyrow <- c("17", "59")
outliercell <- c("17 Resp.C.I", "59 Resp.C.I")

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

print(apply(X.P, 1, sum)[which(rownames(dt) == outlyrow[1])] + apply(X.P, 1, sum)[which(rownames(dt) == outlyrow[2])])
print(X.svd$u[which(rownames(dt) == outlyrow[1]),2]^2 + X.svd$u[which(rownames(dt) == outlyrow[2]),2]^2)
print((X.S^2)[which(matrix_2x3 == outliercell[1], arr.ind = TRUE)]/sum(X.S^2) + (X.S^2)[which(matrix_2x3 == outliercell[2], arr.ind = TRUE)]/sum(X.S^2))


#setting parameter for plot
Flip <- "no flip"
Checkoverlap <- TRUE

xaxis <- 1
yaxis <- 2

savefigure <- "recopcasym"

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


