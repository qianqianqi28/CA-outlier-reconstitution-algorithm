# This is for CA: Figure 1 (d) in the Section Empirical studies/Results
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

#Delete outlying rows and columns
roww <- c("Volvo")
col <- c("Safety")
X <- X[-c(which(rownames(X) %in% roww)),-c(which(colnames(X) %in% col))]
X <- as.matrix(X)

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
Flip <- "no flip"
Checkoverlap <- TRUE

xaxis <- 1
yaxis <- 2

savefigure <- "supbrandscasym"

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


pdf(file = paste0("plots/",savefigure,".pdf"), width = 12,
    height = 12)



if (Checkoverlap == TRUE){
  plot(rowproj[!duplicated(X),c(xaxis)],rowproj[!duplicated(X),c(yaxis)], pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
       ylim=c(axisminy,axismaxy),# xaxs="i", yaxs="i",
       xlab=paste0("Dim", xaxis, ": ",format(round((X.svd$d[xaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[xaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), ylab=paste0("Dim", yaxis, ": ", format(round((X.svd$d[yaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[yaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), fg = "gray", cex.lab=1.5, cex = 1.5)
  
  points(colproj[!duplicated(t(X)),c(xaxis)], colproj[!duplicated(t(X)),c(yaxis)], col = "red", pch = 17, cex = 1.5)
  points(row_sup[,c(xaxis)],row_sup[,c(yaxis)], col = "blue", pch = 23, cex = 1.5)
  points(col_sup[,c(xaxis)],col_sup[,c(yaxis)], col = "red", pch = 23, cex = 1.5)
  abline(v=0, h=0, col="gray")
  tempforlab <- rbind(rowproj[!duplicated(X),c(xaxis, yaxis)], colproj[!duplicated(t(X)),c(xaxis, yaxis)])
  autoLab(tempforlab[,c(xaxis)], tempforlab[,c(yaxis)], labels = c(rownames(rowproj[!duplicated(X), ]),rownames(colproj[!duplicated(t(X)), ])), col = c(rep("blue", nrow(rowproj[!duplicated(X), ])), rep("red", nrow(colproj[!duplicated(t(X)), ]))), cex.lab=1.5, cex = 1.5)
  text(row_sup[,c(xaxis)],row_sup[,c(yaxis)]+0.01*(axismaxy-axisminy), labels = c(roww), col = c(rep("blue", length(roww))), cex.lab=1.5, cex = 1.5)
  text(col_sup[,c(xaxis)],col_sup[,c(yaxis)]+0.01*(axismaxy-axisminy), labels = c(col), col = c(rep("red", length(col))), cex.lab=1.5, cex = 1.5)
  
} else{
  
  plot(rowproj[,c(xaxis)], rowproj[,c(yaxis)], pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
       ylim=c(axisminy,axismaxy),
       xlab=paste0("Dim", xaxis, ": ",format(round((X.svd$d[xaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[xaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), ylab=paste0("Dim", yaxis, ": ", format(round((X.svd$d[yaxis]), 3), nsmall = 3), " (", format(100*round((X.svd$d[yaxis]^2)/sum(X.svd$d^2), 3), nsmall = 1), "%)"), fg = "gray", cex.lab=1.5, cex = 1.5)
  
  points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 1.5)
  points(row_sup[,c(xaxis)],row_sup[,c(yaxis)], col = "blue", pch = 23, cex = 1.5)
  points(col_sup[,c(xaxis)],col_sup[,c(yaxis)], col = "red", pch = 23, cex = 1.5)
  abline(v=0, h=0, col="gray")
  tempforlab <- rbind(rowproj[,c(xaxis, yaxis)], colproj[,c(xaxis, yaxis)])
  autoLab(tempforlab[, c(xaxis)], tempforlab[, c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(rowproj)), rep("red", nrow(colproj))), cex.lab=1.5, cex = 1.5)
  text(row_sup[,c(xaxis)],row_sup[,c(yaxis)]+0.01*(axismaxy-axisminy), labels = c(roww), col = c(rep("blue", length(roww))), cex.lab=1.5, cex = 1.5)
  text(col_sup[,c(xaxis)],col_sup[,c(yaxis)]+0.01*(axismaxy-axisminy), labels = c(col), col = c(rep("red", length(col))), cex.lab=1.5, cex = 1.5)
}
dev.off() 

