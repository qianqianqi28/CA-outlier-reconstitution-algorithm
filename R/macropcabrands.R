# This is for CA: Figure 1 (c) in the Section Empirical studies/Results; and Figure 2

rm(list=ls()) 
library(readxl)
library(xtable)
library(FactoMineR)
library(cellWise)
library("ggpubr")


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

MacroPCApar0 <- list(scale = FALSE, center = rep(0, d), alpha = 0.97)

MacroPCA.out <- MacroPCA(S, k = 2, MacroPCApars = MacroPCApar0)

# DDC

# outlying rows in the DDC
indrows <- MacroPCA.out$DDC$indrows
rownames(X)[indrows]

# outlying cells in the DDC
length(MacroPCA.out$DDC$indcells)

savefigureddc <- "brandsDDCcelloutlier"

pdf(paste0("plots/", savefigureddc, ".pdf"), width = 12,
    height = 12)
ggpmacropca     <- cellMap(R = MacroPCA.out$DDC$stdResid,                        
                           indcells = MacroPCA.out$DDC$indcells,
                           indrows = MacroPCA.out$DDC$indrows,
                           mTitle = "",
                           rowtitle = "",
                           columntitle = "",
                           sizetitles = 2.4)

My_Theme = theme(
  axis.title.x = element_text(size = 10),
  axis.text.x = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 10))
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

Flip <- "no flip"
Checkoverlap <- TRUE

xaxis <- 1
yaxis <- 2

savefigure <- "macbrandscasym"


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


pdf(paste0("plots/", savefigure, ".pdf"), width = 12,
    height = 12)

if (Checkoverlap == TRUE){
  plot(rowproj[!duplicated(X), c(xaxis)], rowproj[!duplicated(X), c(yaxis)], pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
       ylim=c(axisminy,axismaxy), xlab=paste0("Dim", xaxis), ylab=paste0("Dim", yaxis), fg = "gray", cex.lab=1.5, cex = 1.5)
  
  points(colproj[!duplicated(t(X)), c(xaxis)], colproj[!duplicated(t(X)), c(yaxis)], col = "red", pch = 17, cex = 1.5)
  abline(v=0, h=0, col="gray")
  tempforlab <- rbind(rowproj[!duplicated(X), ], colproj[!duplicated(t(X)), ])
  
  autoLab(tempforlab[,c(xaxis)], tempforlab[,c(yaxis)], labels = c(rownames(rowproj[!duplicated(X),]),rownames(colproj[!duplicated(t(X)),])), col = c(rep("blue", nrow(X[!duplicated(X),])), rep("red", ncol(X[!duplicated(t(X)),]))), cex.lab=1.5, cex = 1.5)

  } else {
    
  plot(rowproj[,c(xaxis)], rowproj[,c(yaxis)], pch = 20, asp=1, col = "blue", xlim=c(axisminx,axismaxx),
       ylim=c(axisminy,axismaxy), xlab=paste0("Dim", xaxis), ylab=paste0("Dim", yaxis), fg = "gray", cex.lab=1.5, cex = 1.5)
  
  points(colproj[,c(xaxis)], colproj[,c(yaxis)], col = "red", pch = 17, cex = 1.5)
  
  abline(v=0, h=0, col="gray")
  
  tempforlab <- rbind(rowproj, colproj)
  
  autoLab(tempforlab[, c(xaxis)], tempforlab[, c(yaxis)], labels = c(rownames(rowproj),rownames(colproj)), col = c(rep("blue", nrow(X)), rep("red", ncol(X))), cex.lab=1.5, cex = 1.5)
}

dev.off()







