reconca <- function (X, ncp = 2, ncpsele = FALSE, shrink = FALSE, initials = 0.0001, inisample = FALSE, threshold = 1e-08, maxiter = 100000, outliercell = c(), verbose=FALSE) 
{
  
  shrinkCA <- function(X, ncp = 2, shrink) {
    P <- as.matrix(X/sum(X))
    Rc <- apply(P, 2, sum)
    Rr <- apply(P, 1, sum)
    S <- t(t((P - Rr %*% t(Rc))/sqrt(Rr))/sqrt(Rc))
    svdRes <- svd(S)
    n <- nrow(X)
    p <- ncol(X)
    
    if (shrink == FALSE){
      lambda.shrinked <- svdRes$d[1:ncp]
    } else {
      sigma2 <- sum(svdRes$d[-c(1:ncp)]^2)/((n - 1) * (p - 1) - (n - 1) * ncp - (p - 1) * ncp + ncp^2)
      lambda.shrinked <- (svdRes$d[1:ncp]^2 - n * (p/min(p, (n - 1))) * sigma2)/svdRes$d[1:ncp]
    }
    
    if (ncp == 1){
      recon <- (svdRes$u[, 1] * lambda.shrinked) %*% t(svdRes$v[,  1])
      recon <- sum(X) * (t(t(recon * sqrt(Rr)) * sqrt(Rc)) + Rr %*% t(Rc))
    } else if (ncp > 1 ){
      recon <- svdRes$u[, 1:ncp] %*% (t(svdRes$v[, 1:ncp]) * lambda.shrinked)
      recon <- sum(X) * (t(t(recon * sqrt(Rr)) * sqrt(Rc)) + Rr %*% t(Rc))
    } else if (ncp == 0){
      recon <- sum(X) * (Rr %*% t(Rc))
    }
    
    
    rownames(recon) <- rownames(X)
    colnames(recon) <- colnames(X)
    res <- list(recon = recon)
    return(res)
  }
  
  X <- as.matrix(X)
  
  vector_2x3 <- do.call(paste, expand.grid(rownames(X), colnames(X), stringsAsFactors = FALSE))
  matrix_2x3 <- matrix(vector_2x3, nrow = nrow(X))
  for (i in 1:length(outliercell)){
    print(which(matrix_2x3 == outliercell[i], arr.ind = TRUE))
    X[which(matrix_2x3 == outliercell[i], arr.ind = TRUE)] = NA
  }
  
  if (sum(is.na(X)) == 0) {
    stop("No value is missing/cellwise outlier")
  }
  missing <- which(is.na(X))
  Xhat <- X
  if (inisample == FALSE){
    Xhat[missing] <- initials
  } else {
    Xhat[missing] <- sample(X[-missing], length(missing), replace = TRUE) + 1
  }
  inivalue <- Xhat[missing]
  
  nb.iter <- 1
  old <- Inf
  objective <- 0
  recon <- Xhat
  while (nb.iter > 0) {
    Xhat[missing] <- recon[missing]
    # Xhat[missing][which((Xhat[missing] <= 0))] = 0
    RXhat <- rowSums(Xhat)
    CXhat <- colSums(Xhat)
    recon <- shrinkCA(Xhat, ncp = ncp, shrink)$recon
    diff <- Xhat - recon
    diff[missing] <- 0
    objective <- sum((diff^2))
    criterion <- abs(1 - objective/old)
    if (verbose) print(c("criteria: ", criterion, "number of iterative: ", nb.iter, "xold", Xhat[missing], "xnew", recon[missing]))
    
    old <- objective
    nb.iter <- nb.iter + 1
    if (!is.nan(criterion)) {
      if ((criterion < threshold) && (nb.iter > 50)) 
        nb.iter <- 0
      if ((objective < threshold) && (nb.iter > 50)) 
        nb.iter <- 0
    }
    if (nb.iter > maxiter) {
      nb.iter <- 0
      warning(paste("Stopped after ", maxiter, " iterations"))
    }
  }
  completeObs <- X
  completeObs[missing] <- Xhat[missing]
  newlist <- list(completeObs,inivalue, Xhat[missing])
  return(newlist)
}