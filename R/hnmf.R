#' NMF Iteration
#' 
#' Non-negative matrix factorization with either fixed or variable signature matrix.
#' 
#' If \code{denovo = TRUE}, the full NMF is performed with both \code{W} and \code{H} matrix 
#' determined de novo. If \code{denovo = FALSE}, \code{W} is fixed as input and only \code{H} 
#' is determined.
#' 
#' @param catalog Catalog matrix
#' @param signat Signature matrix
#' @param denovo Full iteration of both signature and exposure; if \code{FALSE} (\code{hnmf}),
#'        only the exposure is fit (requires \code{signat} input).
#' @param K Number of signatures. 
#' @param nrun Number of independent runs to generate
#' @param verbose Verbosity level
#' @param progress.bar Display progress bar.
#' @param Itmax Maximum no. of iteration. 
#' @param Tol Tolerance for checking convergence.
#' @param fix.h Binary flags for exposure to keep fixed during inference (for deviance estimate)
#' @param clustering Perform clustering of local maxima
#' @param logLik.cut Log likelihood fractional cutoff with respect to global maximum for clustering
#'        
#' @return Object of class \code{DeepSig}.
#' 
#' @export
nmf <- function(catalog, sig = NULL, expos = NULL,
                denovo = TRUE, K = NULL, nrun = 10, verbose = 1, 
                progress.bar = FALSE, Itmax = 100000, Tol = 1e-5, a = 10, nprint = 100,
                useC = FALSE, alpha = 1, initializer = 'random', fix.h = NA,
                clustering = FALSE,
                logLik.cut = 0.05){
  
  mat <- catalog
  if(!is.matrix(mat)) mat <- as.matrix(catalog)
  nullr <- sum(Matrix::rowSums(mat)==0)
  nullc <- sum(Matrix::colSums(mat)==0)
# if(nullr>0 & denovo) stop('Input matrix contains empty rows')
  if(nullc>0) stop('Input matrix contains empty columns')
  
  nrow <- dim(mat)[1]
  ncol <- dim(mat)[2]
  
#  if(is.null(sig)) initializer <- 'random'
#  else initializer <- 'pgback'
  
  if(is.null(K)) if(!is.null(sig)) K <- NCOL(sig)
  rank <- K
  fudge <- .Machine$double.eps  # to avoid zero denominator
#  fudge <- 0
     
  if(verbose > 1) progress.bar <- FALSE  # can't show messages within progress.bar
  
  nsample <- NCOL(catalog)
# if(!all(is.na(fix.h))) if(fix.h[1] < 1 | fix.h[1] > K | fix.h[2] < 1 | fix.h[2] > nsample)
  if(!all(is.na(fix.h))){
    is.prefix <- TRUE
    if(NROW(fix.h) != K | NCOL(fix.h) != nsample |sum(!fix.h %in% c(0,1) >0)) stop('Invalid fix.h')
  } else is.prefix <- FALSE
  
  if(progress.bar) pb <- txtProgressBar(style = 3)
  
  rmax <- -Inf
  if(clustering) Ws <- Hs <- list()  # list of converged solutions
    
  llik <- rep(0, nrun)
  for(irun in seq(nrun)){
    if(verbose>1) cat('Run #',irun,':\n')
#   wh <- init(nrow, ncol, mat, rank)
    hyper <- list(aw = 0.0, ah = 0.0, bw = 1.0, bh = 1.0) # ML prior
    wh <- vb_init(nrow = nrow, ncol = ncol, mat = mat, rank = rank, hyper = hyper, initializer = initializer)
    if(initializer %in% c('pgback','restart')){
      wh$h <- expos
      if(irun > 1) wh$h[, ncol] <- gtools::rdirichlet(n = 1, alpha = rep(alpha, K))[1,]
      wh$lh <- wh$eh <- wh$h
    }
    if(!is.null(sig)){
      wh$w <- wh$lw <- wh$ew <- sig
      rownames(wh$eh) <- colnames(wh$ew)
    }
    lkold <- -Inf
    tol <- Inf
    h0 <- rep(1/K, K)
    for(it in seq_len(Itmax)){
      wh0 <- wh
      if(useC){
        wh <- vbnmf_update(mat, wh, hyper, c(fudge))
        wh$lw <- wh$w
        wh$lh <- wh$h
        if(is.prefix) wh$eh <- wh$eh * fix.h
        lk0 <- likelihood(mat, wh$ew, wh$eh)
      }
      else{
        wh <- nmf_updateR(mat, wh$ew, wh$eh, nrow, ncol, rank, denovo = denovo)
        if(is.prefix) wh$eh <- wh$eh * fix.h
        lk0 <- likelihood(mat, wh$ew, wh$eh)
      }
      if(it > 1){ 
        if(initializer=='pgback'){
          if(max(abs(h0-wh$h[,ncol(wh$h)])) < Tol) break
        } else{
           if(abs(lkold - lk0) < Tol * abs(lkold)) break
        }
      }
      h0 <- wh$h[,ncol(wh$h)]
      if(verbose>1) if(it %% nprint == 0)
        cat(it,': likelihood = ', lk0, '\n', sep='')
      lkold <- lk0
    }
    Keff <- K
    if(verbose>1)
      cat('Nsteps =', it, ', likelihood =', lk0, '\n', sep = '')
    if((irun == 1 | lk0 > rmax) & !is.na(lk0)){
      rmax <- lk0
      wmax <- wh$ew
      hmax <- wh$eh
      Kmax <- Keff
    }
    if(clustering){
      llik[irun] <- lk0
      w <- wh$ew
      cs <- colSums(w)
      w <- t(t(w)/cs)
      h <- wh$eh * cs
      h <- t(t(h)/colSums(h))
      Ws[[irun]] <- w
      Hs[[irun]] <- h
    }
    if(verbose>1) cat('Max(likelihood) =',rmax,'\n')
    if(progress.bar) setTxtProgressBar(pb, irun/nrun)
  }
  if(progress.bar) close(pb)
  
  if(clustering){ # find centroids of clusters
    Ws2 <- Hs2 <- list()
    i <- 1
    for(irun in seq(nrun)){
      df <- abs(1-llik[irun]/max(llik))
      if(df < logLik.cut){   # filter out solutions too far away from global maximum
        Ws2[i] <- Ws[irun]
        Hs2[i] <- Hs[irun]
        i <- i + 1
      }
    }
    ctrd <- centroid(K = K, Ws = Ws2, Hs = Hs2, progress.bar = progress.bar)
    wmax <- ctrd$w
    hmax <- ctrd$h
  }
  cs <- colSums(wmax)
  wmax <- t(t(wmax) / cs)
  hmax <- hmax * cs
  W <- wmax
  H <- t(t(hmax)/colSums(hmax))
  rownames(W) <- rownames(sig)
  colnames(W) <- rownames(H) <- rownames(expos)
  colnames(H) <- colnames(expos)
  
  return(list(W = W, H = H, logLik = rmax))
}

# single update step of NMF
nmf_updateR <- function(x, w, h, n, m, r, denovo = TRUE){
  
  x <- as.matrix(x)
  w <- as.matrix(w)
  h <- as.matrix(h)
  
  if(denovo){
    up <- w*((x / (w %*% h)) %*% t(h))
    hj <- rowSums(h)
    down <- matrix(rep(hj, n), nrow = n, ncol = r, byrow = TRUE)
    w2 <- up / down
    w2[w2 < .Machine$double.eps] <- .Machine$double.eps
  }
  
  up <- h*(t(w) %*% (x/(w %*% h)))
  if(denovo) w <- w2
  
  wi <- colSums(w)
  down <- matrix(rep(wi, m), nrow = r, ncol = m, byrow = FALSE)
  h <- up / down
  h[h < .Machine$double.eps] <- .Machine$double.eps

  return(list(ew=w, eh=h))
}

# initialize w and h
init <- function(nrow, ncol, mat, rank, max = 1.0){
  w <- matrix(stats::runif(n = nrow * rank), nrow = nrow, ncol = rank)
  h <- matrix(stats::runif(n = rank * ncol), nrow = rank, ncol = ncol)
  rownames(w) <- rownames(mat)
  if(is.null(colnames(w))) 
    colnames(w) <- rownames(h) <- paste0('S',seq_len(rank))
  colnames(h) <- colnames(mat)
  
  return(list(ew = w, eh = h))
}

likelihood <- function(mat, w, h, per.element = FALSE){
  
  wh <- as.vector(w %*% h)
  amat <- as.vector(mat)
#  x <- sum(amat * log(wh) - wh)
  wh2 <- wh[wh > 0]
  x <- sum(amat[wh > 0] * log(wh2) -wh2)
  z <- amat[amat > 0]
  x <- x + sum(-z * log(z) + z)
  if(per.element) x <- x / nrow(mat) / ncol(mat)
  
  return(x)
}

# cluster (W,H) matrix solutions into k groups
centroid <- function(K, Ws, Hs, itmax = 1000, Tol = 1e-5, progress.bar = FALSE){
  
  kcent <- matrix(runif(96*K), nrow = 96, ncol = K)
  kcent <- t(t(kcent)/colSums(kcent))   # random initial centroids
  rownames(kcent) <- nt <- trinucleotides()
  
  nsamp <- length(Ws)
  N <- NCOL(Hs[[1]])
  df <- rep(Inf, K)
  
  for(it in seq(itmax)){
    cat('Finding centroids: iteration # ',it,'...\n',sep='')
    kcent2 <- matrix(0, nrow = 96, ncol = K)
    hcent <- matrix(0, nrow = K, ncol = N)
    dissim <- list()
    for(k in seq(K)) 
      dissim[[k]] <- matrix(0, nrow=nsamp, ncol=nsamp)
    Idx <- matrix(0, nrow=nsamp, ncol=K)
    if(progress.bar) pb <- txtProgressBar(style=3)
    for(i in seq(nsamp)){
      S <- Ws[[i]]
      rownames(S) <- rownames(kcent) <- nt
      cs <- cosineSimilarity(S, kcent)
      lsap <- clue::solve_LSAP(x = cs, maximum = TRUE)
      idx <- match(seq(K), lsap)
      kcent2 <- kcent2 + Ws[[i]][,idx]
      hcent <- hcent + Hs[[i]][idx,]
      Idx[i,] <- idx
      for(k in seq(K)) for(j in seq(nsamp))
          dissim[[k]][i,j] <- 1 - 
            sum(Ws[[i]][,idx[k]] * Ws[[j]][,idx[k]])/
              sqrt(sum(Ws[[i]][,idx[k]]^2)*sum(Ws[[j]][,idx[k]]^2))
      if(progress.bar) setTxtProgressBar(pb, i/nsamp)
    }
    if(progress.bar) close(pb)
    kcent2 <- kcent2 / nsamp
    hcent <- hcent / nsamp
    df2 <- rep(0, K)
    for(k in seq(K)){
      sil <- cluster::silhouette(x = Idx[,k], dmatrix = dissim[[k]])
      if(!is(sil,'silhouette')) if(is.na(sil)) next()
      df2[k] <- mean(sil[,'sil_width'])
    }
    if(abs(max(df2-df)) < Tol) break()
    kcent <- kcent2
    df <- df2
  }
  
  return(list(w = kcent2, h = hcent))
}
