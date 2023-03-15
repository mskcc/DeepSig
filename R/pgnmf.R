#' piggyback NMF
#' 
#' Enables de novo inference of single sample using pre-inferred exposure 
#' and signature sets for a large data set (stock)
#' @param x Catalog matrix to be analyzed
#' @param method \code{'nmf'} or \code{'bnmf'} for inference engine
#' @param stock \code{DeepSig} object containing the stock solution
#' @param filter Filter exposure using CV cutoff \code{cutoff}
#' @param cutoff CV cutoff matrix for filtering
#' @param discrete Convert proportions into multinomial counts
#' @param alpha Dirichlet alpha for initial guess of exposure vector
#' @param nboot No. of bootstrapping runs
#' @param pg.lr Compute deviance for each signature
#' @param pg.pval Estimate p-values by permutation
#' @param prefilter List of allowed cancer type to be used for prefiltering; 
#'        if \code{NA}, skip prefiltering
#' @return List of \code{mean} and coefficient of variation \code{cv} (sd/mean) of exposure matrices 
#' @export
pgnmf <- function(x, stock = NULL, method = 'bnmf', filter = FALSE, cutoff = NULL, 
                  discrete = FALSE, progress.bar = TRUE, 
                  nrun = 1, verbose = 1, alpha = 1, nboot = 0, pg.dev = TRUE, 
                  pg.pval = FALSE, xtype=NA, prefilter = NA, tspec = FALSE, 
                  s2a = NA , ...){
  
  if(pg.dev & pg.pval) stop('Deviance or p-value by permuation but not both')
  
  if(is.null(stock)){
    fl <- system.file('extdata/pig_K11.rds', package = 'DeepSig')
    if(!file.exists(fl)) stop('Default stock object does not exist')
    stock <- readRDS(fl)
  }
  
  if(sum(!is.na(prefilter))>0) if(sum(is.na(xtype))>0) stop('prefilter requires xtype')
  if(tspec) if(sum(is.na(xtype))>0) stop('prefilter requires xtype')

  if(!method %in% c('bnmf','nmf')) stop(paste0('Uknown method: ',method))
  if(!tspec){
    stock.x <- catalog(stock)
    stock.h <- expos(stock)
    stock.w <- signat(stock)
    sig <- rownames(stock.h)
    K <- nrow(stock.h)
    N <- ncol(stock.x)
  }
  
  m <- colSums(x)
  n <- ncol(x)
  
#  dat <- do.call(rbind, lapply(seq(n), function(k){
  if(tspec) dat <- vdat <- list()
  else dat <- vdat <- NULL
  if(progress.bar) pb <- txtProgressBar(style = 3)
  
  if(!tspec){ 
    Fix.h <- matrix(0, nrow=K, ncol=n)
    colnames(Fix.h) <- colnames(x)
    rownames(Fix.h) <- sig
  }
    
  for(k in seq(n)){
    if(tspec){
      if(!xtype[k] %in% names(stock)){
        stop(paste0('stock object for', xtype[k],' not supplied'))
      }
      stk <- stock[[xtype[k]]]
      stock.x <- catalog(stk)
      stock.h <- expos(stk)
      stock.w <- signat(stk)
      sig <- rownames(stock.h)
      K <- nrow(stock.h)
      N <- ncol(stock.x)
    }
    xcat <- cbind(stock.x, data.frame(xk = x[,k]))
    b <- DeepSig(data = xcat, signat = stock.w)
#    expos(b) <- cbind(stock.h, data.frame(xk = rep(1/K, K)))
    pk <- gtools::rdirichlet(n = 1, alpha = rep(alpha, K))[1,]
#   pk <- 1/K
    expos(b) <- cbind(stock.h, data.frame(xk = pk))
#   tmb(b)[n] <- m[k]

    if(method=='bnmf')
      b <- extractSig(b, method = 'bnmf', nrun = nrun, Kmin = K, Kmax = K, 
                    initializer = 'pgback', verbose = verbose -1 , alpha = alpha, 
                    nboot = nboot, pg.pval = pg.pval, deviance = pg.dev, ...)
    else{
      fix.h <- matrix(1, nrow=K, ncol=N+1)
      if(!tspec & sum(!is.na(prefilter))>0)
        for(l in seq(K)) fix.h[l, N+1] <- xtype[k] %in% prefilter[[l]]
      b <- extractSig(b, method = 'nmf', nrun = nrun, K = K, initializer = 'pgback',
                    verbose = verbose - 1, nboot = nboot, deviance = pg.dev, 
                    fix.h = fix.h, s2a = s2a, ...)
      if(!tspec) Fix.h[, k] <- fix.h[,N+1]
    }
    if(all(is.na(s2a))){  # don't try to rotate if you are aggregating
      cs <- cosineSimilarity(stock.w, signat(b))
      lsap <- clue::solve_LSAP(cs, maximum = TRUE)
      if(!all(lsap==seq(K))){
        warning('Ref. signature scrambled')
        b <- reorderSig(b, order=lsap, label = rownames(cs))
      }
    }
    
    dx <- t(expos(b)[,N+1, drop=F])
    rownames(dx) <- colnames(x)[k]
    if(!tspec) dat <- rbind(dat, dx)
    if(length(misc(b)) > 0){
      if(method=='bnmf'){
        if(nboot==0) dv <- t(misc(b)$dH[[1]][,N+1])
        else dv <- matrix(misc(b)$cv, nrow=1)
      } else dv <- t(misc(b)$dH[[1]][,1, drop=F])
#      } else dv <- t(misc(b)$dH[,1])
      rownames(dv) <- colnames(x)[k]
      if(sum(!is.na(s2a))==0) colnames(dv) <- sig
      if(!tspec) vdat <- rbind(vdat, dv)
      else
        dat[[colnames(x)[k]]] <- data.frame(E=t(dx)[,1], Me=m[k]*t(dx)[,1], Dev=t(dv)[,1])
    }
    if(progress.bar) setTxtProgressBar(pb, value = k/n)
#  }))
  }
  if(progress.bar) close(pb)
  
  if(filter){
#   cat('Filtering exposure...')
    if(is.null(cutoff)){
      fl <- system.file('extdata/pig_K11_cutoff.txt', package = 'DeepSig')
      if(!file.exists(fl)) stop('Default stock object does not exist')
      cutoff <- read.table(fl, header = TRUE, sep = '\t')
    }
    sig <- colnames(cutoff)[-1]
    if(!all(sig == colnames(dat)))
      stop('Signature names in cutoff do not agree with exposure matrix')
    sfit <- function(m, S){
      x <- log10(m)
      if(!S %in% sig) stop(paste0(S, ' not in signature set'))
      g <- smooth.spline(x = log10(cutoff$M), y = cutoff[, S], df = nrow(cutoff))
      gamma <- predict(g, x = x)$y
      return(gamma)
    }
    
    for(i in seq(n)) for(k in seq(K)){
      g <- sfit(m = m[i], S = sig[k])
      if(vdat[i, k] >= g) dat[i, k] <- 0
    }
  }
  
  if(discrete){   # replace proportions by multinomial counts
    for(i in seq(n)){
      r <- rmultinom(n = 1, size = m[i], prob = dat[i,])
      dat[i,] <- r
    }
  }
  
  if(tspec) return(dat)
  else return(list(mean=dat, cv=vdat, fix.h = t(Fix.h)))
  
}
