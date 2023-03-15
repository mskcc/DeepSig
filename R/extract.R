#' Infer Signature Proportions
#'
#' Use known signature list to find the most likely exposures in samples
#'
#' @param object Object of class \code{DeepSig}
#' @param method Refitting method; \code{mle} for maximum likelihood (default) or
#'               \code{mutCone} for mutationalCone.
#' @param itmax Maximum number of iterations for maximum likelihood estimate
#' @param tol Tolerance for convergence
#' @param min.tmb Minimum number of mutations in each sample. If \code{tmb} is less,
#'        \code{NA} will be returned for the sample.
#' @param compute.pval Estimate p-values
#' @param nperm Number of permutations
#' @param progress.bar Display progress bar
#' @param pvtest Algorithm for p-value computation; 
#'            \code{c('permutation','lrt','x.permutation')} for permutation resampling of 
#'            signatures, likelihood ratio test (asymptotic formula), or
#'            permutation of count data.
#' @param deviance Compute deviance for each signature
#' @param ... Other parameters for \code{denovo} with \code{method = 'hnmf'}
#' @import Rcpp
#' @useDynLib DeepSig
#' @examples
#' data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt', package='DeepSig'))
#' b <- DeepSig(data)
#' b <- extractSig(b, progress.bar = TRUE)
#' b_pv <- extractSig(b, compute.pval = TRUE, progress.bar = TRUE)
#' @export
extractSig <- function(object, method = 'mle', itmax = 1000, tol = 1e-4, min.tmb = 2,
                       compute.pval = FALSE, nperm = 1000, progress.bar = FALSE,
                       pvtest = 'permutation', cosmic = FALSE, Kmin = 2, K = 2, alpha = 1,
                       Kmax = 30, nrun =10, useC = TRUE, initializer = 'random', 
                       nboot = 0, verbose = 2, pg.pval = FALSE, deviance = FALSE, 
                       fix.h = NA, s2a = NA, pg.measure = 'cosim', 
                       sig.collapse = NA, ...){

  if(Kmax < 2 ) stop('Kmax must be at least 2')
  if(!is(object, 'DeepSig')) stop('object is not of class DeepSig')
  if(!method %in% c('nmf','hnmf','bnmf','mle','mutcone')) stop('Unknown method')
  if(method != 'mle' & pvtest == 'lrt') stop('Likelihood ratio test is possible only with MLE')
  if(method == 'mutcone' & compute.pval) stop('P-value in mutcone not implemented')

  spectrum <- catalog(object)
  ref <- signat(object)
  S <- colnames(ref)
  nref <- NCOL(ref)
  mut.load <- tmb(object)
  nsample <- length(tmb(object))
  
  is.s2a <- sum(!is.na(s2a))>0
  if(is.s2a) 
    if(!all(S %in% names(s2a))) stop(paste0(S[!S %in% names(s2a)],' not in signature list'))

  nt <- rownames(ref)
  nnt <- length(nt)
  dH <- NULL
  idx <- match(nt, rownames(spectrum))
  if(sum(is.na(idx)) > 0 | sum(duplicated(idx)) > 0)
    stop('Mutation types in spectrum do not match reference.')
  spectrum <- spectrum[idx, , drop = FALSE]  # rearrange rows to match reference
  
  if(method %in% c('nmf','hnmf','bnmf')){
    if(compute.pval) cat('Computing exposures ...\n',sep='')
    if(method=='hnmf'){
      K <- NULL
      sigs <- ref
    } else if(method=='nmf'){
      if(initializer %in% c('pgback','restart')){ 
        sigs <- ref
        expos <- expos(object)
      }
      else sigs <- NULL
    }
    if(method %in% c('nmf','hnmf')){
      nmf <- nmf(catalog = spectrum, denovo = (method=='nmf'), sig = sigs, expos = expos,
                 K = K, progress.bar = progress.bar, nrun = nrun, useC = useC, verbose = verbose -1, 
                 initializer = initializer, fix.h = fix.h, ...)
      if(deviance){
        
        x0 <- nmf$W %*% nmf$H[, nsample, drop=F]
        cs <- cosineSimilarity(spectrum[, nsample], x0, diag=TRUE)
        
        if(is.s2a){
          A0 <- s2a[S]
          A <- A0[!duplicated(A0)]
          nA <- length(A)
          tmp <- aggregate(nmf$H, by=list(A0), FUN=sum)
          Ha <- tmp[,-1]
          rownames(Ha) <- tmp[,1]
        } else{
          nA <- K
          A <- A0 <- S
          Ha <- nmf$H
        }
        dH <- matrix(0, nrow=nA, ncol=nsample)
        
        for(k in seq(nA)){
          if(verbose > 1) cat('Deviance for ', A[k], '...\n',sep='')
          if(initializer=='pgback') sj <- nsample
          else sj <- seq(nsample)
          flag <- A0==A[k]
          for(j in sj){
            fix.h2 <- fix.h
            rownames(fix.h2) <- S
            fix.h2[flag, j] <- 1 - fix.h2[flag, j]
            if(initializer=='pgback')
              nmf1 <- nmf(catalog = spectrum, denovo = (method=='nmf'), sig = sigs, expos = expos,
                     K = K, progress.bar = progress.bar, nrun = nrun, useC = useC, verbose = verbose-1, 
                     initializer = initializer, fix.h = fix.h2, ...)
            else
              nmf1 <- nmf(catalog = spectrum, denovo = (method=='nmf'), sig = nmf$W, expos = nmf$H,
                          K = K, progress.bar = progress.bar, nrun = nrun, useC = useC, verbose = verbose-1, 
                          initializer = initializer, fix.h = fix.h2, ...)
            if(pg.measure=='deviance')
              dH[k, j] <- 2*(nmf$logLik - nmf1$logLik)
            else{
              x1 <- nmf$W %*% nmf1$H[, nsample, drop=F]
              cs1 <- cosineSimilarity(spectrum[, nsample], x1, diag=TRUE)
              dH[k, j] <- cs1 - cs
            }
          }
        }
      } else Ha <- nmf$H
    } else{ # bnmf
      sig <- colnames(signat(object))
      object <- bnmf(object, ranks=seq(Kmin,Kmax), nrun = nrun, useC = useC, 
                     initializer = initializer, verbose = verbose-1, alpha = alpha, 
                     progress.bar = progress.bar, ...)
      W <- signat(object)
      H <- expos(object)
      if(initializer %in% c('restart','pgback'))
        colnames(W) <- rownames(H) <- sig
      else
        sig <- paste0('S',seq(NCOL(W)))
      colnames(W) <- rownames(H) <- sig
      for(k in seq_along(misc(object)$dW)){
        if(!is.null(misc(object)$dW[[k]])) colnames(misc(object)$dW[[k]]) <- sig[seq(1, Kmin+k-1)]
        if(!is.null(misc(object)$dH[[k]])) rownames(misc(object)$dH[[k]]) <- sig[seq(1, Kmin+k-1)]
      }
      signat(object) <- W
      expos(object) <- H
      if(initializer=='pgback'){
        if(!deviance){
          lspec <- spectrum
          px <- spectrum[,nsample]
          if(nboot>0){
            for(l in seq(nboot)){
              if(pg.pval)
                lspec[,nsample] <- lspec[,nsample][sample(NROW(lspec))]
              else  # bootstrapping
                lspec[,nsample] <- rmultinom(n = 1, size=sum(px), prob=px)
              boot <- object
              catalog(boot) <- lspec
              boot <- bnmf(boot, ranks=seq(Kmin,Kmax), nrun=nrun, useC=useC,
                     initializer=initializer, verbose = verbose, alpha=alpha, ...)
              dH <- rbind(dH, expos(boot)[,nsample])
            }
            if(pg.pval)
              dH <- rowMeans(t(dH) >= H[,nsample])
            else
              dH <- apply(dH, 2, sd)/H[,nsample]
          }
        } else{ # pg.dev
          if(Kmin!=Kmax) stop('bnmf with pg.dev requires single K')
          for(k in seq(Kmax)){
            b <- bnmf(object, ranks=Kmax, nrun=nrun, useC=useC, initializer=initializer,
                      verbose=verbose, alpha=alpha, fix.h = k, ...)
            dH[k] <- 2*(misc(object)$measure[1,2] - misc(b)$measure[1,2])*nnt*nsample
          }
        }
        misc(object)$cv <- dH
      }
      return(object)
    }
    
#    H <- t(nmf$H)
    H <- t(Ha)
    if(method != 'hnmf'){ 
      W <- nmf$W
      if(is.null(rownames(W))) rownames(W) <- nt
      if(is.null(colnames(W))) colnames(W) <- S <- paste0('S', seq(ncol(W)))
      if(!is.s2a){
        if(is.null(colnames(H))) colnames(H) <- S
        if(is.null(rownames(H))) rownames(H) <- colnames(spectrum)
      }
      signat(object) <- W
      if(initializer=='pgback' & !deviance){
        lspec <- spectrum
        px <- spectrum[, nsample]
        for(l in seq(nboot)){
          lspec[, nsample] <- rmultinom(n = 1, size = sum(px), prob = px)
          lnmf <- nmf(catalog = lspec, denovo = (method=='nmf'), sig = sigs, expos = expos,
                 K = K, progress.bar = progress.bar, nrun = nrun, useC = useC, ...)
          dH <- rbind(dH, lnmf$H[,nsample])
        }
        dH <- apply(dH, 2, sd)/H[nsample,]
      } 
    }
    expos(object) <- t(H)
    logLik(object) <- nmf$logLik
    if(initializer=='pgback'){
      if(is.s2a) nr <- nA
      else nr <- nref
      dhs <- matrix(dH[,nsample], nrow=nr, ncol=1)
      rownames(dhs) <- A
      misc(object)[['dH']] <- list(dhs)
    } else if(deviance){
      rownames(dH) <- colnames(H)
      colnames(dH) <- rownames(H)
      misc(object)[['dH']] <- list(dH)
    }
    
    if(compute.pval & method == 'hnmf'){   # estimate p-values by permutation resampling
      cat('Estimating p-values ...\n',sep='')
      if(progress.bar) pb <- txtProgressBar(style = 3)
      pv <- matrix(0, nrow = nsample, ncol = nref)
      rownames(pv) <- colnames(spectrum)
      colnames(pv) <- colnames(ref)
      for(l in seq(nperm)){
        if(pvtest == 'x.permutation'){
          spec <- apply(spectrum, 2, function(x){sample(x, size = length(x), replace = FALSE)})
          rownames(spec) <- rownames(spectrum)
          nmf <- nmf(catalog = spec, denovo = FALSE, sig = ref, progress.bar = FALSE, 
                   nrun = 1, verbose = FALSE, ...)
          pv <- pv + (t(nmf$H) >= H)
        } else if(pvtest == 'permutation'){
          rref <- ref
          for(k in seq(nref)){
            rref[, k] <- rref[sample(nnt),k]
            names(rref[,k]) <- nt
            nmf <- nmf(catalog = spectrum, denovo = FALSE, sig = rref, progress.bar = FALSE, 
                       nrun = 1, verbose = FALSE, ...)
            pv[, k] <- pv[, k] + (t(nmf$H)[, k] >= H[, k])
          }
        }
        if(progress.bar) setTxtProgressBar(pb, l/nperm)
      }
      pvalue(object) <- pv/nperm
      if(progress.bar) close(pb)
    }
    return(object)
  }
  
  h <- matrix(0, nrow = nsample, ncol = nref)
  signatures <- colnames(ref)
  colnames(h) <- signatures
  rownames(h) <- colnames(spectrum)
  if(compute.pval) pv <- h
  
  if(progress.bar) pb <- txtProgressBar(style = 3)
  L <- rep(0, nsample)
  names(L) <- colnames(spectrum)
  for(i in seq(nsample)){
    if(mut.load[i] < min.tmb){ 
      h[i, ] <- rep(NA, nref)
      if(compute.pval) pv[i, ] <- rep(NA, nref)
      next()
    }
    spec <- spectrum[, i]
    if(method=='mle'){
      fi <- fitMLE(x = spec, ref = ref, itmax = itmax, tol = tol)
      h[i, ] <- hi <- fi$h
      L[i] <- fi$loglik
    }
    else
      h[i, ] <- hi <- mutationalCone(catalog = spectrum[, i, drop=F], 
                                     signature = ref, normalize = TRUE)
    if(compute.pval){
      if(pvtest == 'lrt'){
        for(k in seq(nref)){
          fm <- fitMLE(x = spec, ref = ref[,-k], itmax = itmax, tol = tol)
          q = 2*mut.load[i]*(L[i] - fm$loglik)
#         pv[i, k] <- pchisq(q, df = 1, lower.tail = FALSE)
          pv[i,k] <- q  # just store the deviance
        }
      } else if(pvtest == 'dcs'){  # Delta CS
          for(k in seq(nref)){
            xh0 <- ref %*% h[i, ]
            cs0 <- sum(spec *xh0)/sqrt(sum(spec^2)*sum(xh0^2))
            sk <- colnames(ref)[k]
            sk.out <- sk.name <- NULL
            if(!any(is.na(sig.collapse))){
              sk.out <- unlist(lapply(sig.collapse, 
                                      function(x){if(sk %in% x) x else NULL}))
              if(!is.null(sk.out) & !any(sk.out %in% colnames(ref)))
                stop(paste0(sk.out[!sk.out %in% colnames(ref)]),' not in signature lists')
              sk.name <- names(sig.collapse)[unlist(lapply(sig.collapse, function(x){sk %in% x}))]
            } 
            if(is.null(sk.out)) sk.out <- colnames(ref)[k]
            sk.in <- colnames(ref)[!colnames(ref) %in% sk.out]
            fh <- fitMLE(x = spec, ref = ref[, sk.in, drop=F], itmax = itmax, tol = tol) # fit without sig k
            xh <- ref[, sk.in, drop=F] %*% fh$h
            cs <- sum(spec * xh)/sqrt(sum(spec^2)*sum(xh^2))
            pv[i, k] <- cs0 - cs   # negative changes in CS upon omission of signature k
            if(length(sk.name) > 0)
              colnames(pv)[colnames(pv) %in% sk.out] <- sk.name
          }
      } else{
        perm <- matrix(0, nrow=nref, ncol=nperm)
        rownames(perm) <- signatures
        for(l in seq(nperm)){   # samples under null hypothesis
          rref <- ref
          if(pvtest == 'x.permutation'){
            rsp <- spec
            rsp <- rsp[sample(length(rsp))]
            names(rsp) <- nt
            fh <- fitMLE(x = rsp, ref = ref, itmax = itmax, tol = tol)
            perm[, l] <- fh$h
            if(progress.bar) if(l %% 100 == 0)
              setTxtProgressBar(pb, ((i - 1)*nperm + l)/nsample/nperm)
          } else if(pvtest == 'permutation'){
            for(k in seq(nref)){
              rref[, k] <- rref[sample(nnt),k]
              names(rref[,k]) <- nt
              fk <- fitMLE(x = spec, ref = rref, itmax = itmax, tol = tol)
              perm[k, l] <- fk$h[k]
              if(progress.bar) if(l %% 10 == 0)
                setTxtProgressBar(pb, ((i - 1)*nperm*nref + (l - 1)*nref + k)/nsample/nperm/nref)
            }
          } else stop('Unknown pvtest')
        }
        pv[i, ] <- rowSums(perm >= h[i, ]) / nperm
      }
    }
    if(progress.bar) setTxtProgressBar(pb, i/nsample)
  }
  if(progress.bar) close(pb)

  expos(object) <- h
  if(compute.pval) pvalue(object) <- pv
  if(method=='mle')
    logLik(object) <- L

  return(object)
}

#' Maximum likelihood inference of signature proportions
#' 
#' @param x Vector of observed proportions of mutation contexts
#' @param ref Matrix of reference signature proportions with mutation types in rows
#'            and signatures in columns
#' @param itmax Maximum number of iterations
#' @param tol Tolerance of convergence
#' @return List of exposure vector \code{p} and \code{loglik}
#' @export
fitMLE <- function(x, ref, itmax = 1000, tol = 1e-4){

  num_sigs <- NCOL(ref)
  num_muts <- sum(x)
  x <- x / num_muts
  x0 <- gtools::rdirichlet(n = 1, alpha = rep(10, num_sigs)) #initial guess
  p <- mlestimate(x, x0, ref, Itmax=itmax, Tol=tol)
  h <- p$x^2/sum(p$x^2)
  names(h) <- colnames(ref)
  
  return(list(h=h, loglik = p$lkh))
}
