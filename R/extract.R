#' Infer Signature Proportions
#'
#' Use known signature list to find the most likely exposures in samples
#'
#' @param object Object of class \code{DeepSig}
#' @param method Refitting method; \code{mle} for maximum likelihood (default) or
#'               \code{mutCone} for mutationalCone.
#' @param screen Logical matrix of dimension \code{nsample x nsignature} whose
#'        rows indicate signatures to be included in refitting for samples
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
#' @import Rcpp
#' @useDynLib DeepSig
#' @examples
#' data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt', package='DeepSig'))
#' b <- DeepSig(data)
#' b <- extractSig(b, progress.bar = TRUE)
#' b_pv <- extractSig(b, compute.pval = TRUE, progress.bar = TRUE)
#' @export
extractSig <- function(object, method = 'mle', screen = NA, 
                       itmax = 1000, tol = 1e-4, min.tmb = 1,
                       compute.pval = FALSE, nperm = 1000, progress.bar = FALSE,
                       pvtest = 'permutation', cosmic = FALSE, Kmin = 2, K = 2, alpha = 1,
                       Kmax = 30, nrun =10, useC = TRUE, initializer = 'random', 
                       nboot = 0, verbose = 2, 
                       fix.h = NA, s2a = NA, 
                       sig.collapse = NA, ...){

  if(Kmax < 2 ) stop('Kmax must be at least 2')
  if(!is(object, 'DeepSig')) stop('object is not of class DeepSig')
  if(!method %in% c('nmf','bnmf','mle','mutcone')) stop('Unknown method')
  if(method != 'mle' & pvtest == 'lrt') stop('Likelihood ratio test is possible only with MLE')
  if(method == 'mutcone' & compute.pval) stop('P-value in mutcone not implemented')

  spectrum <- catalog(object)
  ref <- signat(object)
  S <- colnames(ref)
  nref <- NCOL(ref)
  mut.load <- tmb(object)
  nsample <- length(tmb(object))
  
  if(is(screen,'data.frame')) screen <- as.matrix(screen)
  if(is(screen,'matrix')){ 
    if(NROW(screen)!=nsample | NCOL(screen)!=nref)
      stop('screeen matrix has incorrect dimension')
  } else{
    screen <- matrix(TRUE, nrow=nsample, ncol=nref)
  }
  
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
  
  if(method %in% c('nmf','bnmf')){
    if(compute.pval) cat('Computing exposures ...\n',sep='')
    if(method=='nmf'){
      sigs <- NULL
      nmf <- nmf(catalog = spectrum, denovo = TRUE, sig = sigs, expos = expos,
                 K = K, progress.bar = progress.bar, nrun = nrun, useC = useC, verbose = verbose -1, 
                 initializer = initializer, fix.h = fix.h, ...)
      Ha <- nmf$H
    } else{ # bnmf
      sig <- colnames(signat(object))
      object <- bnmf(object, ranks=seq(Kmin,Kmax), nrun = nrun, useC = useC, 
                     initializer = initializer, verbose = verbose-1, alpha = alpha, 
                     progress.bar = progress.bar, ...)
      W <- signat(object)
      H <- expos(object)
      sig <- paste0('S',seq(NCOL(W)))
      colnames(W) <- rownames(H) <- sig
      for(k in seq_along(misc(object)$dW)){
        if(!is.null(misc(object)$dW[[k]])) colnames(misc(object)$dW[[k]]) <- sig[seq(1, Kmin+k-1)]
        if(!is.null(misc(object)$dH[[k]])) rownames(misc(object)$dH[[k]]) <- sig[seq(1, Kmin+k-1)]
      }
      signat(object) <- W
      expos(object) <- H
      return(object)
    }
    
    H <- t(Ha)
    expos(object) <- t(H)
    logLik(object) <- nmf$logLik
    
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
      sci <- screen[i,]
      refi <- ref[, sci, drop=F]
      if(sum(sci)>1)
        fi <- fitMLE(x = spec, ref = refi, itmax = itmax, tol = tol)
      else
        fi <- list(h=1, loglik=NA)
      if(all(sci))
        hi <- fi$h
      else{
        hi <- rep(0, nref)
        names(hi) <- S
        hi[sci] <- fi$h
      }
      h[i, ] <- hi
      L[i] <- fi$loglik
    }
    else
      h[i, ] <- hi <- mutationalCone(catalog = spectrum[, i, drop=F], 
                                     signature = ref, normalize = TRUE)
    if(compute.pval){
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
