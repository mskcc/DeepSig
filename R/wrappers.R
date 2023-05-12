#' Wrapper Function
#' 
#' One-step execution of exposure inference, filtering, and output
#' 
#' @param catalog Catalog input file
#' @param signat Reference signature choices, \code{c('v2','SA','SP')} for sCOSMIC v2, 
#'        SigAnalyzer, and sigProfiler, respectively
#' @param nperm No. of permutations 
#' @param alpha False positive rate threshold
#' @param output Output file
#' @export
run_tSig <- function(catalog, signat = 'v2', nperm = 1000, alpha = 0.05, output = 'exposure.txt'){
  
  if(signat=='SA') ref <- system.file('extdata', 'cosmic_SigAnalyzer_SBS_signatures.txt', package = 'DeepSig')
  else if(signat=='SP') ref <- system.file('extdata', 'cosmic_sigProfiler_SBS_signatures.txt', package = 'DeepSig')
  else if(signat=='v2') ref <- system.file('extdata', 'cosmic_snv_signatures_v2.txt', package = 'DeepSig')
  else stop('Unknown choice for signat')
  b <- DeepSig(catalog, signat = ref)
  b <- extractSig(b, compute.pval = TRUE, nperm = nperm, progress.bar = TRUE)
  e <- expos(b)
  pv <- pvalue(b)
  for(i in seq(NROW(e))) for(j in seq(1, NCOL(e))){
    if(is.na(e[i,j])) next()
    if(pv[i,j] >= alpha) e[i,j] <- 0
  }
  e <- cbind(data.frame(Sample.Name = rownames(e), Number.of.Mutations = tmb(b)), e)
  write.table(e, file = output, row.names = F, sep = '\t', quote=F)
}

#' Run multiple decomposition with bootstrapping
##' @export
#runBootstrap <- function(catalog, method = 'nmf', verbose = 2, nboot = 100, K = 2, nrun = 100, 
#                         progress.bar = TRUE, Kmin = 2, Kmax = 2, reorder = TRUE, ...){
#  
#  nt <- trinucleotides()
#  nsamp <- NCOL(catalog)
#  M <- colSums(catalog)
#  b <- list()
#  for(i in seq(nboot)){
#    if(verbose >= 2)
#      cat('Bootstrapped run #', i, '...\n', sep='')
#    xcat <- matrix(0, nrow = 96, ncol = nsamp)
#    rownames(xcat) <- nt
##    colnames(xcat) <- colnames(catalog)
#    for(j in seq(nsamp)){
#      tmp <- sample(factor(nt, levels = nt), size = M[j], replace = TRUE, prob = catalog[,j])
#      xcat[,j] <- table(tmp)
#    }
#    b[[i]] <- DeepSig(data = xcat)
##    if(method == 'nmf')
#     b[[i]] <- extractSig(b[[i]], method = method, K = K, nrun = nrun, verbose = verbose - 1, 
#                         progress.bar = progress.bar, ...)
#    else if(method == 'bnmf')
#      b[[i]] <- extractSig(b[[i]], method = method, Kmin = Kmin, Kmax = Kmax, nrun = nrun, verbose = verbose-1,
#                           progress.bar = progress.bar, ...)
#    else stop(paste0('Method ',method,' not allowed'))
##    if(reorder & i > 1){
#      s1 <- signat(b[[1]])
#      si <- signat(b[[i]])
#      cs <- cosineSimilarity(si, s1, diag = FALSE)
#      idx <- clue::solve_LSAP(cs, maximum = TRUE)
#      if(method=='bnmf') K <- Kmax
#      idx <- match(seq(K), idx)
#      b[[i]] <- reorderSig(b[[i]], order = idx)
#    }
#  }
#  
#  return(b)
#}
