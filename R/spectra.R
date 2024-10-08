#' Simulate mutational spectra
#'
#' Generate simulated data using zero-inflated Poisson model
#'
#' See Omichessan et al. DOI: 10.1371/journal.pone.0221235
#'
#' @param W Input signature matrix
#' @param pzero Proportion of extra zero counts (scalar)
#' @param nmut Mean number of mutations per sample (scalar)
#' @param h Vector of average proportion of mutations in each process;
#'          length equal to the number of columns in \code{W}.
#' @param N No. of samples
#' @param dilute.ultra Dipute ultra-mutated samples to reduce bias
#' @param min.mut Minimum mutation load
#' @param distr Underlying distribution of exposure counts; 
#'        \code{c('pois', 'nbinom')}, either Poisson or negative binomial.
#' @param vmr Variance to mean ratio of negative binomial.
#' @return List of components \code{X}: simulated mutation counts;
#'         \code{W}, \code{H}: signature and exposure matrices.
#' @examples
#' set.seed(135)
#' K <- 5 # no. of processes
#' W <- t(gtools::rdirichlet(n = K, alpha = rep(1, 96)))
#' rownames(W) <- trinucleotides()
#' h <- gtools::rdirichlet(n = 1, alpha = rep(5, K))
#' names(h) <- colnames(W) <- pase0('S',seq(K))
#' x <- simulateSpectra(W = W, h = h, N = 100)
#' 
#' h <- read.table(system.file('extdata', 'hmean_breast_lung_skin_pancr.txt',
#'                             package = 'DeepSig'), header = TRUE, sep = '\t')
#' h <- t(h)['Breast',]
#' x <- simulateSpectra(W = 'SA', h = h, N = 100)
#' @export
simulateSpectra <- function(W = NULL, pzero = 0, nmut = 100, h, N = 10, dilute.ultra = FALSE, 
                            min.mut = 1, distr = 'pois', vmr = 2, Ntry = 10000){

  if(is.null(W)) W <- 'cosmic_v2'
  if(is.character(W)){
    if(!W %in% c('cosmic_v2','SA','SP')) stop(paste0('Unknown signature set: ',W))
    if(W=='cosmic_v2')
      z <- read.table(system.file('extdata', 'cosmic_snv_signatures_v2.txt',
                                  package = 'DeepSig'), header = TRUE, sep = '\t')
    if(W=='SA')
      z <- read.table(system.file('extdata', 'cosmic_SigAnalyzer_SBS_signatures.txt',
                                package = 'DeepSig'), header = TRUE, sep = '\t')
    if(W=='SP')
      z <- read.table(system.file('extdata', 'cosmic_sigProfiler_SBS_signatures.txt',
                                package = 'DeepSig'), header = TRUE, sep = '\t')
    W <- z
  }
  
  if(!is(W, 'matrix')) W <- as.matrix(W)
  m <- NROW(W)
  h <- h[h > 0]
  if(is.null(names(h))) stop('h must have element names')
  W <- W[, colnames(W) %in% names(h)]
  idx <- match(colnames(W), names(h))
  if(sum(is.na(idx)) > 0 | sum(duplicated(idx)) > 0) 
    stop('names of h do not match reference signature names')
  h <- h[idx]
  K <- length(h)
  
  if(pzero < 0 | pzero > 1) stop('Invalid pzero value')
  lsk <- h*nmut/(1 - pzero)  # process-dependent Poisson mean
  H <- matrix(0, nrow = K, ncol = N)
  for(i in seq(N)){
    itry <- 0
    while(TRUE){
      for(k in seq(1, K)){
        if(pzero > 0) p <- rbinom(n = 1, size = 1, prob = 1 - pzero)
        else p <- 1
        if(p==0) next()
        if(distr=='delta')
          H[k, i] <- lsk[k]
        else if(distr=='pois')
          H[k, i] <- rpois(n = 1, lambda = lsk[k])  # poisson
        else if(distr=='nbinom'){
          if(vmr <= 1) stop('Invalid vmr in nbinom')
          H[k, i] <- rnbinom(n = 1, size = lsk[k]/(vmr-1), mu = lsk[k])
        }
        else stop('Unknown distribution')
      }
      if(sum(H[,i]) >= min.mut) break()
      itry <- itry + 1
      if(itry > Ntry) stop('Maximum count generation tries exceeded')
    }
  }
  rownames(H) <- colnames(W)

  xmean <- W %*% H
#  while(TRUE){
#    X <- matrix(rpois(n = m*N, lambda = xmean), nrow = m, ncol= N, byrow=FALSE)
#    if(all(colSums(X) >= min.mut)) break()  # repeat until all columns are non-empty
#  }
#  rownames(X) <- rownames(W)
#  colnames(X) <- colnames(H) <- paste0('Sample_',seq(N))
  
  X <- matrix(rpois(n = m*N, lambda= xmean), nrow = m, ncol = N, byrow = FALSE)
  bad <- colSums(X) < min.mut
  X <- X[, !bad]
  H <- H[, !bad]
  rownames(X) <- rownames(W)
  colnames(X) <- colnames(H) <- paste0('Sample_',seq(NCOL(X)))
  
  if(dilute.ultra)
    X <- diluteUltraMutated(X)

  x <- list(X=X, W=W, H=H)
  return(x)
}

# Dilute ultra-mutated samples (Kim et al. DOI: 10.1038/hg.3557)

diluteUltraMutated <- function(X, maxiter=100){

  if(is.null(colnames(X)))
    colnames(X) <- seq(1, NCOL(X))
  for(i in seq(1,maxiter)){
    snv <- colSums(X)
    q1 <- quantile(snv, prob=1/4)
    q3 <- quantile(snv, prob=3/4)
    s.ultra <- colnames(X)[snv > (median(snv) + 1.5*(q3 - q1))]
    if(length(s.ultra)==0) break
    m.ultra <- as.matrix(X[, (colnames(X) %in% s.ultra)])
    colnames(m.ultra) <- s.ultra
    m.normal <- X[, !(colnames(X) %in% s.ultra)]
    m.ultra1 <- m.ultra2 <- apply(m.ultra, 2, function(x) x/2)
    colnames(m.ultra1) <- paste(colnames(m.ultra1), 1, sep='__')
    colnames(m.ultra2) <- paste(colnames(m.ultra2), 2, sep='__')
    X <- cbind(m.normal, m.ultra1, m.ultra2)
  }
  if(i>=maxiter) warning('Max. iteration reached in diluteUltraMutated')

  return(X)
}

#' Reformat spectra to column-major
#'
#' Transform spectral matrix with mutation contexts in rows to its transpose
#' with columns renamed for Mutation-signatures input
#'
#' @export
#'
row2Column <- function(X){

  nt <- c("ACAA", "ACAC", "ACAG", "ACAT", "CCAA", "CCAC", "CCAG", "CCAT",
          "GCAA", "GCAC", "GCAG", "GCAT", "TCAA", "TCAC", "TCAG", "TCAT",
          "ACGA", "ACGC", "ACGG", "ACGT", "CCGA", "CCGC", "CCGG", "CCGT",
          "GCGA", "GCGC", "GCGG", "GCGT", "TCGA", "TCGC", "TCGG", "TCGT",
          "ACTA", "ACTC", "ACTG", "ACTT", "CCTA", "CCTC", "CCTG", "CCTT",
          "GCTA", "GCTC", "GCTG", "GCTT", "TCTA", "TCTC", "TCTG", "TCTT",
          "ATAA", "ATAC", "ATAG", "ATAT", "CTAA", "CTAC", "CTAG", "CTAT",
          "GTAA", "GTAC", "GTAG", "GTAT", "TTAA", "TTAC", "TTAG", "TTAT",
          "ATCA", "ATCC", "ATCG", "ATCT", "CTCA", "CTCC", "CTCG", "CTCT",
          "GTCA", "GTCC", "GTCG", "GTCT", "TTCA", "TTCC", "TTCG", "TTCT",
          "ATGA", "ATGC", "ATGG", "ATGT", "CTGA", "CTGC", "CTGG", "CTGT",
          "GTGA", "GTGC", "GTGG", "GTGT", "TTGA", "TTGC", "TTGG", "TTGT")
  a <- rownames(X)
  b <- sapply(a, function(x){paste(strsplit(x, split = '')[[1]][c(1,3,5,7)],
                                   collapse = '')})
  idx <- match(nt, b)
  sid <- colnames(X)
  if(sum(duplicated(idx)) > 0 | sum(is.na(idx)) > 0)
    stop('Error in rownames(X)')
  X <- t(X[idx, ])
  rownames(X) <- sid
  colnames(X) <- nt
  X2 <- data.frame(Tumor_Sample_Barcode = sid, X)

  return(X2)
}

#' Generate simulated mutational catalog
#' @param W Signature matrix with rows 96 contexts and columns signatures
#' @param n Number of catalogs to generate
#' @param min.M Minimum mutation count
#' @param max.M Maximum mutation count
#' @return Catalog matrix
#' @export
simulateCatalog <- function(W, n, min.M, max.M){
  
  if(!is.matrix(W)) W <- as.matrix(W)
  if(!all(rownames(W)==trinucleotides())) stop('Rows of signat do not match trinucleotides')
  K <- NCOL(W)
  if(any(abs(colSums(W)-rep(1,K)) > 1e-5)) stop('Signature proportions not normalized')
  
  M <- exp(runif(n = n, min = log(min.M), max = log(max.M)))
  H <- matrix(runif(n = n*K), nrow = K)
  H <- t(M*t(H)/colSums(H))
  
  X <- W %*% H
  X <- apply(X, 1:2, function(x){rpois(n=1, lambda = x)})
  rownames(X) <- rownames(W)
  colnames(X) <- colnames(H) <- paste0('Sample_',seq(n))
  rownames(H) <- colnames(W)

  return(list(X=X, W=W, H=H))
}
