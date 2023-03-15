#' Use DL-trained Model to predict HRD
#' @param catalog Matrix of mutational catalog with samples in rows and 96 trinucleotides 
#'                in columns
#' @param cancer.type Cancer type
#' @param model.path Directory path for trained models
#' @param ref.sig File path for reference signatures
#' @param threshold File path for thresholds
#' @param verbose Verbosity level
#' @param progress.bar Progress bar
#' @export

DL.call <- function(catalog, cancer.type = 'breast', model.path = NA, ref.sig = NA, threshold = NA, verbose = 1,
                    mbins = NA, progress.bar = TRUE, min.attr = 1, min.attr.pole = 10, alpha = NA){
  
  tf <- reticulate::import('tensorflow')
  pd <- reticulate::import('pandas')
  np <- reticulate::import('numpy')
  
  if(is.character(catalog))
    catalog <- read.table(catalog, header=TRUE, sep='\t')
  catalog <- as.matrix(catalog)
  if(NCOL(catalog)!=96) stop('No. of columns in catalog does not match 96 trinucleotides')
  
  cancer.type <- tolower(cancer.type) # no upper casing necessary
  known.types <- c('breast', 'ovarian', 'prostate', 'pancreatic', 'bladder','colorectal','germ_cell', 
                   'melanoma','cns','lung','head_neck','renal_cell','endometrial', 'pan_cancer')
  # default alpha
  alpha0 <- c('breast' = 0.05, 'bladder' = 0.05, 'colorectal' = 0.05, 'endometrial' = 0.05,
              'germ_cell'= 0.05, 'melanoma' = 0.05, 'head_neck' = 0.05, 'ovarian' = 0.05,
              'prostate' = 0.05, 'pancreas' = 0.05, 'cns' = 0.05, 'lung' = 0.05, 
              'renal_cell' = 0.05, 'pan_cancer' = 0.05)
  
  if(!cancer.type %in% known.types){
    if(is.na(model.path) | is.na(threshold))
    stop(paste0('Unknown cancer type ', cancer.type, 'requires model input'))
  }
  
  if(is.na(alpha))
    alpha <- alpha0[cancer.type]
  
  if(is.na(model.path))
    mfile <- system.file(paste0('extdata/dlsig/current/', cancer.type,'/models/'), package = 'DeepSig')
  else mfile <- model.path
  if(!dir.exists(mfile)) stop(paste0('Model ', mfile, ' could not be found'))
  
  sid <- rownames(catalog)
  M <- rowSums(catalog)
  nsamp <- length(M)

  if(is.na(ref.sig))
    frefsig <- system.file(paste0('extdata/dlsig/v0.9/', cancer.type, '/refsig.txt'), package = 'DeepSig')
  else
    frefsig <- ref.sig
  if(!file.exists(frefsig)) stop(paste0(frefsig, ' does not exist'))
  refsig <- read.table(frefsig, header = TRUE, sep = '\t')
  
  if(is.na(threshold))
    fthr <- system.file(paste0('extdata/dlsig/current/', cancer.type, '/mthreshold_a',alpha,'.txt'), package = 'DeepSig')
  else
    fthr <- threshold
  if(!file.exists(fthr)) stop(paste0(fthr, ' does not exist'))
  thr <- read.table(fthr, header = TRUE, sep = '\t')
  
  if(is.na(mbins))
    mbins <- system.file('extdata/dlsig/current/', cancer.type, '/mbins.txt', package = 'DeepSig')
  if(!file.exists(mbins)) stop(paste0(mbins, ' does not exist'))
  mbins <- read.table(mbins, header = TRUE, sep = '\t')

  S <- colnames(refsig)
  nsig <- length(S)
  
  pr <- matrix(0, nrow = nsamp, ncol = nsig)  # predicted score
  rownames(pr) <- sid
  colnames(pr) <- S
  bcall <- matrix(FALSE, nrow = nsamp, ncol = nsig)  # binary calls
  rownames(bcall) <- sid
  colnames(bcall) <- S

  if(verbose > 0) cat('Making binary DL calls...\n')
  xa <- np$array(catalog)
  
  if(!all(is.na(mbins)))
    mi <- vapply(M, 
      function(x){z <- which(mbins$Mmin <= x & x < mbins$Mmax); if(length(z)>1) NA else z}, integer(1))

  for(s in S){
    if(verbose > 0) cat(paste0(s,'...\n'))
#    fl <- paste0(mfile,'/', s, '/', s, '_1a')
    fl <- paste0(mfile,'/', s, '_1a')
    model <- tf$keras.models$load_model(fl)
    pr[, s] <- model$predict(xa, verbose = 0)
    if(!all(is.na(mbins))){
      tmp <- thr[thr$S==s, c('mi','thr')]
      z <- tmp$thr
      names(z) <- tmp$mi
      pthr <- z[mi]
    } else
      pthr <- thr[thr$S==s,'thr']
    if(length(pthr)==0 | !is.numeric(pthr)) 
      stop('Threshold for ', s, ' in model ', m, 'cannot be found')
    bcall[, s] <- (pr[, s] >= pthr)
  }

  if(verbose > 0) cat('Fitting catalogs to ref. sigs...\n')
  b <- DeepSig(data = t(catalog), signat = refsig)
  b <- extractSig(b, method = 'mle', progress.bar = progress.bar)
  E0 <- cbind(data.frame(sid = sid, M = M), expos(b))
  A <- E0[,-2:-1]*E0$M
  E1 <- E0
  E1[,-2:-1] <- E1[,-2:-1]*bcall*(A >= min.attr)
  pole <- c('S10','S10a','S10b','S10d','SBS10a','SBS10b','SBS10d','SBS10')
    
  # apply extra filtering to POLE signatures based on minimum attribution
  E1[,colnames(E1) %in% pole] <- E1[,colnames(E1) %in% pole]*(A[,colnames(A) %in% pole] >= min.attr.pole)
       
  B <- cbind(E1[,1:2], bcall)
  
  x <- list(exposure.raw = E0, exposure.fltrd = E1, binary.call = B, ref.sig = refsig, score = pr)
  return(x)
}
