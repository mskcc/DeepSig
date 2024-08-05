#' Use DL-trained Model to Make Signature Calls
#' @param catalog Matrix of mutational catalog with samples in rows and 96 trinucleotides 
#'                in columns
#' @param cancer.type Cancer type
#' @param model.path Directory path for trained models, If `NA`, defaults to
#'          model corresponding to `cancer.type` and will query github API to download
#' @param min.M Minimum no. of mutations
#' @param min.attr Mininum attribution
#' @param verbose Verbosity level
#' @param progress.bar Progress bar
#' @param ... Other parameters to `modelFetch`
#' @examples 
#' data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt',package='DeepSig'))
#' z <- DL.call(catalog = t(data), cancer.type = 'breast')
#' head(z$exposure.fltrd)
#' 
#' @export

DL.call <- function(catalog, cancer.type = 'pancancer', model.path = './.DeepSig', 
                    mode = 'catalog', verbose = 1,  progress.bar = TRUE, 
                    min.M = 1,  min.attr = 1, ...){
  
  tf <- reticulate::import('tensorflow')
  pd <- reticulate::import('pandas')
  np <- reticulate::import('numpy')
  
  if(is.character(catalog))
    catalog <- read.table(catalog, header=TRUE, sep='\t', check.names = F)
  catalog <- as.matrix(catalog)
  if(mode=='catalog')
    if(NCOL(catalog)!=96) stop('No. of columns in catalog does not match 96 trinucleotides')
  
  cancer.type <- tolower(cancer.type) # no upper casing necessary
  known.types <- c('breast', 'ovary', 'prostate', 'pancreas', 'bladder','colorect',
                   'gct','skin','cns','head_neck','kidney','uterus', 'nsclc','liver',
                   'sclc', 'pancancer')
  
  if(!cancer.type %in% known.types){  # Convert Oncotree code into cancer.type
    cancer.type <- oncoTree(onco = cancer.type)
  }
#  if(!cancer.type %in% known.types){
#    if(is.na(model.path) | is.na(threshold))
#    stop(paste0('Unknown cancer type ', cancer.type, 'requires model input'))
#  }
  
  mfile <- paste(model.path, cancer.type, sep = '/')
  if(!dir.exists(mfile)){
     if(verbose) cat('Model for ', cancer.type, ' cannot be found. Downloading...\n')
     if(!dir.exists(model.path)) dir.create(model.path)
     if(!dir.exists(mfile)) dir.create(mfile)
     modelFetch(cancer.type = cancer.type, verbose = verbose, 
                         model.path = mfile, ...)
  }
  
  x0 <- catalog
  sid0 <- rownames(catalog)
  M0 <- rowSums(catalog)

  sid.bad <- sid0[M0 < min.M]
  sid <- sid0[M0 >= min.M]
  M <- M0[M0 >= min.M]
  nsamp <- length(M)
  catalog <- catalog[sid,,drop=F]

  frefsig <- paste(mfile,'refsig.txt',sep='/')
  if(!file.exists(frefsig)) stop(paste0(frefsig, ' does not exist'))
  refsig <- read.table(frefsig, header = TRUE, sep = '\t')
      
  fthr <- paste(mfile, 'threshold_cut.txt', sep='/')
  if(!file.exists(fthr)) stop(paste0(fthr, ' does not exist'))
  thr <- read.table(fthr, header = TRUE, sep = '\t')
  
  S <- colnames(refsig)
  nsig <- length(S)
  
  tcall <- matrix('N', nrow = nsamp, ncol = nsig)  # ternary calls
  rownames(tcall) <- sid
  colnames(tcall) <- S
  
  if(verbose > 0) cat('Making ternary DL calls...\n')
  if(mode=='catalog')
    xa <- np$array(catalog)
  else{
    z <- expos(b)
    z <- z*tmb(b)
    xa <- np$array(z)
  }
  pr <- matrix(0, nrow = nsamp, ncol = nsig)  # predicted score
  rownames(pr) <- sid
  colnames(pr) <- S
  
  for(s in S){
    if(verbose > 0) cat(paste0(s,'...\n'))
    fl <- paste0(mfile,'/',s)
    if(!dir.exists(fl))
    fl <- paste0(mfile, '/',s, '/', s)
    if(!dir.exists(fl)) stop(paste0('Trained model for ',s,' cannot be found'))
    model <- tf$keras.models$load_model(fl)
    mp <- model$predict(xa, verbose = 0)
    pr[, s] <- mp
    tmp <- thr[thr$S==s, c('Threshold', 'Precision.cutoff')]
    if(NROW(tmp)!=2) stop('Less/more than 2 cutoff exists in threshold file')
    z1 <- tmp[tmp$Precision.cutoff == min(tmp$Precision.cutoff),'Threshold']
    z2 <- tmp[tmp$Precision.cutoff == max(tmp$Precision.cutoff),'Threshold']
    tcall[pr[, s] > z1, s] <- 'I'
    tcall[pr[, s] > z2, s] <- 'P'
  }

  if(verbose > 0) cat('Fitting catalogs to ref. sigs...\n')
  b <- DeepSig(data = t(catalog), signat = refsig)
  screen <- as.matrix(tcall!='N')
  b <- extractSig(b, method = 'mle', screen = screen, min.tmb = min.M, progress.bar = progress.bar)
  E0 <- cbind(data.frame(sid = sid, M = M), expos(b))

  A <- E0[,-2:-1]*E0$M
  E1 <- E0
  E1[,-2:-1] <- E1[,-2:-1]*(A >= min.attr)

  B <- cbind(E1[,1:2], tcall) # ternary calls
  pr <- cbind(E1[,1:2], pr)   # scores
  
  if(length(sid.bad)>0){
    na <- matrix(NA, nrow=length(sid.bad), ncol=nsig)
    rownames(na) <- sid.bad
    colnames(na) <- S
    E1 <- rbind(E1, cbind(data.frame(sid = sid.bad, M = M0[sid.bad]), na))
    E1 <- E1[sid0, ]
    B <- rbind(B, cbind(data.frame(sid=sid.bad, M = M0[sid.bad], na)))
    B <- B[sid0, ]
    pr <- rbind(pr, cbind(data.frame(sid=sid.bad, M = M0[sid.bad], na)))
    pr <- pr[sid0, ]
  }
  
  x <- list(score = pr,
            ternary.call = B, 
            exposure = E1, 
            ref.sig = refsig, 
            threshold=thr)
  return(x)
}
