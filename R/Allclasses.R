#' Class \code{DeepSig} for storing catalog and exposure
#'
#' \code{S4} class with storage for input catalog, reference signatures,
#' and output exposure inferred.
#'
#' @slot catalog Matrix of mutation count with samples in columns
#'               and mutation types in rows
#' @slot signat Matrix of reference signatures; each column
#'                 stores individual signatures proportions with mutation
#'                 types in rows
#' @slot tmb Vector of total mutation loads for each sample
#' @slot expos Matrix of inferred proportions of signatures
#'                exposure; each row gives proportions of each
#'                sample with signatures in columns
#' @slot pvalue P-values of exposures estimated from permutation tests;
#'            same dimension as \code{expos}
#' @slot logLik Log likelihood
#' @slot misc List of factorization results and measures for de novo runs
#' @return Object of class \code{DeepSig}
#' @export DeepSig
setClass('DeepSig',
         slots = c(catalog = 'matrix',
                   signat = 'matrix',
                   tmb = 'vector',
                   expos = 'matrix',
                   pvalue = 'matrix',
                   logLik = 'vector',
                   misc = 'list'
         ))

#' Create \code{DeepSig} object
#'
#' @param data Input data of mutation catalog; samples in columns
#'             and mutation types in rows.
#' @param signat Reference signatures table; default is COSMIC
#'                  version 3 signatures. It can be a matrix of the table, the file name
#'                  containing the table, or character code \code{c('v2', 'SA','SP')},
#'                  signifying \code{COSMIC v2}, \code{sigAnalyzer}, or
#'                  \code{SigProfiler} references (the latter two are COSMIC v3).
#' @return Object of class \code{DeepSig}.
#' @examples
#' set.seed(130)
#' data <- rmultinom(n = 5, size = 100, prob = rep(1,96))
#' rownames(data) <- trinucleotides()
#' colnames(data) <- paste0('S',seq(5))
#' s <- DeepSig(data)
#' data <- read.table(system.file('extdata', 'tcga-brca_catalog.txt',
#'                   package='DeepSig'), header=TRUE, sep='\t')
#' b <- DeepSig(data)
#' @export
DeepSig <- function(data, signat = 'v2'){

  if(!is(data, 'matrix')) data <- as.matrix(data)
  if(is(signat, 'data.frame')) signat <- as.matrix(signat)
  else if(!is(signat, 'matrix')){
    if(!is(signat, 'character')) stop('Unusuable input of signat')
    if(is.null(signat) | signat == 'v2')
      signat <- as.matrix(read.table(system.file('extdata/cosmic_snv_signatures_v2.txt', 
                                                 package = 'DeepSig')))
    else if(signat == 'SA')
      signat <- as.matrix(read.table(system.file('extdata/cosmic_SigAnalyzer_SBS_signatures.txt', 
                    package = 'DeepSig')))
    else if(signat == 'SP')
      signat <- as.matrix(read.table(system.file('extdata/cosmic_sigProfiler_SBS_signatures.txt', 
                    package = 'DeepSig')))
    else{
      if(!file.exists(signat)) stop(paste0('File ',signat,' does not exist'))
      signat <- as.matrix(read.table(signat))
    }
  }
  
  signat <- t(t(signat) / colSums(signat)) # renormalize just in case (also to remove truncation errors)
  if(!all(abs(colSums(signat) - 1) < 1e-4)) stop('Signature list not normalized')
  nts <- rownames(signat)
  ntd <- rownames(data)
  if(sum(is.na(ntd)) > 0) stop('Data must have explicit mutation type names')
  idx <- match(nts, ntd)
  if(sum(is.na(idx)) > 0 | sum(duplicated(idx)) > 0)
    stop('Row names in data do not match reference signatures')
  data <- data[idx, , drop = FALSE]
  
  x <- new('DeepSig', catalog = data, signat = signat)
  x@tmb <- colSums(data)

  return(x)
}

#' Display \code{DeepSig} object
#'
#' Display the class and dimension
#'
#' @param object Object of class \code{DeepSig}
#' @return \code{NULL}
#' @export
setMethod('show', signature = 'DeepSig',
          definition = function(object){
            cat('Class:', class(object), '\n')
            ca <- catalog(object)
            cat('Number of mutation types = ', NROW(ca), '\n')
            cat('  ')
            print(utils::head(rownames(ca)))
            cat('Number of samples = ', NCOL(ca), '\n')
            sg <- signat(object)
            cat('Number of signatures = ', NCOL(sg), '\n')
            cat('  ')
            print(utils::head(colnames(sg)))
          })

#' @export
setGeneric('catalog', function(object) standardGeneric('catalog'))
#' Accessor for catalog
#'
#' @param object Object containing catalog
#' @return \code{catalog} data frame
#' @export
setMethod('catalog', signature = 'DeepSig',
          function(object){
            object@catalog
          }
)
#' @export
setGeneric('catalog<-', function(object, value) standardGeneric('catalog<-'))
#' @export
setMethod('catalog<-', signature = 'DeepSig',
          function(object, value){
            object@catalog <- as.matrix(value)
            if(validObject(object)) return(object)
          })
#' @export
setGeneric('signat', function(object) standardGeneric('signat'))
#' Accessor for signature
#'
#' @param object Object containing signatures
#' @return \code{signat} data frame
#' @export
setMethod('signat', signature = 'DeepSig',
          function(object){
            object@signat
          }
)
#' @export
setGeneric('signat<-', function(object, value) standardGeneric('signat<-'))
#' @export
setMethod('signat<-', signature = 'DeepSig',
          function(object, value){
            object@signat <- as.matrix(value)
            if(validObject(object)) return(object)
          })
#' @export
setGeneric('tmb', function(object) standardGeneric('tmb'))
#' Accessor for TMB
#'
#' @param object Object containing \code{tmb}
#' @return \code{tmb} data frame
#' @export
setMethod('tmb', signature = 'DeepSig',
          function(object){
            object@tmb
          }
)
#' @export
setGeneric('tmb<-', function(object, value) standardGeneric('tmb<-'))
#' @export
setMethod('tmb<-', signature = 'DeepSig',
          function(object, value){
            object@tmb <- value
            if(validObject(object)) return(object)
          })
#' @export
setGeneric('expos', function(object) standardGeneric('expos'))
#' Accessor for exposure
#'
#' @param object Object containing \code{expo}
#' @return \code{expos} data frame
#' @export
setMethod('expos', signature = 'DeepSig',
          function(object){
            object@expos
          }
)
#' @export
setGeneric('expos<-', function(object, value) standardGeneric('expos<-'))
#' @export
setMethod('expos<-', signature = 'DeepSig',
          function(object, value){
            object@expos <- as.matrix(value)
            if(validObject(object)) return(object)
})
#' @export
setGeneric('pvalue', function(object) standardGeneric('pvalue'))
#' Accessor for p-values
#'
#' @param object Object containing \code{pvalue}
#' @return \code{pvalue} data frame
#' @export
setMethod('pvalue', signature = 'DeepSig',
          function(object){
            object@pvalue
          }
)
#' @export
setGeneric('pvalue<-', function(object, value) standardGeneric('pvalue<-'))
#' @export
setMethod('pvalue<-', signature = 'DeepSig',
          function(object, value){
            object@pvalue <- as.matrix(value)
            if(validObject(object)) return(object)
          })
#' Accessor for logLik
#'
#' @param object Object containing \code{logLik}
#' @return \code{logLik} value
#' @export
setMethod('logLik', signature = 'DeepSig',
          function(object){
            object@logLik
          }
)
#' @export
setGeneric('logLik<-', function(object, value) standardGeneric('logLik<-'))
#' @export
setMethod('logLik<-', signature = 'DeepSig',
          function(object, value){
            object@logLik <- value
            if(validObject(object)) return(object)
          })
#' @export
setGeneric('misc', function(object) standardGeneric('misc'))
#' Accessor for misc
#'
#' @param object Object containing \code{misc}
#' @return List \code{misc}
#' @export
setMethod('misc', signature = 'DeepSig',
          function(object){
            object@misc
          }
)
#' @export
setGeneric('misc<-', function(object, value) standardGeneric('misc<-'))
#' @export
setMethod('misc<-', signature = 'DeepSig',
          function(object, value){
            object@misc <- value
            if(validObject(object)) return(object)
          })
