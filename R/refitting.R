#' Cosine similarity
#' 
#' Row names of two data matrices (mutation contexts) are used for comparison.
#' If the row names do not match completely, the larger set is used as reference 
#' with the smaller set expanded with zero-padding.
#' 
#' @param A Test matrix of dimension \code{(m,a)} 
#' @param B Reference matrix of dimension \code{(m, b)}
#' @param diag Only compare column \code{A[, i]} with column \code{B[, i]} 
#'        where \code{i=1, ..., ncol(A)=ncol(B)}.
#' @return If \code{diag = TRUE}, vector of overlap between columns of \code{A}
#'         and columns of \code{B} in one-to-one mapping; if \code{diag = FALSE},
#'         matrix of dimension \code{(a,b)}, whose elements give overlap of column 
#'         \code{a} in matrix \code{A} with column \code{b} in matrix \code{B}.
#' @export
cosineSimilarity <- function(A, B, diag = FALSE){
  
  if(any(is.na(rownames(A)))) stop('rownames(A) not complete')
  if(any(is.na(rownames(B)))) stop('rownames(B) not complete')
  if(!is.matrix(A)) A <- as.matrix(A)
  if(!is.matrix(B)) B <- as.matrix(B)
  a <- NCOL(A)
  b <- NCOL(B)
  if(diag & a != b) stop('diag requires same dimension of A and B')
  
  if(diag){ 
    cos <- rep(0, a)
    names(cos) <- colnames(A)
  } else{ 
    cos <- matrix(0, nrow=a, ncol=b)
    rownames(cos) <- colnames(A)
    colnames(cos) <- colnames(B)
  }
    
  for(i in seq(1,a)) for(j in seq(1,b)){
    if(diag & i!=j) next()
    xa <- A[,i]
    xb <- B[,j]
    if(length(xa) < length(xb)){
      xtest <- xa
      xref <- xb
    } else{
      xtest <- xb
      xref <- xa
    }
    xtest2 <- rep(0, length(xref))
    names(xtest2) <- names(xref)
    if(sum(is.na(names(xtest)) > 0) | 
       sum(!names(xtest) %in% names(xref)) > 0) stop('Names mismatch')
    xtest2[names(xtest)] <- xtest
    if(sum(xtest2)==0 | sum(xref)==0) x <- 0
    else x <- sum(xtest2*xref) / sqrt(sum(xtest2^2) * sum(xref^2))
    if(diag) cos[i] <- x
    else cos[i,j] <- x
  }
  return(cos)
}
