#' Generate Cosine Similarity Matrix Between De novo and Reference Signatures
#'
#' @param sig Input de novo signature matrix
#' @param ref.sig Reference signature set
#' @return Cosine similarity matrix
#' @examples
#' set.seed(135)
#' @export
denovoCosim <- function(sig, ref.sig = 'cosmic'){
  
  if(is.character(ref.sig)){
    if(ref.sig=='cosmic'){
       ref.sig <- read.table(
         system.file('extdata','cosmic_sigProfiler_SBS_signatures_v3.1.txt', 
                  package = 'DeepSig'), header = TRUE, sep = '\t')
    }else{
      if(!file.exists(ref.sig)) stop(paste0(ref.sig,' does not exist'))
      ref.sig <- read.table(ref.sig, header = TRUE, sep = '\t')
    }
  }
  if(!class(ref.sig) %in% c('data.frame','matrix')) 
    stop('ref.sig not of correct type')
  
  if(NROW(sig)!=96) stop('Inut sig does not have 96 rows')
  
  x <- cosineSimilarity(A = sig, B = ref.sig, diag = FALSE)
  return(x)
}

#' Cosine Similarity Matrix Heatmap
#'
#' @param cs Cosine similarity matrix
#' @param show.lsap Display LSAP assignment as rectangles
#' @return NULL
#' @examples
#' set.seed(135)
#' @export
csHeatmap <- function(z, pal=NA, zmin=0, zmax=1, cex=1, ddx=0.03, ddy=0.03, 
                      mar=c(1,1,5,2), grid.col='gray', na.col='gray',legend=NA, 
                      rescale=FALSE,lwd=0.5, col.break = 30, show.lsap = TRUE){
  
    if(any(is.na(pal))) pal <- colorRampPalette(RColorBrewer::brewer.pal(9,'OrRd'))(10)
    nstep <- length(pal)
    N <- nrow(z)
    K <- ncol(z)
    nK <- ceiling(K/col.break)
    Kmax <- min(col.break, K)
    par(mfrow=c(nK, 1), mar=mar)
    dx <- 1/Kmax
    dy <- 1/N
    xgrid <- seq(0,1,by=dx)
    ygrid <- seq(0,1,by=dy)
    if(rescale) zs <- z /(zmax - zmin)
    else zs <- z
    zs <- apply(zs, 1:2, function(x){if(is.na(x)) NA else if(x< zmin) zmin else x})
    zs <- apply(zs, 1:2, function(x){if(is.na(x)) NA else if(x> zmax) zmax else x})
    dz <- (zmax - zmin)/(nstep-1)
    bz <- apply(zs,1:2, function(x){as.integer((x-zmin)/dz)+1})
    bsize.x <- dx*(nstep-bz+2)*ddx
    bsize.y <- dy*(nstep-bz+2)*ddy
    
    lsap <- sigLSAP(z)
    ref.sig <- colnames(z)
      
    for(k in seq(nK)){
      plot(NA, xlim=c(0,1), ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
      for(i in seq(Kmax)){
        i1 <- i+(k-1)*Kmax
        if(i1 > K) next()
        for(j in seq(N)){
          if(is.na(bz[N-j+1,i1])){
            if(is.na(na.col)) next()
            ddx <- ddy <- 0
            col <- na.col
          }else{
            ddx <- bsize.x[N-j+1,i1]
            ddy <- bsize.y[N-j+1,i1]
            col <- pal[bz[N-j+1,i1]]
          }
          if(show.lsap & ref.sig[i1] == lsap[N-j+1]){
            lty <- 2
            bwd <- 2*lwd
            bcol <- 'blue'
          } else{
            lty <- 1
            bwd <- 1*lwd
            bcol <- grid.col
          }
          rect(xleft = xgrid[i], xright = xgrid[i+1], ybottom=ygrid[j], 
               ytop=ygrid[j+1], col='white', border=bcol,
               lwd=bwd, lty=lty)
          rect(xleft = xgrid[i] + ddx, xright = xgrid[i+1]- ddx,
            ybottom=ygrid[j] + ddy, ytop=ygrid[j+1] - ddy, col=col, 
            border=NA,lwd=lwd)
        }
      }
      text(x=-dx/5,y=(seq(N)-0.5)*dy, label=rev(rownames(z)),cex=cex, adj=1,xpd=NA)
      idx <- seq((k-1)*Kmax+1, min(k*Kmax, K))
      
      text(x=((idx-(k-1)*Kmax)-0.5)*dx,y=1+dy/3,label=colnames(z)[idx], cex=cex, adj=0,srt=45,xpd=NA)
    
      if(k==1){
        legend.height=0.3
        ldy <- legend.height/nstep
        ytop <- seq(1,1-legend.height, length.out=nstep)
        ybottom <- ytop-ldy
        rect(xleft=max(xgrid)+dx*0.3,xright=max(xgrid)+dx*0.7, ytop=ytop,
             ybottom=seq(1-ldy,1-legend.height-ldy, length.out=nstep), col=rev(pal), 
             xpd=NA,lwd=lwd/2)
        if(any(is.na(legend))){
          z0 <- round(zmin,digits=2)
          if(z0==0) z0 <- '0.0'
          z0.5 <- round((zmax-zmin)/2,digits=2)
          z1 <- round(zmax,digits=2)
          if(z1==1) z1 <- '1.0'
        } else{
          z0 <- legend[1]
          z0.5 <- legend[2]
          z1 <- legend[3]
        }
        text(x=max(xgrid)+dx*0.8, y=min(ybottom)+legend.height/nstep*0.5, 
             label=z0, cex=cex, xpd=NA,adj=0)
        text(x=max(xgrid)+dx*0.8, y=mean(c(max(ybottom),min(ybottom))) +
               legend.height/nstep*0.5, label=z0.5, cex=cex, xpd=NA,adj=0)
        text(x=max(xgrid)+dx*0.8, y=max(ybottom)+legend.height/nstep*0.5, 
             label=z1, cex=cex, xpd=NA,adj=0)
      }
    }
}

#' Automatic Signature Annotation
#' 
#' @param x Signature cosine similarity matrix (rows: de novo, columns = reference)
#' @return Vector of assigned reference signature names for de novo signatures
#' @export

sigLSAP <- function(x){
  
  S <- rownames(x)
  Sref <- colnames(x)
  z <- clue::solve_LSAP(x, maximum = TRUE)
  z <- Sref[as.integer(z)]
  names(z) <- S
  return(z)
}
#' Plot De Novo Signatures and Closest Reference Signatures
#' 
#' @param object Object whose signature set is to be plotted
#' @param ref.sig Reference signature set
#' @export
denovoProfiles <- function(object, ref.sig = 'cosmic'){
  
  if(ref.sig=='cosmic'){
    ref.sig <- read.table(
      system.file('extdata','cosmic_sigProfiler_SBS_signatures_v3.1.txt', 
                  package = 'DeepSig'), header = TRUE, sep = '\t')
  }else if(is.character(ref.sig)){
    if(!file.exists(ref.sig)) stop(paste0(ref.sig,' does not exist'))
    ref.sig <- read.table(ref.sig, header = TRUE, sep = '\t')
  }
  
  sig <- signat(object)
  cosim <- denovoCosim(sig, ref.sig)
  
  K <- NCOL(sig)
  old.par <- par(mfrow=c(K,3), mar=c(1.5,3,2,1),lwd = 0.1, cex.axis = 0.6, 
      cex.lab = 0.6, tck=-0.05, mgp = c(2,0.5,0), cex.main = 1)
  S <- colnames(sig)
  
  for(k in seq(K)){
    sigplot(t(sig)[k, ])
    x <- cosim[S[k], ]
    x <- x[order(x, decreasing = TRUE)]
    c1 <- names(x)[1]
    c2 <- names(x)[2]
    sk <- S[k]
    title(main = paste0(sk))
    sigplot(t(ref.sig)[c1, ])
    title(main = paste0(c1, ' (', round(x[c1], digits = 2), ')'))
    sigplot(t(ref.sig)[c2, ])
    title(main = paste0(c2, ' (', round(x[c2], digits = 2), ')'))
  }
  par(old.par)
  return(invisible(object))
}