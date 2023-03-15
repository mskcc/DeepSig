#' Generate Newick format tree string from tree list object 
#' 
#' @param tree Tree list object from \code{\link{build_tree}}
#' @param parent Parent ID 
#' @param string Newick string of parent tree
#' @return String of newick tree
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50,ncol=100,rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s,ranks=seq(2,8),nrun=5)
#' tree <- build_tree(s,rmax=5)
#' nw <- newick(tree=tree)
#' nw
#' @export
newick <- function(tree, parent='1.1',string=''){
  
  if(string=='') root=TRUE
  else root=FALSE
  string <- paste0(string,'(')
  for(i in seq_len(length(tree))){
    if(length(tree[[i]])==1){
      bl <- as.numeric(strsplit(tree[[i]],split='[.]')[[1]][1]) -
        as.numeric(strsplit(parent,split='[.]')[[1]][1])
      string <- paste0(string,tree[[i]],':',toString(bl))
    }
    else{
      sub <- newick(tree[[i]], parent=names(tree)[i], string)
      bl <- as.numeric(strsplit(names(tree)[i],split='[.]')[[1]][1]) -
        as.numeric(strsplit(parent,split='[.]')[[1]][1])
      string <- paste0(sub, names(tree)[i],':',toString(bl))
    }
    if(i < length(tree)) string <- paste0(string,',')
    else string <-paste0(string,')')
  }
  if(root) string <- paste0(string,';')
  return(string)
}

# Recursive function for branching a tree
branch.tree <- function(tree, parent.id, progenies){
  for(i in seq_len(length(tree))){
    if(length(tree[[i]])==1){
      if(tree[[i]] == parent.id){
        tree[[i]] <- progenies
        names(tree)[i] <- parent.id
        return(tree)
      }
    } else{
      tree[[i]] <- branch.tree(tree[[i]], parent.id, progenies)
    }
  }
  return(tree)
}

update.tree <- function(tree, parent.id, progenies){
  for(i in seq_len(length(tree))){
    if(length(tree[[i]])==1){
      if(tree[[i]] %in% parent.id)
       tree[[i]] <- progenies[tree[[i]]==parent.id]
    } else
       tree[[i]] <- update.tree(tree[[i]], parent.id, progenies)
  }
  return(tree)
}

#' Build tree connecting clusters at different ranks
#' 
#' @param object Object of class \code{scNMFSet}
#' @param rmax Maximum rank at which tree branching stops
#' @return List containing the tree structure
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50,ncol=100,rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s,ranks=seq(2,8),nrun=5)
#' tree <- build_tree(s,rmax=5)
#' tree
#' @export
build_tree <- function(object, rmax){
  
  ranks <- as.integer(gsub('K','',names(object$W)))
  r0 <- 3
  if(missing(rmax)) rmax <- ranks[length(ranks)]
#  id <- seq(which(ranks == r0-1),
#            which(ranks == rmax))
#  cluster <- NULL
#  for(i in id){
#    h <- object$H[[i]]
#    cid <- apply(h,2,which.max)
#    cluster <- cbind(cluster,cid)
#    colnames(cluster)[ncol(cluster)] <- paste0('r', ranks[i])
#  }
  
  rse <- seq(r0,rmax,1)
  y <- vector('list',length(rse))
  tree <- list('2.1','2.2') # initial split for rank 2
  names(y) <- rse
  for(k in rse){
    rank0 <- paste0('K',toString(k-1))
    rank1 <- paste0('K',toString(k))
    x <- cosineSimilarity(A=object$W[[rank0]], B=object$W[[rank1]])
    z <- apply(x, 2, which.max)
    z <- vapply(z, function(x){if(length(x)>1) x[1] else x}, integer(1))
        # break ties
    names(z) <- seq(k)
    y[[toString(k)]] <- z
    
    w <- names(table(z))[table(z) > 1]  # parent id that splitted
    p <- mapply(function(x){names(z)[z==x]},w, SIMPLIFY=FALSE)
    
    for(i in seq_len(length(w))){
      np <- paste(k,p[[i]],sep='.')
      split <- as.list(np)
      tree <- branch.tree(tree, paste(k-1,w[i],sep='.'), split)
    }
    
    w0 <- names(table(z))[table(z)==1]
    p0 <- names(z)[match(w0,z)]
    rename <- paste(k,p0,sep='.')
    tree <- update.tree(tree, paste(k-1,w0,sep='.'), rename)
  }
  return(tree)
}

#' Rename tips of trees with cell types
#' 
#' @param tree List containing tree
#' @param tip.labels Vector of new names for tips
#' @param rank Rank value of which tip names are to be replaced
#' @return List containing tree with updated tip labels 
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50,ncol=100,rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s,ranks=seq(2,8),nrun=5)
#' tree <- build_tree(s,rmax=5)
#' tree <- rename_tips(tree,rank=5,tip.labels=letters[seq_len(5)])
#' tree 
#' @export
rename_tips <- function(tree, rank, tip.labels){
  
  for(i in seq_len(length(tree))){
    if(length(tree[[i]])==1){  # a tip
      x <- strsplit(tree[[i]], split='[.]')[[1]]
      if(as.numeric(x[1])==rank){
        cid <- as.numeric(x[2])
        tree[[i]] <- paste(toString(rank),tip.labels[cid],sep='.')
      }
    }else tree[[i]] <- rename_tips(tree[[i]], rank, tip.labels)
  }
  return(tree)
}

#' Plot cluster tree
#' 
#' Visualize the output of \code{\link{build_tree}} as a dendrogram.
#' 
#' Uses \code{\link[ape]{plot.phylo}} to visualize cluster tree.
#' 
#' @param tree List containing tree structure. Output from 
#'   \code{\link{build_tree}}
#' @param direction \code{c('rightwards','downwards')}; 
#'   the direction of dendrogram
#' @param cex Font size of edge/tip labels
#' @param ... Other parameters to \code{\link{plot.phylo}}
#' @return \code{NULL}
#' @examples
#' set.seed(1)
#' x <- simulate_whx(nrow=50,ncol=100,rank=5)
#' s <- scNMFSet(x$x)
#' s <- vb_factorize(s,ranks=seq(2,8),nrun=5)
#' tree <- build_tree(s,rmax=5)
#' plot_tree(tree)
#' @export
plot_tree <- function(tree, direction='rightwards', cex = 0.7, ...){
  
  nwk <- newick(tree)
  phylo <- ape::read.tree(text=nwk)
  plot(x=phylo, direction=direction, cex=cex, ...)
     
  return(invisible(tree))
}
