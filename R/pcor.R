#' Wrapper for estimating posterior probabilities with bdgraph
#' 
#' @param data An n x p matrix of data
#' @param mpl Use true for marginal pseudo-likehlihood
#' @param ... One or more arguments passed to bdgraph
#' 
#' @return A bdgraph object
#' 
#' @import BDgraph
#' 
#' @export
bdg.estimate <- function(data, mpl=FALSE, ...) {
    if (mpl) {
        return(BDgraph::bdgraph.mpl(data=data, method="ggm", ...))
    }
    else {
        return(BDgraph::bdgraph(data=data, method="ggm", ...))
    }
}

#' Estimate connected components independently
#' 
#' @param data An n x p matrix of data
#' @param prior A matrix of prior edge probabilities
#' @param mpl Use true for marginal pseudo-likehlihood
#' @param ... One or more arguments passed to bdgraph
#' 
#' @return A matrix of posterior edge probabilities
#' 
#' @import BDgraph
#' 
#' @export
bdg.islands <- function(data, prior, mpl=FALSE, ...) {
    genes <- colnames(prior)
    posterior <- matrix(0, ncol=length(genes), nrow=length(genes), dimnames=list(genes, genes))
    
    # Find connected components
    ig <- igraph::graph_from_adjacency_matrix(prior, mode="undirected", weighted=TRUE)
    connected <- igraph::components(ig)    

    # Find islands of minimum size
    islands.size <- table(connected$membership)
    islands <- names(islands.size)[islands.size >= 3]
    
    # Estimate each island independently
    for (isle in islands) {
        isle.genes <- names(connected$membership[connected$membership == isle])
        
        # Subset data and prior
        isle.data <- data[,isle.genes]
        isle.prior <- prior[isle.genes, isle.genes]
        
        # Check
        stopifnot(isle.genes == colnames(isle.data))
        stopifnot(isle.genes == colnames(isle.prior))
        
        # Estimate component
        isle.plinks <- bdg.estimate(data=isle.data, mpl=mpl, g.prior=isle.prior) %>%
                       BDgraph::plinks(round=10) %>%
                       magrittr::set_rownames(colnames(.))
        
        # Fill in posterior links
        posterior[isle.genes, isle.genes] <- isle.plinks[isle.genes, isle.genes]
    }
    
    return(posterior)
}

#' Estimate a meta-network from module eigengenes
#' 
#' @param eigs Module eigengenes
#' @param cut Threshold for selected edges from posterior probabilities
#' @param mpl Use true for marginal pseudo-likehlihood
#' @param ... One or more arguments passed to bdgraph
#' 
#' @return A list of data pertaining to resulting meta-network
#' 
#' @import BDgraph
#' @importFrom igraph V graph_from_adjacency_matrix as_edgelist
#' 
#' @export
bdg.metanet <- function(eigs, cut=0.5, mpl=FALSE, ...) {
    bdg <- bdg.estimate(data=eigs, mpl=mpl, ...)
    adj <- BDgraph::select(bdg, cut=cut)
    ig <- igraph::graph_from_adjacency_matrix(adj, mode="undirected", diag=FALSE)
    edges <- igraph::as_edgelist(ig)
    nodes <- names(igraph::V(ig))

    return(list(bdg=bdg,
                adj=adj,
                ig=ig,
                edges=edges,
                nodes=nodes))
}

#' Create a new blanket covering the entire graph search space
#' 
#' @param genes The genes of the graph
#' @param val Initial matrix values
#' 
#' @return A gene by gene matrix filled with zeroes
#' 
#' @export
blanket.new <- function(genes, val=0) {
    genes <- as.character(genes)
    blanket <- matrix(val, ncol=length(genes), nrow=length(genes), dimnames=list(genes, genes))
    return(blanket)
}

#' Lifts the blanket from within or between modules
#' 
#' @param blanket A blanket matrix
#' @param mods A list of co-expression modules to fill in
#' @param pairs A dataframe of pairs of co-expression modules to fill in
#' @param val Value for lifted links
#' 
#' @return A blanket matrix
#' 
#' @export
blanket.lift <- function(blanket, mods, pairs=NULL, val=0.5) {
    # Within modules
    for (m in mods) {
        blanket[m, m] <- val
    }
    # Between modules
    if (!is.null(pairs)) {
        for(row in seq(nrow(pairs))) {
            m1 <- mods[[pairs[row,1]]]
            m2 <- mods[[pairs[row,2]]]
            blanket[m1, m2] <- val
            blanket[m2, m1] <- val
        }
    }
    return(blanket)
}

#' Computes complexity reduction of the blanket
#' 
#' @param blanket A blanket matrix
#' 
#' @return Fold complexity reduction
#' 
#' @export
blanket.cred <- function(blanket) {
    ut <- blanket[upper.tri(blanket, diag=FALSE)]
    length(ut) / sum(ut != 0)
}

#' Adds prior information where the blanket is lifted
#' 
#' @param blanket A blanket matrix
#' @param prior.edges A matrix of prior edge probabilities
#' @param prior.blanket The blanket used in estimating prior edge probabilities
#' @param val Value for when a prior blanket is included
#' @return A blanket matrix
#' 
#' @export
blanket.inform <- function(blanket, prior.edges, prior.blanket=NULL, val=0.5) {
    blanket.informed <- blanket
    blanket.informed[blanket != 0] <- prior.edges[blanket != 0]
    if (!is.null(prior.blanket)) {
        blanket.informed[blanket != 0 & prior.blanket == 0] <- val
    }
    return(blanket.informed)
}
