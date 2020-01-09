#' Wrapper for estimating posterior probabilities with bdgraph
#' 
#' @param data An n x p matrix of data
#' @param mpl Use true for marginal pseudo-likehlihood
#' @param ... One or more arguments passed to bdgraph
#' @return A bdgraph object
#' 
#' @import BDgraph
#' 
#' @export
bdg.estimate <- function(data, mpl=FALSE, ...) {
    if (mpl) {
        return( BDgraph::bdgraph.mpl(data=data, method="ggm", ...) )
    }
    else {
        return( BDgraph::bdgraph(data=data, method="ggm", ...) )
    }
}

#' Estimate a meta-network from module eigengenes
#' 
#' @param eigs Module eigengenes
#' @param cut Threshold for selected edges from posterior probabilities
#' @param mpl Use true for marginal pseudo-likehlihood
#' @param ... One or more arguments passed to bdgraph
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

#' Estimate connected components independently
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

#' Create a new blanket covering the entire graph
#' 
#' @param genes The genes of the graph
#' @return A gene by gene matrix filled with zeroes
#' 
#' @export
blanket.new <- function(genes, val=0) {
    genes <- as.character(genes)
    blanket <- matrix(val, ncol=length(genes), nrow=length(genes), dimnames=list(genes, genes))
    return(blanket)
}

# A good test
#blanket <- blanket.new(wgcna$genes, val=0)
#blanket <- blanket.add.mods(blanket, wgcna$mods, val=1)
#table(blanket)
#c <- components(igraph::graph_from_adjacency_matrix(blanket, mode="undirected"))
#sort(table(wgcna$colors))
#sort(table(c$membership))

#' Lifts the blanket from single modules
#' 
#' @param blanket A blanket matrix
#' @param mods A list of co-expression modules to fill in
#' @return A blanket matrix
#' 
#' @export
blanket.add.mods <- function(blanket, mods, val=0.5) {
    for (m in mods) {
        blanket[m, m] <- val
    }
    return(blanket)
}

#' Lifts the blanket from pairs of modules
#' 
#' @param blanket A blanket matrix
#' @param mods A list of co-expression modules
#' @param modpairs A dataframe of pairs of co-expression modules to fill in
#' @return A blanket matrix
#' 
#' @export
blanket.add.modpairs <- function(blanket, mods, modpairs, val=0.5) {
    for(row in seq(nrow(modpairs))) {
        m1 <- mods[[modpairs[row,1]]]
        m2 <- mods[[modpairs[row,2]]]
        blanket[m1, m2] <- val
        blanket[m2, m1] <- val
    }
    return(blanket)
}

blanket.cred <- function(blanket) {
    ut <- blanket[upper.tri(blanket, diag=FALSE)]
    length(ut) / sum(ut != 0)
}

#' Adds prior information where the blanket is lifted
#' 
#' @param blanket A blanket matrix
#' @param prior A matrix of probabilities with the same dimensions as the blanket
#' @return A blanket matrix
#' 
#' @export
blanket.inform <- function(blanket, prior) {
    blanket.informed <- blanket
    blanket.where.not.zero <- blanket.informed != 0
    blanket.informed[blanket.where.not.zero] <- prior[blanket.where.not.zero]
    return(blanket.informed)
}
