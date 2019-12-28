#' Plot curve for soft thresholding
#' 
#' @param sft Table of values for scale-free fit
#' @param powers The power vector tested
#' @return A plot
#'  
#' @keywords internal
plt.soft <- function(sft, powers) {
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",
         type="n", 
         main=paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers, 
         cex=1, 
         col="black")
    abline(h=0.85, col="red") 
}

#' Check and choose a soft threshold
#' 
#' @param sft Table of values for scale-free fit
#' @return A theshold
#' 
#' @keywords internal
check.sft <- function(sft) {
    beta <- sft$powerEstimate
    if (is.na(beta)) {
        beta <- 6 # Default
        cat("Using the following power:", beta, "\n")
    } else {
        cat("Optimal power selected:", beta, "\n")
    }
    return(beta)
}

#' Wrapper for weighted gene co-expression analysis
#' 
#' @param eset An expression set object
#' @param min.size Minimum module size
#' @param min.sft Minimum acceptable scale-free fit when choosing soft threshold
#' @param cores Number of cpus to use
#' @param do.plot Use true to see plots
#' @return A list of data pertaining to resulting co-expression modules
#' 
#' @import WGCNA
#' @import flashClust
#' @import dynamicTreeCut
#' @import doParallel
#' @import Biobase
#' @import magrittr
#' 
#' @export
mods.get <- function(eset, min.size=10, min.sft=0.85, cores=1, do.plot=TRUE) {
    registerDoParallel(cores=cores)

    dat.expr <- t(exprs(eset))

    powers <- c(seq(1, 10, by = 1), seq(12, 20, by = 2))
    sft <- pickSoftThreshold(dat.expr, corFnc="bicor", RsquaredCut=min.sft, powerVector=powers)

    if (do.plot) plt.soft(sft, powers)

    # Check selected power
    beta <- check.sft(sft)

    adjacency <- WGCNA::adjacency(dat.expr,
                                  power=beta, 
                                  corFnc="bicor", 
                                  type="unsigned", 
                                  corOptions=list(pearsonFallback="individual"))

    dissTOM <- TOMdist(adjacency, TOMType="unsigned") 
    geneTree <- flashClust(as.dist(dissTOM), method="average")

    # All genes that are not significantly co-expressed within a module are 
    # summarized in an additional grey module (M0) for further analysis.

    # Module identification using dynamic tree cut algorithm
    modules <- cutreeDynamic(dendro=geneTree,
                             method="hybrid",
                             distM=dissTOM, 
                             deepSplit=4, 
                             pamRespectsDendro=FALSE,
                             minClusterSize=min.size)

    # Assign module colours
    modules.colored <- labels2colors(modules, zeroIsGrey=TRUE)

    # Plot the dendrogram and corresponding colour bars underneath
    if (do.plot) plotDendroAndColors(geneTree, modules.colored, "Modules", dendroLabels=FALSE, addGuide=TRUE)

    # Define modules
    genes <- rownames(eset)
    mods <- list()
    for(i in modules.colored) {
        mods[[i]] <- genes[modules.colored == i]
    }

    # Calculate eigengenes
    eigengenes <- moduleEigengenes(dat.expr, modules.colored)$eigengenes %>%
                  magrittr::set_colnames(substr(colnames(.), 3, 500))

    return(list(dat=dat.expr,
                genes=genes,
                colors=modules.colored,
                mods=mods,
                eigengenes=eigengenes))
}

#' Build a meta-network from module eigengenes
#' 
#' @param eigs Module eigengenes
#' @param cut Threshold for selected edges from posterior probabilities
#' @param mpl Use true for marginal pseudo-likehlihood
#' @param ... One or morge arguments passed to bdgraph
#' @return A list of data pertaining to resulting meta-network
#' 
#' @import BDgraph
#' @importFrom igraph V graph_from_adjacency_matrix as_edgelist
#' 
#' @export
metanet.build <- function(eigs, cut=0.5, mpl=FALSE, ...) {
    if (mpl) {
        metanet.bdg <- bdgraph.mpl(data=eigs, method="ggm", ...)
    }
    else {
        metanet.bdg <- bdgraph(data=eigs, method="ggm", ...)
    }
    
    metanet.adj <- BDgraph::select(metanet.bdg, cut=cut)
    metanet.ig  <- igraph::graph_from_adjacency_matrix(metanet.adj, mode="undirected", diag=FALSE)
    metanet.edges <- igraph::as_edgelist(metanet.ig)
    metanet.nodes <- names(igraph::V(metanet.ig))

    return(list(metanet.bdg   = metanet.bdg,
                metanet.adj   = metanet.adj,
                metanet.ig    = metanet.ig,
                metanet.edges = metanet.edges,
                metanet.nodes = metanet.nodes))
}

#' Create a new blanket covering the entire graph
#' 
#' @param genes The genes of the graph
#' @return A gene by gene matrix filled with zeroes
#' 
#' @export
blanket.new <- function(genes) {
    genes <- as.character(genes)
    blanket <- matrix(0, ncol=length(genes), nrow=length(genes), dimnames=list(genes, genes))
    return(blanket)
}

#' Lifts the blanket from single modules
#' 
#' @param blanket A blanket matrix
#' @param mods A list of co-expression modules to fill in
#' @return A blanket matrix
#' 
#' @export
blanket.add.mods <- function(blanket, mods) {
    for (m in mods) {
        blanket[m, m] <- 0.5
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
blanket.add.modpairs <- function(blanket, mods, modpairs) {
    for(row in seq(nrow(modpairs))) {
        m1 <- mods[[modpairs[row,1]]]
        m2 <- mods[[modpairs[row,2]]]
        blanket[m1, m2] <- 0.5
        blanket[m2, m1] <- 0.5
    }
    return(blanket)
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
