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
#' @param cor.fn Method for calculation co-expression similarity
#' @param powers A vector of values to test for soft thresholding
#' @param merging Merge similar modules by eigengene
#' @param merging.cut Maximum dissimilarity that qualifies modules for merging
#' @param do.plot Use true to see plots
#' @return A list of data pertaining to resulting co-expression modules
#' 
#' @import WGCNA
#' @import flashClust
#' @import dynamicTreeCut
#' @import doParallel
#' @import Biobase
#' 
#' @export
mods.detect <- function(eset,
                        min.size=10, 
                        min.sft=0.85, 
                        cores=1,
                        cor.fn=c("bicor", "cor"),
                        powers=c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                        merging=FALSE,
                        merging.cut=0.2,
                        do.plot=TRUE) {
    
    # Handle arguments
    args <- as.list(environment())
    cor.fn <- match.arg(cor.fn)
    
    # Correlation options
    if (cor.fn == "cor") cor.options = list(use="p")
    if (cor.fn == "bicor") cor.options = list(pearsonFallback="individual")
    
    # Set parallel computing environment
    doParallel::registerDoParallel(cores=cores)
    
    # Format expression set into sample x gene matrix
    dat <- t(Biobase::exprs(eset))
    
    # Pick soft threshold via scale-free fit
    sft <- WGCNA::pickSoftThreshold(data=dat, 
                                    corFnc=cor.fn, 
                                    RsquaredCut=min.sft, 
                                    powerVector=powers)
    
    if (do.plot) plt.soft(sft, powers)
    
    # Check selected power
    beta <- check.sft(sft)
    
    # Construct co-expression similarity
    adj <- WGCNA::adjacency(datExpr=dat, 
                            power=beta, 
                            corFnc=cor.fn, 
                            type="unsigned", 
                            corOptions=cor.options)
    
    # Topological overlap dissimilarity transformation
    dis <- WGCNA::TOMdist(adjMat=adj, TOMType="unsigned")
    
    # Fast hierarchical clustering of dissimilarity
    dendro <- flashClust::flashClust(d=as.dist(dis), method="average")
    
    # Module identification using dynamic tree cut algorithm
    modules <- dynamicTreeCut::cutreeDynamic(dendro=dendro,
                                             method="hybrid",
                                             distM=dis, 
                                             deepSplit=4,
                                             pamRespectsDendro=FALSE,
                                             minClusterSize=min.size)
    
    # Assign module colours
    colors <- WGCNA::labels2colors(labels=modules, zeroIsGrey=TRUE)
    
    # Merging close modules
    if (merging) {
        
        merged <- WGCNA::mergeCloseModules(exprData=dat,
                                           colors=colors,
                                           corFnc=cor.fn,
                                           corOptions=cor.options,
                                           cutHeight=merging.cut)
        
        if (do.plot) WGCNA::plotDendroAndColors(dendro=dendro, 
                                                colors=cbind(colors, merged$colors), 
                                                groupLabels=c("Original Modules", "Merged Modules"), 
                                                dendroLabels=FALSE,
                                                addGuide=TRUE)
        # Merged data
        colors <- merged$colors
        eigengenes <- merged$newMEs
        
    }
    else {

        eigengenes <- WGCNA::moduleEigengenes(expr=dat, colors=colors)$eigengenes
        
        if (do.plot) WGCNA::plotDendroAndColors(dendro=dendro,
                                                colors=colors, 
                                                groupLabels="Modules", 
                                                dendroLabels=FALSE, 
                                                addGuide=TRUE)
    }

    # Formatting
    colnames(eigengenes) <- substr(colnames(eigengenes), 3, 500)

    # Define modules
    genes <- rownames(eset)
    mods <- list()
    for(i in unique(colors)) {
        mods[[i]] <- genes[colors == i]
    }
    
    # Check
    stopifnot(sort(table(colors)) == sort(unlist(lapply(mods, length))))
    
    # For downstream analysis
    return(list(dat=dat,
                beta=beta,
                genes=genes,
                colors=colors,
                mods=mods,
                dendro=dendro,
                eigengenes=eigengenes,
                args=args))
}

# Eigengene for one module from the kth principal component
#' 
#' @export
get.eigengene <- function(dat, mod, k=1) {
    #prcomp(dat[, mod], scale=TRUE, center=TRUE)$x[,k]
    svd(t(dat[, mod]))$v[,k]
}

# Eigengene for all modules from the kth principal component 
#' 
#' @export
get.eigengenes <- function(dat, mods, k=1) {
    lapply(mods, function(mod) get.eigengene(dat, mod, k)) %>%
    do.call(cbind, .) %>%
    as.data.frame() %>%
    magrittr::set_rownames(rownames(dat))
}

# Module membership by correlation with eigengene for one module 
#' 
#' @export
fuzzy.membership <- function(dat, me, fn=bicor) {
    stopifnot(rownames(dat) == rownames(me))
    abs(fn(dat, me))
}

# Extend modules by a quantile of multiple module eigengenes
#' 
#' @export
fuzzy.mods <- function(dat, mods, q=0.5, k=1, fn=bicor) {
    me.x <- lapply(seq(k), function(x) get.eigengenes(dat, mods, x))
    fm.x <- lapply(me.x, function(me) fuzzy.membership(dat, me, fn=fn))
    fm.sum <-  Reduce('+', fm.x)
    
    mapply(function(mod.name, mod.members) {
        fm.mod <- fm.sum[,mod.name]
        fm.members <- fm.mod[mod.members]
        fm.members.fuzzy <- fm.mod[fm.mod >= quantile(fm.members, q)]
        union(mod.members, names(fm.members.fuzzy))
    }, names(mods), mods, SIMPLIFY=FALSE)
}
