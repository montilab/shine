#' Plot curve for soft thresholding
#' 
#' @param sft Table of values for scale-free fit
#' @param powers The power vector tested
#' 
#' @return A plot
#'  
#' @keywords internal
sft.plot <- function(sft, powers) {
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
#' 
#' @return A theshold
#' 
#' @keywords internal
sft.check <- function(sft) {
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
#' 
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
    
    if (do.plot) sft.plot(sft, powers)
    
    # Check selected power
    beta <- sft.check(sft)
    
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

#' Module eigengene from the kth principal component for a module
#'
#' @param dat Expression data
#' @param mod A module of genes
#' @param k The principal component to use
#' 
#' @return An eigengene for the kth principle component
#' 
#' @keywords internal
eig.get <- function(dat, mod, k=1) {
    svd(t(dat[, mod]))$v[,k]
}

#' Variance explained by eigengen from the kth principal component for a module
#'
#' @param dat Expression data
#' @param mod A module of genes
#' @param k The principal component to use
#' 
#' @return Variance explained by the kth principle component
#' 
#' @keywords internal
eig.exp <- function(dat, mod, k=1) {
    d <- svd(t(dat[, mod]))$d
    prop <- d^2/sum(d^2)
    prop[k]
}

#' Module eigengenes from the kth principal component for all modules
#'
#' @param dat Expression data
#' @param mods A list of modules
#' @param k The principal component to use
#' 
#' @return A matrix of module eigengenes
#' 
#' @importFrom magrittr set_rownames
#' 
#' @export
me.get <- function(dat, mods, k=1) {
    lapply(mods, function(mod) eig.get(dat, mod, k)) %>%
    do.call(cbind, .) %>%
    as.data.frame() %>%
    magrittr::set_rownames(rownames(dat))
}

#' Module memberships by correlation with eigengene for a module 
#'
#' @param dat Expression data
#' @param mods A list of modules
#' @param k The principal component to use
#' @param fn A correlation function
#' 
#' @return A matrix of module memberships
#' 
#' @importFrom WGCNA bicor
#' 
#' @export
mm.get <- function(dat, mods, k=1, fn=bicor) {
    me <- me.get(dat, mods, k)
    abs(fn(dat, me))
}

#' Wrapper to merge module memberships for the first two eigengenes
#'
#' @param dat Expression data
#' @param mods A list of modules
#' @param fn A correlation function
#' 
#' @return A dataframe of module memberships
#' 
#' @importFrom WGCNA bicor
#' @importFrom reshape2 melt
#' @importFrom magrittr set_colnames
#' 
#' @export
mm.merged <- function(dat, mods, fn=bicor) {
    mm.1 <- mm.get(dat, mods, k=1, fn=fn)
    mm.2 <- mm.get(dat, mods, k=2, fn=fn)
    
    mm.1.melt <- reshape2::melt(as.matrix(mm.1))
    mm.2.melt <- reshape2::melt(as.matrix(mm.2))
    
    merge(mm.1.melt, mm.2.melt, by.x=c("Var1", "Var2"), by.y=c("Var1", "Var2")) %>%
    magrittr::set_colnames(c("gene", "module", "MM1", "MM2"))
}

#' Plot membership of all genes for one module
#'
#' @param dat Expression data
#' @param mods A list of modules
#' @param mod A module name
#' @param colors Gene colors
#' @param size Size of data points
#' @param fn A correlation function
#' 
#' @return A ggplot object
#' 
#' @import ggplot2
#' @importFrom dplyr %>% mutate filter
#' @importFrom scales percent
#' 
#' @export
mod.plot <- function(dat, mods, mod, colors, size=1, fn=bicor) {
    
    # Variance explained
    ex.1 <- eig.exp(dat, mods[[mod]], k=1)
    ex.2 <- eig.exp(dat, mods[[mod]], k=2)

    mm.x.melt <- mm.merged(dat, mods, fn) %>% 
                 dplyr::mutate(primary = colors[gene]) %>%
                 dplyr::mutate(membership = as.integer(primary == module)) %>%
                 dplyr::filter(module == mod)
    
    color.values <- setNames(unique(colors), unique(colors))
    ggplot(mm.x.melt, aes(x=MM1, y=MM2, fill=primary)) +
    geom_point(size=size, shape=21, color="black") +
    scale_fill_manual(values=color.values) +
    ggtitle(mod) +
    xlab(paste("MM1 ", "(", percent(ex.1), ")", sep="")) + 
    ylab(paste("MM2 ", "(", percent(ex.2), ")", sep="")) +
    theme_bw()
}

#' Plot membership of all genes for all modules
#' 
#' @param dat Expression data
#' @param mods A list of modules
#' @param colors Gene colors
#' @param size Size of data points
#' @param ncol Number of facet columns
#' @param legend Use true to include legend
#' @param fn A correlation function
#' 
#' @return A ggplot object
#' 
#' @import ggplot2
#' @importFrom dplyr %>% mutate filter
#' 
#' @return A ggplot object
#' 
#' @export
mods.plot <- function(dat, mods, colors, size=1, ncol=round(sqrt(length(mods))), legend=FALSE, fn=bicor) {
    
    mm.x.melt <- mm.merged(dat, mods, fn) %>% 
                 magrittr::set_colnames(c("gene", "module", "MM1", "MM2")) %>%
                 dplyr::mutate(primary = colors[gene]) %>%
                 dplyr::mutate(membership = as.integer(primary == module))
    
    color.values <- setNames(unique(colors), unique(colors))
    p <- ggplot(mm.x.melt, aes(x=MM1, y=MM2, fill=primary)) +
         geom_point(size=size, shape=21, color="black") +
         scale_fill_manual(values=color.values) +
         facet_wrap(~module, ncol=ncol, scales="free") +
         theme_bw()

    if (!legend) p <- p + theme(legend.position="none")
    return(p)
}

#' Classify fuzzy membership for a module with quadratic discriminant analysis
#' 
#' @param dat Expression data
#' @param mods A list of modules
#' @param mod A module name
#' @param fn A correlation function
#' 
#' @return The model and a dataframe of classifications
#' 
#' @import ggplot2
#' @importFrom dplyr %>% mutate
#' @importFrom reshape2 melt
#' @importFrom MASS qda
#' 
#' @export
fuzzy.predict <- function(dat, mods, mod, fn=bicor) {

    mm.1 <- mm.get(dat, mods, k=1, fn=fn)
    mm.2 <- mm.get(dat, mods, k=2, fn=fn)    

    data <- data.frame(MM1=mm.1[,mod], MM2=mm.2[,mod], stringsAsFactors=F) %>%
            dplyr::mutate(gene = rownames(mm.1)) %>%
            dplyr::mutate(membership = gene %in% mods[[mod]]) %>%
            dplyr::mutate(membership = as.factor(as.integer(membership)))
    
    # Quadratic Discriminant Analysis
    dat <- data[,c("MM1", "MM2", "membership")]
    model <- MASS::qda(membership ~ ., data=dat, prior=c(0.5, 0.5))
    predictions <- predict(model, dat, type="class")
    df <- cbind(data, predictions$posterior[,"1"], predictions$class)
    colnames(df) <- c("MM1", "MM2", "gene", "membership", "posterior", "prediction")
    
    return(list(df=df, model=model))
}

#' Visualize fuzzy module classification
#' 
#' @param dat Expression data
#' @param mods A list of modules
#' @param mod A module name
#' @param size Size of data points
#' @param resolution Resolution of heatmap
#' @param fn A correlation function
#' 
#' @return A ggplot object
#' 
#' @import ggplot2
#' @importFrom dplyr %>% mutate case_when
#' 
#' @export
fuzzy.plot <- function(dat, mods, mod, size=2.5, resolution=100, fn=bicor) {
    
    fp <- fuzzy.predict(dat, mods, mod, fn=fn)
    
    # Extract model and data
    df <- fp$df
    model <- fp$model

    # Generate values for the entire map
    r <- sapply(df[,1:2], range, na.rm=TRUE)
    x <- seq(r[1,1], r[2,1], length.out=resolution)
    y <- seq(r[1,2], r[2,2], length.out=resolution)
    g <- cbind(rep(x, each=resolution), rep(y, time=resolution))
    colnames(g) <- colnames(r)
    g <- as.data.frame(g)
    
    # Predict probabilities of those values
    p <- predict(model, g)
    z <- matrix(p$posterior[,"1"], nrow=resolution, byrow=TRUE)
    mat <- data.frame(x=rep(x, length(y)), y=rep(y, each=length(x)), posterior=as.vector(z))
    
    memberships <- c("Primary Member" = "#2F9599", # Teal
                     "Fuzzy Member"   = "#F7DB4F", # Yellow
                     "Lost Member"    = "#EC2049", # Red
                     "Outsider"       = "#474747") # Grey

    df %<>% dplyr::mutate(status = dplyr::case_when(membership == 1 & prediction == 1 ~ "Primary Member",
                                                    membership == 0 & prediction == 1 ~ "Fuzzy Member",  
                                                    membership == 1 & prediction == 0 ~ "Lost Member",  
                                                    TRUE ~ "Outsider")) %>%
            dplyr::mutate(status = factor(status, levels=names(memberships)))
            
    ggplot(df, aes(x=MM1, y=MM2)) +
    geom_tile(data=mat, aes(x=x, y=y, fill=posterior), alpha=0.90) +
    scale_fill_gradient(low="#FFF9F1", high="#BACBD2") +
    geom_point(size=size, aes(color=status)) +
    scale_color_manual(values=memberships) + 
    theme_classic() +
    ggtitle(mod)
}

#' Classify fuzzy modules with quadratic disciminant analysis
#' 
#' @param dat Expression data
#' @param mods A list of modules
#' @param cores The number of cores for parallel execution
#' @param fn A correlation function
#' 
#' @return A list of fuzzy modules
#' 
#' @importFrom dplyr %>% filter pull
#' @importFrom parallel mcmapply
#' 
#' @export
fuzzy.mods <- function(dat, mods, cores=1, fn=bicor) {
    
    mcmapply(function(mod.name, mod.members) {
        
        fp <- fuzzy.predict(dat, mods, mod=mod.name, fn=fn)
        fuzzy.members <- fp$df %>%
                         dplyr::filter(prediction == 1) %>%
                         dplyr::pull(gene)
        
        union(mod.members, fuzzy.members)
        
    }, names(mods), mods, SIMPLIFY=FALSE, mc.cores=cores)
}
