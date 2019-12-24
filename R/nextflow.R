#' Plot a hierarchy description
#' 
#' @param hierarchy A string describing the hierarchy in dot language notation
#' 
#' @importFrom magrittr %>%
#' @importFrom DiagrammeR grViz
#' @importFrom stringr str_interp
#' 
#' @export
plt.hierarchy <- function(hierarchy) {
    "digraph dot {
     graph [layout = dot]
     node [shape = circle,
           style = filled,
           color = black,
           style = filled,
           fillcolor = AliceBlue,
           fontname = Futura]
     edge [color = black]
     ${hierarchy}
    }" %>%
    str_interp(list(hierarchy=hierarchy)) %>%
    grViz()
}

#' Generate a nextflow workflow
#' 
#' @param hierarchy A string describing the hierarchy in dot language notation
#' @param condition The column in pData where subtypes are described
#' @param eset An expression set object
#' @param blanket A gene x gene constraint matrix
#' @param outdir A directory where the workflow will be generated
#' @param iter The number of iterations for the sampling algorithm 
#' @param cores The number of cores for each network to use for parallel execution
#' 
#' @importFrom Biobase pData
#' @importFrom reticulate import_from_path
#' 
#' @export
build.workflow <- function(hierarchy,
                           condition,
                           eset,
                           blanket=NULL,
                           outdir=".",
                           iter=5000,
                           cores=1) {
    
    # Check eset
    stopifnot(is(eset, "ExpressionSet"))
    stopifnot(nrow(eset) >= 3)
    stopifnot(condition %in% colnames(Biobase::pData(eset)))
    
    # Check blanket
    if (!is.null(blanket)) {
        stopifnot(is(blanket, "matrix"))
        stopifnot(colnames(blanket) == rownames(blanket))
        stopifnot(colnames(blanket) == rownames(eset))
    }

    # Check hierarchy
    steps <- unlist(strsplit(gsub(" ", "", hierarchy), "\n"))
    steps <- steps[sapply(steps, function(x) grepl(pattern="->", x))]
    nets <- unname(unlist(sapply(steps, function(x) strsplit(x, "->"))))
    subtypes <- unique(unlist(sapply(nets, function(x) strsplit(x, "_"))))
    subtypes.available <- unique(Biobase::pData(eset)[,condition])
    subtypes.missing <- !(subtypes %in% subtypes.available)
    if (any(subtypes.missing)) stop(paste("The following are missing from eset$", condition, " : ",
                                    paste(subtypes[subtypes.missing], collapse=" "), sep=""))
    
    # Package data
    data <- list(eset=eset, 
                 blanket=blanket,
                 condition=condition,
                 iter=iter,
                 cores=cores)
    
    # Create workflow directory
    dir.create(outdir, showWarnings=FALSE)
    file.copy(system.file("nextflow", package="shine"), outdir, recursive=TRUE)
    
    # Save data package
    dir.create(file.path(outdir, "nextflow/data"), showWarnings=FALSE)
    saveRDS(data, file.path(outdir, "nextflow/data/data.rds"))

    # Build pipeline
    nf <- reticulate::import_from_path("nextflow", system.file("python", package="shine"))
    nf$nf_build(hierarchy,
                data=file.path(outdir, "nextflow/data/data.rds"),
                outdir=file.path(outdir, "nextflow"))
}
