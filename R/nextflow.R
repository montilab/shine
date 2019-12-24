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
#' @param path.eset A full path to an expression set object saved as an rds file
#' @param path.genes A full path to a character vector of genes saved as an rds file
#' @param path.blanket A full path to a gene x gene constraint matrix saved as an rds file
#' @param outdir A directory where the workflow will be generated
#' @param iters The number of iterations for the sampling algorithm 
#' @param cores The number of cores for each network to use for parallel execution
#' 
#' @importFrom Biobase pData
#' @importFrom reticulate import_from_path
#' 
#' @export
build.workflow <- function(hierarchy,
                           condition,
                           path.eset,
                           path.genes,
                           path.blanket=NULL,
                           outdir=".",
                           iters=5000,
                           cores=1) {
    
    # Checks
    stopifnot(file.exists(path.eset))
    stopifnot(file.exists(path.genes))
    
    # Required paths
    eset <- readRDS(path.eset)
    genes <- readRDS(path.genes)
    
    # Check fData and pData
    if (!(is(eset, "ExpressionSet"))) stop("Not a valid expression set object")
    if (any(!(genes %in% rownames(eset)))) stop("All genes must be in expression set rownames")
    if (!(condition %in% colnames(pData(eset)))) stop("Cannot find eset$", condition, sep="")
    
    # Hiearchy parsing
    steps <- unlist(strsplit(gsub(" ", "", hierarchy), "\n"))
    steps <- steps[sapply(steps, function(x) grepl(pattern="->", x))]
        
    nets <- unname(unlist(sapply(steps, function(x) strsplit(x, "->"))))
    subtypes <- unique(unlist(sapply(nets, function(x) strsplit(x, "_"))))
    
    subtypes.available <- unique(pData(eset)[,condition])
    subtypes.missing <- !(subtypes %in% subtypes.available)
    if (any(subtypes.missing)) stop(paste("The following are missing from eset$", condition, " : ",
                                    paste(subtypes[subtypes.missing], collapse=" "), sep=""))
    
    # Blanket
    if (is.null(path.blanket)) {
        path.blanket <- "/"
    } else {
        stopifnot(file.exists(path.blanket))
        blanket <- readRDS(path.blanket)
        stopifnot(is(blanket, "matrix"))
        if (any(rownames(blanket) != genes) | any(colnames(blanket) != genes)) stop("Blanket must be a gene x gene matrix")
    }
    
    # Output
    dir.create(outdir, showWarnings=FALSE)
    file.copy(system.file("nextflow", package="shine"), outdir, recursive=TRUE)
    
    # Build pipeline
    nf <- reticulate::import_from_path("nextflow", system.file("python", package="shine"))
    nf$nf_build(outdir=file.path(outdir, "nextflow"),
                hierarchy=hierarchy,
                condition=condition,
                path_eset=path.eset,
                path_genes=path.genes,
                path_blanket=path.blanket,
                iters=iters,
                cores=cores)
}
