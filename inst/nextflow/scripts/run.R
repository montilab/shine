library(shine)

# Parse args
args = commandArgs(trailingOnly=TRUE)

# Argument #1
path.to.data <- args[1]

stopifnot(file.exists(path.to.data))
data <- readRDS(path.to.data)

# Static variables across networks
eset <- data$eset
blanket <- data$blanket
condition <- data$condition
iter <- data$iter
cores <- data$cores

# Argument #2
path.to.prior <- args[2]
if (path.to.prior != "/") {
    bdg <- readRDS(file=path.to.prior)
    prior <- BDgraph::plinks(bdg, round=10)
    rownames(prior) <- colnames(prior)
} else {
    prior <- NULL
}

# Arguments #3/4
label <- args[3]
include <- args[4:length(args)]

# Log data
log <- paste(label, ".log", sep="")
sink(log, append=FALSE, split=TRUE)

cat(paste("Parameters"))
cat(paste("\nLabel:", label))
cat(paste("\nData:", path.to.data))
cat(paste("\nPrior:", path.to.prior))
cat(paste("\nCondition:", condition))
cat(paste("\nInclude:", paste(include, collapse=", ")))
cat(paste("\nIter:", iter))
cat(paste("\nCores:", cores))

# Setup constraints and prior
if (!is.null(blanket) & !is.null(prior)) {
    g.prior <- blanket.inform(blanket, prior) 
} else if (!is.null(blanket) & is.null(prior)) {
    g.prior <- blanket
} else if (is.null(blanket) & !is.null(prior)) {
    g.prior <- prior
} else {
    g.prior <- 0.5
}

# Subset on samples
eset.sub <- eset[, Biobase::pData(eset)[,condition] %in% include]

cat(paste("\n\nData:\n"))
print(dim(eset.sub))

cat(paste("\nSamples:"))
print(table(Biobase::pData(eset.sub)[,condition]))

cat(paste("\nEstimate\n"))
# Estimate
bdg <- BDgraph::bdgraph.mpl(data=t(Biobase::exprs(eset.sub)), 
                            g.prior=g.prior, 
                            method="ggm", 
                            iter=iter, 
                            print=round(iter/25), 
                            cores=cores-1, 
                            save=FALSE)
sink()

# Output graph to current directory
saveRDS(bdg, paste(label, ".rds", sep=""))
