library(shine)

args = commandArgs(trailingOnly=TRUE)
print(args)

# File paths
path.to.eset    <- args[1]
path.to.genes   <- args[2]
path.to.prior   <- args[3]
path.to.blanket <- args[4]

# Sample selection
iters <- as.numeric(args[5])
cores <- as.numeric(args[6])
condition <- args[7]
label <- args[8]
include <- args[9:length(args)]

# Log data
log <- paste(label, ".log", sep="")
write(paste("Parameters...", file=log))
write(paste("Eset:", path.to.eset), file=log, append=TRUE)
write(paste("Genes:", path.to.genes), file=log, append=TRUE)
write(paste("Prior:", path.to.prior), file=log, append=TRUE)
write(paste("Blanket:", path.to.blanket), file=log, append=TRUE)
write(paste("Iters:", iters), file=log, append=TRUE)
write(paste("Cores:", cores), file=log, append=TRUE)
write(paste("Condition:", condition), file=log, append=TRUE)
write(paste("Label:", label), file=log, append=TRUE)
write(paste("Include:", paste(include, collapse=" ")), file=log, append=TRUE)

# Print data
cat("Parameters...", "\n")
cat("Eset:", path.to.eset, "\n")
cat("Genes:", path.to.genes, "\n")
cat("Prior:", path.to.prior, "\n")
cat("Blanket:", path.to.blanket, "\n")
cat("Iters:", iters, "\n")
cat("Cores:", cores, "\n")
cat("Condition:", condition, "\n")
cat("Label:", label, "\n")
cat("Include:", paste(include, collapse=" "), "\n")

# Start inference
eset.raw <- readRDS(file=path.to.eset)
genes <- readRDS(file=path.to.genes)

# Subset the data
which.genes <- rownames(eset.raw) %in% genes
which.include <- Biobase::pData(eset.raw)[,condition] %in% include 
eset.sub <- eset.raw[which.genes, which.include]
print(dim(eset.sub))
print(table(Biobase::pData(eset.sub)[,condition]))
edat <- t(Biobase::exprs(eset.sub))

# Check prior and constraints
if (path.to.blanket != "/") {
    blanket <- as.matrix(readRDS(file=path.to.blanket))
    rownames(blanket) <- colnames(blanket)
} else {
    blanket <- NULL
}
if (path.to.prior != "/") {
    bdg <- readRDS(file=path.to.prior)
    prior <- BDgraph::plinks(bdg, round=10)
    rownames(prior) <- colnames(prior)
} else {
    prior <- NULL
}

if (!is.null(blanket) & !is.null(prior)) {
    g.prior <- blanket.inform(blanket, prior) 
} else if (!is.null(blanket) & is.null(prior)) {
    g.prior <- blanket
} else if (is.null(blanket) & !is.null(prior)) {
    g.prior <- prior
} else {
    g.prior <- 0.5
}

cat("Estimating network...\n")

bdg <- BDgraph::bdgraph.mpl(data=t(Biobase::exprs(eset.sub)), 
                            g.prior=g.prior, 
                            method="ggm", 
                            iter=iters, 
                            print=round(iters/25), 
                            cores=cores-1, 
                            save=F)

# Output graph to current directory
saveRDS(bdg, paste(label, ".rds", sep=""))
