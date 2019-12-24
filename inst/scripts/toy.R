library(BDgraph)
library(Biobase)

# Simulate Scale-Free Graph
set.seed(1)
graph <- sim.graph(15)
plot(graph, vertex.color="white", vertex.size=15)

# Simulate Expression Data
set.seed(1)
bdg <- bdgraph.sim(n=30, graph=as.matrix(igraph::get.adjacency(graph)), type="Gaussian", vis=F)

# Expression set
edat <- t(bdg$data)
rownames(edat) <- paste("G", seq(nrow(edat)), sep="")
colnames(edat) <- paste("S", seq(ncol(edat)), sep="")
eset <- ExpressionSet(edat)
eset$subtype <- c(rep("A", 10), rep("B", 10), rep("C", 10))

# Check
dim(eset)
rownames(eset)
colnames(eset)
rownames(fData(eset))
pData(eset)

# Genes
genes <- rownames(eset)

# Constraints
blanket <- matrix(0.5, nrow=15, ncol=15, dimnames=list(genes, genes))
print(blanket)
M1 <- c(7,8,9,15)
M2 <- c(5,6,10,11)
blanket[M1, M2] <- 0
blanket[M2, M1] <- 0
print(blanket)

# Save
saveRDS(eset, file.path(system.file("extdata", package="shine"), "eset.rds"))
saveRDS(genes, file.path(system.file("extdata", package="shine"), "genes.rds"))
saveRDS(blanket, file.path(system.file("extdata", package="shine"), "blanket.rds"))
