library(BDgraph)
library(Biobase)
library(usethis)

# Simulate Scale-Free Graph
set.seed(1)
graph <- sim.graph(100)

# Simulate Expression Data
set.seed(1)
bdg <- bdgraph.sim(n=30, graph=as.matrix(igraph::get.adjacency(graph)), type="Gaussian", vis=F)

# Expression set
edat <- t(bdg$data)
rownames(edat) <- paste("G", seq(nrow(edat)), sep="")
colnames(edat) <- paste("S", seq(ncol(edat)), sep="")
toy <- ExpressionSet(edat)
toy$subtype <- c(rep("A", 10), rep("B", 10), rep("C", 10))

usethis::use_data(toy)
