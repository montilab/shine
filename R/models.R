#' Repeatable vector of distinct colors
#' 
#' @keywords internal
rcolors <- function(reps=1) {
    rep(c("#BDBE4D","#E8EE77","#CA874D","#945DCF","#76DDEC",
          "#5BDE78","#D6A0E8","#CCADB5","#719DEB","#819F9E",
          "#8671B3","#E6ECDC","#68C495","#EBEAB7","#DDEC3C",
          "#6AA03A","#D782E9","#E29B91","#ECBD42","#E294B6",
          "#E176C6","#EFD7E8","#6D79E5","#AF9D7A","#97E088",
          "#6C9A61","#56EB4D","#61F0BA","#B750E5","#9DB6E3",
          "#7037E4","#6FEDD7","#785867","#D0B6E2","#E438E9",
          "#DCC27C","#E05E38","#A5E753","#56B7E4","#D7EA97",
          "#C0D0E0","#B7E7DE","#6283A0","#EAC5A8","#D73D7B",
          "#61BDBB","#CE636D","#E34CC6","#B4EDBF","#ABC298"), reps)
}

#' Rename igraph vertices as characters
#' 
#' @keywords internal
model.rename <- function(ig, prefix="G") {
    ids <- paste(prefix, igraph::V(ig), sep="")
    igraph::V(ig)$name <- ids
    igraph::V(ig)$label <- ids
    return(ig)
}

#' Simulate network via modular preferential attachment model
#' 
#' @export
model.mpa <- function(p=300, m=6, q.hub=0.95, m.links=6, power=1.7, z.appeal=1, ...) {
    # Calculate module sizes
    m.size <- as.integer( (p - m.links) / m)
    
    # Create scale-free modules
    graphs <- lapply(seq_len(m), function(x) {
        g.x <- sample_pa(m.size, 
                         power=power, 
                         zero.appeal=z.appeal, 
                         directed=FALSE, 
                         ...)
        
        igraph::V(g.x)$module <- x
        igraph::V(g.x)$hub <- igraph::degree(g.x) >= quantile(igraph::degree(g.x), q.hub)
        igraph::V(g.x)$target <- !V(g.x)$hub
        igraph::V(g.x)$color <- rcolors(10)[x]
        
        return(g.x)
    })
    
    # Merge modules into a single graph
    g <- Reduce("+", graphs)
    
    pairs <- sapply(seq_len(m.links), function(x) {
        
        # Randomly choose two different modules to connect
        random.modules <- sample(unique(igraph::V(g)$module), size=2, replace=F)
        
        m1 <- as.numeric(random.modules[1])
        m2 <- as.numeric(random.modules[2])
        
        # Identify hub node of module 1
        m1.hubs <- as.character(igraph::V(g)[igraph::V(g)$hub & igraph::V(g)$module == m1])
        m2.targets <- as.character(igraph::V(g)[igraph::V(g)$target & igraph::V(g)$module == m2])

        c(sample(m1.hubs, size=1), sample(m2.targets, size=1))
    })
    
    # Add in new links colored red
    igraph::E(g)$color <- "black"
    e <- igraph::edges(pairs)
    e$color <- "red"
    
    g <- g + e
    g <- model.rename(g)
    return(g)   
}

#' Simulate network via modular erdos-renyi model
#' 
#' @export
model.mer <- function(p=300, m=6, m.prob=0.07, n.inter=1, ...) {
    igraph::sample_islands(islands.n=m, 
                           islands.size=p/m,
                           islands.pin=m.prob,
                           n.inter=n.inter,
                           ...) %>%
    model.rename()
}

#' Simulate network via preferential attachment model
#' 
#' @export
model.pa <- function(p=300, power=1, z.appeal=1, ...) {
    igraph::sample_pa(n=p, 
                      power=power,
                      zero.appeal=z.appeal,
                      directed=FALSE, 
                      ...) %>%
    model.rename()
}

#' Simulate network via small world model
#' 
#' @export
model.sw <- function(p=300, dim=1, nei=3, rw.prob=0.05, ...) {
    igraph::sample_smallworld(size=p, 
                              dim=dim,
                              nei=nei,
                              p=rw.prob,
                              ...) %>%
    model.rename()
}

#' Simulate network via forest fire model
#' 
#' @export
model.ff <- function(p=300, fw.prob=0.3, bw.factor=0.5, ...) {
    igraph::sample_forestfire(nodes=p, 
                              fw.prob=fw.prob, 
                              bw.factor=bw.factor,
                              directed=FALSE, 
                              ...) %>%
    model.rename()
}

#' Simulate network via erdos-renyi model
#' 
#' @export
model.er <- function(p=300, e.prob=0.01, ...) {
    igraph::sample_gnp(n=p, 
                       p=e.prob,
                       directed=FALSE,
                        ...)
}

#' Plot model
#' 
#' @param ig An igraph object
#' @param seed A number to seed random layouts
#' @return A graph visualization
#' 
#' @export
model.plot <- function(ig, seed=1, colors=FALSE) {
    set.seed(seed)
    if (colors) {
        plot(ig, 
            vertex.size=4, 
            vertex.label=NA)
    } 
    else {
        plot(ig, 
            vertex.size=4, 
            vertex.color="orange",
            edge.color="black", 
            vertex.label=NA)
    }
}

#' Plot model degree distribution
#' 
#' @param ig An igraph object
#' @return A ggplot object
#' 
#' @import ggplot2
#' 
#' @export
model.kd <- function(ig) {
    k <- igraph::degree(ig)
    k.nz <- k[k > 0]
    df <- data.frame(k=seq(0, max(k.nz)), pk=igraph::degree_distribution(ig))
    df <- df[df$pk > 0, ]
    ggplot(df) +
    geom_point(aes(x=k, y=pk)) + 
    scale_y_continuous(trans='log10') +
    scale_x_continuous(trans='log10') +
    ylab("P(k)") +
    xlab("k") +
    ggtitle("Degree Distribution") +
    theme_bw()
}

#' Simulate multivariate gaussian data for a model
#' 
#' @export
model.sim <- function(n, ig, seed=1) {
    # Extract adjacency matrix
    adj <- as.matrix(igraph::get.adjacency(ig))

    # Simulate multivariate gaussian
    set.seed(seed)
    bdg <- BDgraph::bdgraph.sim(n=n, graph=adj, type="Gaussian", vis=F)

    # Wrap data in expression set object
    eset <- Biobase::ExpressionSet(t(bdg$data))
    colnames(eset) <- paste("S", colnames(eset), sep="")
    rownames(eset) <- colnames(adj)
    return(list(bdg=bdg, eset=eset))
}

#' Get a string-representation of igraph edges
#' 
#' @param ig An igraph object
#' @param rev Use true to return reversed edges
#' @return A character vector of edges
#' 
#' @importFrom igraph as_data_frame
#' @importFrom dplyr mutate %>%
#' 
#' @keywords internal
string.edges <- function(ig, rev=FALSE) {
    df <- igraph::as_data_frame(ig, what="edges") %>%
          dplyr::mutate(fwd=paste(from, to, sep="-")) %>%
          dplyr::mutate(rev=paste(to, from, sep="-"))
    
    if (rev) c(df$fwd, df$rev) else df$fwd
}

#' Plot multiple igraph objects
#' 
#' @param ... One or more igraph objects
#' @return A graph visualization
#' 
#' @importFrom igraph E
#' 
#' @export
plt.graphs <- function(...) {
    graphs <- list(...)
    edges <- lapply(graphs, function(x) string.edges(x, rev=TRUE))
    shared <- Reduce(intersect, edges)
    par(mfrow=c(1,length(graphs)))
    for (ig in graphs) {
        igraph::E(ig)$width <- 2
        igraph::E(ig)$color <- "#FC4840"
        igraph::E(ig)$color[string.edges(ig) %in% shared] <- "#113D5C"      
        plt.graph(ig)
    }
}

#' Szymkiewicz–Simpson coefficient of edges
#' 
#' @param ig1 An igraph object
#' @param ig2 An igraph object
#' @return An edge similarity
#' 
#' @importFrom igraph graph_from_adjacency_matrix intersection ecount
#' @importFrom dplyr %>%
#' 
#' @export
graph.similarity <- function(ig1, ig2) {
    if (is(ig1, "matrix")) ig1 <- igraph::graph_from_adjacency_matrix(ig1, mode="undirected", diag=FALSE)
    if (is(ig2, "matrix")) ig2 <- igraph::graph_from_adjacency_matrix(ig2, mode="undirected", diag=FALSE)
    
    igraph::intersection(ig1, ig2, byname=FALSE) %>%
    igraph::ecount() %>%
    ( function(x) x / min(ecount(ig1), ecount(ig2)) )
}

#' Barabasi-Albert algorithm for scale-free graphs
#' 
#' @param p The number of vertices
#' @param seed The starting igraph for the preferential attachment algorithm
#' @param ... One or morge arguments passed to the algorithm
#' @return An igraph object
#' 
#' @importFrom igraph sample_pa
#' 
#' @export
sim.graph <- function(p, seed=NULL, ...) {
    igraph::sample_pa(n=p, directed=FALSE, algorithm="psumtree", start.graph=seed, ...)
}

#' Simulate a range of scale-free graphs of different similarity
#' 
#' @param n The number of graphs
#' @param p.graphs The number of desired vertices
#' @param p.start The number of vertices in the seed graph for the least similar graphs
#' @param p.step The step size between the seed graphs for the least and most similar grpahs
#' @return A list of graphs with decreasing similarity
#' 
#' @importFrom dplyr %>%
#' 
#' @export
sim.graphs <- function(n, p.graphs, p.start=p.graphs/2, p.step=(p.start-p.graphs)/10, ...) {
    
    # Decreasing length of diverging preferential attachment
    iter <- seq(p.graphs, p.start, as.integer(p.step))
    names(iter) <- paste("P", iter, sep="")
    
    simulation <- lapply(iter, function(p) {
        
        # Generate seed
        seed <- sim.graph(p)
        
        # Extend graphs
        graphs <- lapply(seq(n), function(x) {
            sim.graph(p=p.graphs, seed=seed, ...)
        })
        
        # Calculate pairwise similarity of graphs
        similarity <- combn(graphs, 2, function(x) {
            graph.similarity(x[[1]], x[[2]])
        }) %>%
        mean()
        
        list(seed=seed, graphs=graphs, similarity=similarity)
    })
    return(simulation)
}
