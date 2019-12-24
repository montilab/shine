#' Plot an igraph object
#' 
#' @param ig An igraph object
#' @param seed A number to seed random layouts
#' @return A graph visualization
#' 
#' @export
plt.graph <- function(ig, seed=1) {
    set.seed(seed)
    plot(ig, vertex.color="grey", vertex.size=4, vertex.label=NA)
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
