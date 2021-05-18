#' Spring-embedded layout
#' 
#' @param ig An igraph object
#' 
#' @return A graph layout_with_graphopt
#' 
#' @keywords internal
model.layout <- function(ig) {
    igraph::layout_with_graphopt(ig, 
                                 start=NULL, 
                                 niter=500,
                                 charge=0.001,
                                 mass=30, 
                                 spring.length=0,
                                 spring.constant=1, 
                                 max.sa.movement=5)
}

#' Removes loops, multi-edges, and isolated vertices from graph
#' 
#' @param ig An igraph object
#' 
#' @return An igraph object
#' 
#' @importFrom igraph simplify delete_vertices degree
#' 
#' @keywords internal
model.simplify <- function(ig) {
    ig <- igraph::simplify(ig, remove.multiple=TRUE, remove.loops=TRUE)
    ig <- igraph::delete_vertices(ig, which(igraph::degree(ig) == 0))
    return(ig)
}

#' Get a string-representation of igraph edges
#' 
#' @param ig An igraph object
#' @param rev Use true to return reversed edges
#' 
#' @return A character vector of edges
#' 
#' @importFrom igraph as_data_frame
#' @importFrom dplyr mutate %>%
#' 
#' @keywords internal
model.edges <- function(ig, rev=FALSE) {
    df <- igraph::as_data_frame(ig, what="edges") %>%
          dplyr::mutate(fwd=paste(from, to, sep="-")) %>%
          dplyr::mutate(rev=paste(to, from, sep="-"))
    
    if (rev) c(df$fwd, df$rev) else df$fwd
}

#' Rename igraph vertices as characters
#' 
#' @param ig An igraph object
#' @param prefix A prefix for vertex labels
#' 
#' @return An igraph object
#' 
#' @importFrom igraph V set.vertex.attribute
#' 
#' @keywords internal
model.rename <- function(ig, prefix="G") {
    ids <- paste(prefix, igraph::V(ig), sep="")
    ig <- igraph::set.vertex.attribute(ig, "name", value=ids)
    ig <- igraph::set.vertex.attribute(ig, "label", value=ids)
    return(ig)
}

#' Szymkiewicz–Simpson coefficient of edges
#' 
#' @param ig1 An igraph object
#' @param ig2 An igraph object
#' 
#' @return An edge similarity
#' 
#' @importFrom igraph graph_from_adjacency_matrix intersection ecount
#' @importFrom dplyr %>%
#' 
#' @export
model.similarity <- function(ig1, ig2) {
    if (is(ig1, "matrix")) ig1 <- igraph::graph_from_adjacency_matrix(ig1, mode="undirected", diag=FALSE)
    if (is(ig2, "matrix")) ig2 <- igraph::graph_from_adjacency_matrix(ig2, mode="undirected", diag=FALSE)
    
    igraph::intersection(ig1, ig2, byname=FALSE) %>%
    igraph::ecount() %>%
    ( function(x) x / min(igraph::ecount(ig1), igraph::ecount(ig2)) )
}

#' Edge similarity of two or more graphs
#' 
#' @param igs A list of igraph objects
#' 
#' @return An edge similarity
#' 
#' @importFrom dplyr %>%
#' 
#' @export
models.similarity <- function(igs) {
    combn(igs, 2, function(x) {
        model.similarity(x[[1]], x[[2]])
    }) %>%
    mean()
}

#' Simulate diverging models through preferential attachment
#' 
#' @param n The number of divergent graphs to create
#' @param p The number of vertices in the divergent graph 
#' @param seed An igraph object to start algorithm from
#' @param ... Additional arguments passed to /code{igraph::sample_pa}
#' 
#' @return A list of igraph objects
#' 
#' @importFrom igraph sample_pa
#' 
#' @export
model.hpa <- function(n, p, seed, ...) {
    lapply(seq(n), function(x) {
        igraph::sample_pa(n=p, directed=FALSE, start.graph=seed, ...)
    })
}

#' Simulate network via modular preferential attachment model
#' 
#' @param p The number of vertices in the graph
#' @param m The number of modules
#' @param q.hub The quantile for degree cutoff for defining modular hubs  
#' @param m.links The number of links between modules
#' @param power The power of the preferential attachment
#' @param z.appeal The attractiveness of the vertices with no adjacent edges
#' @param ... Additional arguments passed to /code{igraph::sample_pa}
#' 
#' @return An igraph object
#' 
#' @importFrom igraph V E edges degree sample_pa
#' 
#' @export
model.mpa <- function(p=300, m=6, q.hub=0.95, m.links=6, power=1.7, z.appeal=1, ...) {
    # Calculate module sizes
    m.size <- as.integer( (p - m.links) / m)
    
    # Create scale-free modules
    graphs <- lapply(seq_len(m), function(x) {
        g.x <- igraph::sample_pa(m.size, 
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
    return(g)   
}

#' Simulate network via modular erdos-renyi model
#' 
#' @param p The number of vertices in the graph
#' @param m The number of islands in the graph
#' @param m.prob The probability to create each possible edge into each island
#' @param n.inter The number of edges to create between two islands
#' @param ... Additional arguments passed to /code{igraph::sample_islands}
#' 
#' @return An igraph object
#' 
#' @importFrom igraph sample_islands
#' 
#' @export
model.mer <- function(p=300, m=6, m.prob=0.07, n.inter=1, ...) {
    igraph::sample_islands(islands.n=m, 
                           islands.size=p/m,
                           islands.pin=m.prob,
                           n.inter=n.inter,
                           ...)
}

#' Simulate network via preferential attachment model
#' 
#' @param p The number of vertices in the graph
#' @param power The power of the preferential attachment
#' @param z.appeal The attractiveness of the vertices with no adjacent edges
#' @param ... Additional arguments passed to /code{igraph::sample_pa}
#' 
#' @return An igraph object
#' 
#' @importFrom igraph sample_pa
#' 
#' @export
model.pa <- function(p=300, power=1, z.appeal=1, ...) {
    igraph::sample_pa(n=p, 
                      power=power,
                      zero.appeal=z.appeal,
                      directed=FALSE, 
                      ...)
}

#' Simulate network via small world model
#' 
#' @param p The number of vertices in the graph
#' @param dim The dimension of the starting lattice
#' @param nei The neighborhood within which the vertices of the lattice will be connected
#' @param rw.prob The probability for rewiring
#' @param ... Additional arguments passed to /code{igraph::sample_smallworld}
#' 
#' @return An igraph object
#' 
#' @importFrom igraph sample_smallworld
#' 
#' @export
model.sw <- function(p=300, dim=1, nei=3, rw.prob=0.05, ...) {
    igraph::sample_smallworld(size=p, 
                              dim=dim,
                              nei=nei,
                              p=rw.prob,
                              ...)
}

#' Simulate network via forest fire model
#' 
#' @param p The number of vertices in the graph
#' @param fw.prob The probability for foward burning
#' @param bw.factor The backward burning ratio
#' @param ... Additional arguments passed to /code{igraph::sample_forestfire}
#' 
#' @return An igraph object
#' 
#' @importFrom igraph sample_forestfire
#' 
#' @export
model.ff <- function(p=300, fw.prob=0.3, bw.factor=0.5, ...) {
    igraph::sample_forestfire(nodes=p, 
                              fw.prob=fw.prob, 
                              bw.factor=bw.factor,
                              directed=FALSE, 
                              ...)
}

#' Simulate network via erdos-renyi model
#' 
#' @param p The number of vertices in the graph
#' @param e.prob A probability for drawing an edge between two arbitrary vertices
#' @param ... Additional arguments passed to /code{igraph::sample_gnp}
#' 
#' @return An igraph object
#' 
#' @importFrom igraph sample_gnp
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
#' @param vertex.size Vertex size(s)
#' @param vertex.label Vertex labels(s)
#' @param colors Use true to include vertex colors in plot
#' @param ... Additional arguments passed to /code{igraph::plot}
#' 
#' @return A graph visualization
#' 
#' @export
model.plot <- function(ig, seed=1, vertex.size=4, vertex.label=NA, colors=FALSE, ...) {
    set.seed(seed)
    if (colors) {
        plot(ig, 
             vertex.size=vertex.size, 
             vertex.label=vertex.label,
             layout=model.layout(ig),
             ...)
    } 
    else {
        plot(ig, 
             vertex.size=vertex.size,
             vertex.label=vertex.label,
             vertex.color="grey",
             edge.color="black", 
             layout=model.layout(ig),
             ...)
    }
}

#' Plot one or more models highlighting differences
#' 
#' @param igs A list of igraph objects
#' 
#' @return A graph visualization
#' 
#' @export
models.plot <- function(igs) {
    edges <- lapply(igs, function(x) model.edges(x, rev=TRUE))
    shared <- Reduce(intersect, edges)
    par(mfrow=c(1,length(igs)))
    for (ig in igs) {
        igraph::E(ig)$width <- 2
        igraph::E(ig)$color <- "#FC4840" # blue
        igraph::E(ig)$color[model.edges(ig) %in% shared] <- "#113D5C" # red
        plot(ig, vertex.size=4, vertex.color="grey", vertex.label=NA, layout=model.layout(ig))
    }
}

#' Plot model degree distribution
#' 
#' @param ig An igraph object
#' 
#' @return A plot
#' 
#' @importFrom igraph degree
#' 
#' @export
model.kd <- function(ig) {
    k <- igraph::degree(ig)
    pk <- prop.table(table(k))
    df <- data.frame(k=as.numeric(names(pk)), pk=as.numeric(pk))
    df <- df[df$k > 0,]
    plot(x=df$k, y=df$pk, log="xy",
         main="Degree Distribution",
         ylab="P(k)",
         xlab="k")
}

#' Simulate multivariate Gaussian data for a model
#' 
#' @param n The number of samples to generate data for
#' @param ig An igraph object for the graph structure
#' @param seed A number to seed the random data
#' 
#' @return The true graph structure and multivariate gaussian data
#' 
#' @importFrom BDgraph bdgraph.sim
#' @importFrom Biobase ExpressionSet
#' @importFrom igraph get.adjacency
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

    return(list(bdg=bdg, eset=eset))
}
