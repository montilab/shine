#' A push/pop capable vector
#'
#' @importFrom R6 R6Class
#' 
#' @keywords internal
pvector <- R6Class("pvector", list(
    #' @field values A vector of values
    values = NULL,

    #' @description
    #' Create a pvector
    #' @param values A vector of values
    #' @return A new pvector  
    initialize = function(values=c()) {
        self$values <- values
    },
    #' @description
    #' Print pvector
    #' @return NULL  
    print = function() {
        base::print(self$values)
        invisible(self)
    },
    #' @description
    #' Get length of pvector
    #' @return An integer  
    length = function() {
        base::length(self$values)
    },
    #' @description
    #' Pop vector
    #' @return Popped value 
    pop = function() {
        if (length(self$values) > 0) {
            popped.value <- self$values[1]
            self$values <- self$values[-1]
            return(popped.value)   
        }
    },
    #' @description
    #' Push values
    #' @param pushed.values A vector of values
    #' @return NULL
    push = function(pushed.values) {
        self$values <- c(self$values, pushed.values)
    }
))

#' Filter variables with zero variance in one or more subtypes
#' 
#' @param eset An expression set object
#' @param column The column in pData where subtypes are described
#' @param subtypes One or more unique subtypes
#' @param fn A function required to be non-zero
#' 
#' @importFrom Biobase exprs pData
#'  
#' @export
filter.var <- function(eset, column, subtypes, fn=var) {
    type <- pData(eset)[,column]
    variable.genes <- list()
    for (i in subtypes) {
        eset.sub <- eset[,type == i]
        sub.exprs.mat <- Biobase::exprs(eset.sub)
        sub.genes.keep <- apply(sub.exprs.mat, 1, fn) != 0
        variable.genes[[i]] <- rownames(sub.exprs.mat)[sub.genes.keep]
    }
    Reduce(intersect, variable.genes)
}

#' Select variables by median absolute deviation across one or more subtypes
#' 
#' @param eset An expression set object
#' @param column The column in pData where subtypes are described
#' @param subtypes One or more unique subtypes
#' @param limit Number of genes to select
#' @param genes Allowed genes
#' @param fn A function to rank variables by
#' 
#' @importFrom Biobase exprs pData
#'  
#' @export
select.var <- function(eset, column, subtypes, limit=2500, genes=rownames(eset), fn=mad) {
    type <- pData(eset)[,column]
    ranked.genes <- list()
    for (i in subtypes) {
        eset.sub <- eset[,type == i]
        sub.exprs.mat <- Biobase::exprs(eset.sub[genes,])
        ranked.genes[[i]] <- pvector$new(names(sort(apply(sub.exprs.mat, 1, fn), decreasing=TRUE)))
    }
    
    genes.selected <- c()
    i <- 1
    while (TRUE) {
        popped <- ranked.genes[[i]]$pop()
        if (popped %in% genes.selected) next
        
        genes.selected <- c(genes.selected, popped)
        i <- i+1
        
        if (i == length(ranked.genes)+1) i <- 1
        if (length(genes.selected) >= limit) break
    }
    
    return(genes.selected)
}
