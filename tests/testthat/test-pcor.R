test_that("blanket is working", {
    # Initalizing a blanket
    genes <- c("G1", "G2", "G3", "G4", "G5", "G6")
    
    # Empty matrices
    blanket <- blanket.new(genes, val=0)
    expect_equal(colnames(blanket), genes)
    expect_equal(rownames(blanket), genes)
    expect_equal(sum(blanket), 0)
    
    # Matrices with values
    blanket <- blanket.new(genes, val=0.5)
    expect_equal(sum(blanket == 0.5), length(genes)^2)
    
    # Filling in modules
    mods <- list("M1"=c("G1", "G2"), "M2"=c("G3", "G4", "G5"), "M3"=c("G6"))
    
    blanket <- blanket.new(genes, val=0)
    blanket <- blanket.lift(blanket, mods, val=0.5)
    expect_equal(sum(blanket == 0.5), 14)
    
    # Randome module pairs
    pairs <- data.frame("M1", "M2", stringsAsFactors=FALSE)
    blanket <- blanket.new(genes, val=0)
    blanket <- blanket.lift(blanket, mods, pairs, val=0.5)
    expect_equal(sum(blanket == 0.5), 5^2+1^2)
    pairs <- data.frame("M1", "M3", stringsAsFactors=FALSE)
    blanket <- blanket.new(genes)
    blanket <- blanket.lift(blanket, mods, pairs)
    expect_equal(sum(blanket == 0.5), 3^2+3^2)
    
    # Complexity reduction
    x <- blanket.cred(blanket)
    tri <- blanket[upper.tri(blanket, diag=F)]
    expect_equal(x, length(tri)/sum(tri != 0))
    
    # Inform blanket with prior
    prior <- blanket.new(genes, val=0)
    prior[TRUE] <- runif(length(genes)^2, 0, 1)
    diag(prior) <- 0
    
    mods <- list("M1"=c("G1", "G2"), "M2"=c("G3", "G4", "G5"), "M3"=c("G6"))
    blanket <- blanket.new(genes, val=0)
    blanket <- blanket.lift(blanket, mods, val=0.5)
    blanket.informed <- blanket.inform(blanket, prior)
    
    # Check
    expect_equal(blanket.informed[blanket == 0], rep(0, sum(blanket == 0)))
    expect_equal(blanket.informed[blanket != 0], prior[blanket != 0])
    
    # Parent blanket
    prior <- blanket.informed
    prior.blanket <- blanket
    
    # Child blanket
    mods <- list("M1"=c("G1", "G2"), "M2"=c("G3", "G4"), "M3"=c("G5", "G6"))
    blanket <- blanket.new(genes, val=0)
    blanket <- blanket.lift(blanket, mods, val=0.5)
    
    blanket.informed <- blanket.inform(blanket, prior, prior.blanket)
    expect_equal(blanket.informed[blanket == 0], rep(0, sum(blanket == 0)))
    expect_equal(blanket.informed[blanket != 0 & prior.blanket != 0], prior[blanket != 0 & prior.blanket != 0])
    expect_equal(blanket.informed[blanket != 0 & prior.blanket == 0], rep(0.5, sum(blanket != 0 & prior.blanket == 0)))
})
