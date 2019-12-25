
<!-- README.md is generated from README.Rmd. Please edit that file -->

# shine

[![](https://img.shields.io/badge/platforms-linux%20%7C%20osx%20-2a89a1.svg)]()
[![](https://img.shields.io/badge/lifecycle-stable-4ba598.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![](https://img.shields.io/github/last-commit/montilab/shine.svg)](https://github.com/montilab/shine/commits/master)

Bayesian **S**tructure Learning for **Hi**erarchical **Ne**tworks  
*A package to aid in the Bayesian estimation of hierarchical biological
regulatory networks via Gaussian graphical modeling*

## Documentation

Please visit <https://montilab.github.io/shine/> for comprehensive
documentation.

## Requirements

We suggest R 3.6.0 but R (\>= 3.5.0) is required to install directly
from Github. For workflows, Nextflow can be used on any POSIX compatible
system (Linux, OS X, etc) and requires BASH and Java 8 (or higher) to be
installed. Alternatively, check out usage with Docker.

## Installation

Install the development version of the package from Github.

``` r
devtools::install_github("montilab/shine")
```

## Basic Overview

``` r
library(shine)
```

### Example Data

``` r
# Toy expression set object
data(toy)
```

``` r
dim(toy)
```

    #> Features  Samples 
    #>      150       30

``` r
table(toy$subtype)
```

    #> 
    #>  A  B  C 
    #> 10 10 10

### Variable Selection

``` r
# Filter out non-varying genes
genes.filtered <- filter.var(toy, column="subtype", subtypes=c("A", "B", "C"))

# Select top genes by median absolute deviation
genes.selected <- select.var(toy, column="subtype", subtypes=c("A", "B", "C"), genes=genes.filtered, limit=150)

# Subset toy dataset
eset <- toy[genes.selected,]
```

### Structural Constraints

``` r
# Find constraints for structure with zero-order correlations
mods <- mods.get(eset, min.size=5, cores=3, do.plot=FALSE)
meta <- metanet.build(mods$eigengenes, cut=0.5, mpl=TRUE, iter=20000, cores=3)

# Set constraints
blanket <- blanket.new(mods$genes)
blanket <- blanket.add.mods(blanket, mods$mods)
blanket <- blanket.add.modpairs(blanket, mods$mods, meta$metanet.edges)
```

### Hierarchical Network Workflows

``` r
condition <- "subtype"

# Learn networks as a hierarchy
hierarchy <- "
A_B_C -> A_B
A_B_C -> C
A_B -> A
A_B -> B
"

# Generate workflow
build.workflow(hierarchy,
               condition,
               eset,
               blanket,
               iter=10000,
               cores=3)
```

``` bash
./nextflow workflow.nf -c workflow.config -profile local
```

``` md
N E X T F L O W  ~  version 19.10.0
Launching `workflow.nf` [mad_northcutt] - revision: e7aff2bf48
-
W O R K F L O W ~ Configuration
===============================
data      : /Users/anthonyfederico/Downloads/nextflow/data/data.rds
output    : /Users/anthonyfederico/Downloads/nextflow
-------------------------------

Hierarchy
A_B_C -> A_B
A_B_C -> C
A_B -> A
A_B -> B

-

executor >  local (5)
[4d/715ff1] process > A_B_C [100%] 1 of 1 ✔
[87/c2d03b] process > C     [100%] 1 of 1 ✔
[6b/544fdd] process > A_B   [100%] 1 of 1 ✔
[96/916564] process > A     [100%] 1 of 1 ✔
[43/c599d4] process > B     [100%] 1 of 1 ✔

Completed at: 24-Dec-2019 23:01:05
Duration    : 1m 45s
CPU hours   : 0.1
Succeeded   : 5
```

### Network Visualization

Check out <https://github.com/montilab/netviz> to explore your networks.
