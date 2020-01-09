FROM rocker/tidyverse:3.6.0

MAINTAINER Anthony Federico <anfed@bu.edu>

RUN Rscript -e \
    'install.packages("devtools"); \
     install.packages("BiocManager"); \
     BiocManager::install("impute"); \
     BiocManager::install("preprocessCore"); \
     BiocManager::install("GO.db"); \
     BiocManager::install("AnnotationDbi"); \
     BiocManager::install("Biobase"); \
     devtools::install_github("montilab/shine");'
