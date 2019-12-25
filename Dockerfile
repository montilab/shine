FROM rocker/tidyverse:3.6.0

MAINTAINER Anthony Federico <anfed@bu.edu>

RUN R -e 'install.packages("BiocManager");\
BiocManager::install("Biobase");\
BiocManager::install("impute");\
BiocManager::install("preprocessCore");\
BiocManager::install("GO.db");\
BiocManager::install("AnnotationDbi");\
devtools::install_github("montilab/shine");'
