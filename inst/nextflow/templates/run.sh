#!/bin/sh

Rscript $PWD/scripts/run.R \\
"${params.path.eset}" \\
"${params.path.genes}" \\
"${prior}" \\
"${params.path.blanket}" \\
"${params.iters}" \\
"${params.cores}" \\
"${params.condition}" \\
"${label}" \\
\$(echo ${include});