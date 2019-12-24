#!/bin/sh

SCRIPT="cat(system.file('scripts/run.R', package='shine'))";
RUN=\$(Rscript -e "\$SCRIPT");

Rscript \$RUN \\
"${params.path.eset}" \\
"${params.path.genes}" \\
"${prior}" \\
"${params.path.blanket}" \\
"${params.iters}" \\
"${params.cores}" \\
"${params.condition}" \\
"${label}" \\
\$(echo ${include});
