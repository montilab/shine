#!/bin/sh

Rscript $run \\
"${data}" \\
"${prior}" \\
"${label}" \\
\$(echo ${include});
