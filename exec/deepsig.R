#!/usr/bin/env Rscript
library('argparse')
library('DeepSig')
parser <- ArgumentParser(description = 'Extract mutational signatures using the DeepSig algorithm')
parser$add_argument('-i', '--input', dest = 'catalog', metavar = 'CATALOG', required = TRUE,
                    help = 'input catalog data file')
parser$add_argument('-o', '--output', dest = 'output', metavar = 'OUTPUT', required = TRUE,
                    help = 'output directory')
parser$add_argument('-c', '--cancer-type', dest = 'cancer', metavar = 'CANCER_TYPE',
                    required = TRUE, help = 'use a model for the specified cancer type')
parser$add_argument('-q', '--quiet', dest = 'quiet', action = 'store_true', default=FALSE,
                    help = 'Run quietly')
args <- parser$parse_args()

# Set verbosity
verbose <- 1
if(args$quiet) {
  verbose <- 0
}

# Read input
data <- read.table(args$catalog, header=TRUE, sep = '\t', check.names = FALSE)
# Extract signatures
output <- DL.call(catalog = t(data),
                  cancer.type = args$cancer,
                  verbose = verbose)
# Write each table to a file in output dir
sapply(names(output),
 function (x) write.table(output[[x]], file=file.path(args$output, x))
)