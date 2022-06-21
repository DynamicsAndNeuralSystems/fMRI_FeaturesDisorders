# Parse arguments
library(argparse)
parser <- ArgumentParser(description = "Define data paths and feature set")

parser$add_argument("--github_dir", default="/project/hctsa/annie/github/")
parser$add_argument("--rdata_path", default="/project/hctsa/annie/data/scz/UCLA/Rdata/")
parser$add_argument("--pydata_path", default="/project/hctsa/annie/data/scz/UCLA/pydata/")
parser$add_argument("--feature_set", default="pyspi_19")

# Parse input arguments
args <- parser$parse_args()
github_dir <- args$github_dir
rdata_path <- args$rdata_path
pydata_path <- args$pydata_path
feature_set <- args$feature_set

