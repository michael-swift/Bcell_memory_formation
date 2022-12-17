#!/usr/bin/env Rscript
library("optparse", quietly=TRUE)

# create parser object
parser <- OptionParser()
parser <- add_option(parser, c("-i", "--input"), default="./",
                     help="path to input directory")
parser <- add_option(parser,c("-o", "--output"), default="./",
                    help="path to output directory")
parser <- add_option(parser, c("-n", "--name"), default="sample")

parser <- add_option(parser, c("-m", "--mode"), default="raw")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parse_args(parser)

SAMPLENAME = args$n
INPUT_DIR = args$i
OUTPUT_DIR = args$o
MODE = args$m
library("singleCellTK", quietly=TRUE)

message(paste0("Processing sample ", SAMPLENAME, "..."))
  
sce <-singleCellTK::importCellRanger(sampleDirs = INPUT_DIR, dataType = MODE)

#remove null cells
rm.ix <- which(colSums(assay(sce, 'counts')) == 0)
if(length(rm.ix) > 0){
  sce <- sce[,-rm.ix]
}
# run decontX
sce <- singleCellTK::runDecontX(sce)
message("Exporting...")
singleCellTK::exportSCEtoAnnData(sce,
                                 useAssay = "decontXcounts", 
                                 outputDir = OUTPUT_DIR, 
                                 prefix= SAMPLENAME)
message("Done!")