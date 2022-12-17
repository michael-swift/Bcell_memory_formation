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

parser <- add_option(parser, c("-l", "--labels"), default = "celltypist")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parse_args(parser)


INPUT_COUNTS = args$c
INPUT_ANNO = args$a
INPUT_FEAT = args$f
OUTPUT_DIR = args$o
SAMPLENAME = args$n
MODE = args$m
library("singleCellTK", quietly=TRUE)

message(paste0("Processing sample ", SAMPLENAME, "..."))
  
sce <-singleCellTK::importFromFiles(assayFile = INPUT_COUNTS, annotFile = INPUT_ANNO, featureFile = INPUT_FEAT, annotFileHeader = True, annotFileSep = ',', featureSep = ",", featureHeader = True)

#remove null cells
rm.ix <- which(colSums(assay(sce, 'counts')) == 0)
if(length(rm.ix) > 0){
  sce <- sce[,-rm.ix]
}

# run decontX
sce <- singleCellTK::runDecontX(sce, useAssay = "counts", z = "celltypist")
message("Exporting...")
singleCellTK::exportSCEtoAnnData(sce,
                                 useAssay = "decontXcounts", 
                                 outputDir = OUTPUT_DIR, 
                                 prefix= SAMPLENAME)
message("Done!")