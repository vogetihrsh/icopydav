library(GenomicAlignments)
library(rtracklayer)
args <- commandArgs(trailingOnly=TRUE)
bamFile <- args[1]
binFile <- args[2]

#Import Bam file
param <- ScanBamParam(mapqFilter=1)
bamfile = readGAlignments(bamFile,param=param)

#Import window bed file
bedfile = import.bed(binFile)

#count numbe rof overlaps for each element of the bedfile object
overlap.counts <- countOverlaps(bedfile,bamfile, type="any")
col = cat(overlap.counts,sep="\n")
col

