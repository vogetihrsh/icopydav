library(bamsignals)
library(GenomicRanges)


# read commandline arguments
args <- commandArgs(TRUE)
bamFile = args[1]
outputPrefix = args[2]
windowSize = as.integer(args[3])
chromosome = args[4]
chromosomeLength = as.double(args[5])

# construct the query
query = GRanges(chromosome,IRanges(breakInChunks(chromosomeLength,windowSize)))
counts = as.character(bamCount(bamFile,query))

# create bincor file
startPos = as.character(start(query)-1)
endPos = as.character(end(query))
fileConnection = file(paste(outputPrefix,"_preproc_intermediate.txt",sep=""))
writeLines(paste(startPos,endPos,counts,sep=" "),con=fileConnection)
close(fileConnection)

# create bincounts file


