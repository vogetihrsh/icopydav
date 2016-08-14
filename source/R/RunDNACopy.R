

library("DNAcopy");
library("parallel");
library("ParDNAcopy");

args <- commandArgs(trailingOnly = TRUE)
inputFile=args[1];
numOfDataPoints=args[2];
numOfProcesses=strtoi(args[3]);

dataPoint <- scan(inputFile,what=double());

chrom <- rep(3,strtoi(numOfDataPoints));

maploc=1:strtoi(numOfDataPoints);
ranseed=sample(1:123,1);
obj=CNA(dataPoint,chrom,maploc,data.type="logratio",presorted=TRUE);
segmented_output<-parSegment(obj,ranseed=ranseed,distrib="Rparallel",njobs=numOfProcesses,undo.splits="sdundo");
segmented_output
