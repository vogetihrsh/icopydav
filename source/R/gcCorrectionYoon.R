# GC-cORRECTION USING YOON ET AL METHOD
# CURRENTLY DIVIDES THE GC CONTENT DOMAIN INTO 1000 PARTS 
# gc correction code input: read-depth gc-content(in fractions)

plotGC <- function(uncorrectedRD,correctedRD,gcContent,fileName){

# Remove all NA values  
	uncorrectedRD = uncorrectedRD[!is.na(gcContent)];
	correctedRD = correctedRD[!is.na(gcContent)];
	gcContent=gcContent[!is.na(gcContent)];
	
#
	medURD = median(uncorrectedRD);
	medRD = median(correctedRD);

# write the plots to a PNG file
	png(filename = paste(fileName,".png",sep = ""), width =1920, height =1080,units = "px", pointsize=20,bg = "white", res = 72);

# start plotting now 
	par(mfrow = c(2,1));
	plot(gcContent,uncorrectedRD,main = "Biased RD vs GC-content",col="violet",ylim=c(0,3*medURD),pch='.',xlab="GC fraction",ylab="Read Depth");
	lines(lowess(gcContent,uncorrectedRD),col="green",lwd=2);
	plot(gcContent,correctedRD,main = "Corrected RD vs GC-content",col="violet",ylim=c(0,3*medRD),pch='.',xlab="GC fraction",ylab="Read Depth");
	lines(lowess(gcContent,correctedRD),col="green",lwd=2);

}

args <- commandArgs(trailingOnly = TRUE)
x=scan(file=args[1],what=double()) # intial scanning 
type=args[2]
n=length(x)								
n1=n/2;					# number of reads
y=matrix(x,ncol=2,byrow=T) # y[n1,2] accessing elements 2d matrix containing read depth and gc content (
rm(x)

uncorrectedRD = y[,1];
gcContent = y[,2];

# remove all 0 RD value bins
nonZeroIndices = which(uncorrectedRD!=0);
nonZeroRD = uncorrectedRD[nonZeroIndices];
nonZeroRDGCContent= gcContent[nonZeroIndices];
med = median(nonZeroRD);

uniqueGCValues = sort(unique(nonZeroRDGCContent));
num = length(uniqueGCValues);
for (i in 1:num){
	indices  = which(nonZeroRDGCContent==uniqueGCValues[i]);
	counts = nonZeroRD[indices];
	gcMedian = median(counts);
	if(gcMedian==0)
		gcMedian=1;
	counts = (counts*med)/gcMedian;
	nonZeroRD[indices] = counts;
}
correctedRD = uncorrectedRD;
correctedRD[nonZeroIndices]=nonZeroRD;
for(i in 1:n1)
{
	cat(correctedRD[i]);cat("\n");
}

if(type=="gc")
	plotGC(correctedRD,uncorrectedRD,gcContent,"GC-plot")
