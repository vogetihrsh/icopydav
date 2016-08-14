## CNV visualization
## This script assumes to have a header with labels: Chromosome \t Start \t End \t Call
## Ploidy is assumed to be 2
## Pls change 'xlim' for user defined coordinates, else It will create for complete chromosome

args <- commandArgs(trailingOnly = TRUE);
file=args[1];
dataTable <-read.table(file, header=TRUE);
ratio<-data.frame(dataTable)
ploidy <- 2 # amplification > 2, deletion < 2

##

png(filename = paste(file,".png",sep = ""), width =1268, height =720,
units = "px", pointsize=20,bg = "white", res = 72)
plot(1:10)

##In case of many chromosome data
#op <- par(mfrow = c(5,5))
##

chrom = unique(ratio$Chromosome)
for (i in (chrom)) {
    region <- which(ratio$Chromosome==i)

   if (length(region)>0) {
    plot(ratio$Start[region],ratio$Call[region],xlim = c(0,max(ratio$End[region])),ylim = c(0,max(ratio$Call)),xlab = paste ("Chromosome", i, "(position)"),ylab = "CNV",pch = ".",col = "black")
    region <- which(ratio$Chromosome==i  & ratio$Call>ploidy )
    segments(ratio$Start[region],ratio$Call[region],ratio$End[region],ratio$Call[region],col = "red",lwd=2)
    region <- which(ratio$Chromosome==i  & ratio$Call<ploidy )
    segments(ratio$Start[region],ratio$Call[region],ratio$End[region],ratio$Call[region],col = "blue",lwd=2)
    region <- which(ratio$Chromosome==i  & ratio$Call==ploidy)
    segments(ratio$Start[region],ratio$Call[region],ratio$End[region],ratio$Call[region],col = "darkgreen",lwd=2)
   }
   #dev.off()
}
##
dev.off()
