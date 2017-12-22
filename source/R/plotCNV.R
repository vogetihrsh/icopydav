## CNV visualization
## This script assumes to have a header with labels: Chromosome \t Start \t End \t Call
## Ploidy is assumed to be 2


args <- commandArgs(trailingOnly = TRUE);
file=args[1];
ref_gen=args[2];
u_start=args[3];
u_end=args[4];
user_def=c(u_start,u_end)
#print(ref_gen)
x_axis=0

dataTable <-read.table(file, header=TRUE);
ratio<-data.frame(dataTable)
ploidy <- 2 # amplification > 2, deletion < 2

library(quantsmooth)
chrom=ratio[,1]
chr = gsub("^chr","",chrom) # removes 'chr' from col 1

ratio1 = cbind(chr,ratio[,2:4])
ratio = ratio1

png(filename = paste(file,".png",sep = ""), width =1268, height =720,units = "px", pointsize=20,bg = "white", res = 72)
par(mfrow=c(2,1))
for (i in (chr)) {
	region <- which(ratio$chr==i)
	chrom = unique(ratio$Chromosome)
        if (length(region)>0) {
		xscale = c(0,(lengthChromosome(i,"bases")))
		yscale=c(0,0.5)
		ideo.width<-0.2
		ideo.ypos<-ideo.width
		ideo.bleach<-0.0
		plot(xscale,yscale,bty="n",type="n",main="",xlab= "",ylab="",xaxt="n", yaxt="n")
		paintCytobands(i,pos=c(0,ideo.ypos),units=ref_gen,width=ideo.width,legend=TRUE,bleach=ideo.bleach)
		region <- which(ratio$chr==i)
		if (is.na(u_start)){
			ratio$Start <- ratio$Start
			ratio$End <- ratio$End
		} else	{
			ratio$Start[ratio$Start < as.numeric(u_start)] <- NA
			ratio$End[ratio$Start < as.numeric(u_start)] <-NA
			ratio$Start[ratio$Start > as.numeric(u_end)] <- NA
			ratio$End[ratio$Start > as.numeric(u_end)] <- NA
		}
		region <- which(ratio$chr==i  & ratio$Call>ploidy )
		pos = seq(0.10, 0.10, length.out=length(region))
		points(ratio$Start[region],pos,pch="|",cex=3,col = "red")
		region <- which(ratio$chr==i  & ratio$Call<ploidy )
		pos1 = seq(0.10, 0.10, length.out=length(region))
		points(ratio$Start[region],pos1,pch="|",cex=3,col = "blue")
	
	}

	if (length(region)>0) {
		if (is.na(u_start)){
			x_axis <- xscale
		} else {
			x_axis <- as.numeric(user_def)
		}
		ratio$Call[ratio$Call > 10] <- 10
		plot(ratio$Start[region],ratio$Call[region],xlim = x_axis,ylim = c(0,max(ratio$Call)),xlab = paste ("Chromosome", i, "(position)"),ylab = "Absolute Copy Number",pch = "-",col = "black")
		region <- which(ratio$chr==i  & ratio$Call>ploidy )
		segments(ratio$Start[region],ratio$Call[region],ratio$End[region],ratio$Call[region],col = "red",pch='-',lwd=6)
		region <- which(ratio$chr==i  & ratio$Call<ploidy )
		segments(ratio$Start[region],ratio$Call[region],ratio$End[region],ratio$Call[region],col = "blue",pch='-',lwd=6)
		region <- which(ratio$chr==i  & ratio$Call==ploidy)
		segments(ratio$Start[region],ratio$Call[region],ratio$End[region],ratio$Call[region],col = "darkgreen",lwd=6)
   }
   #dev.off()
}
##
dev.off()
