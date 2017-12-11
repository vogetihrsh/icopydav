
## GC score file preparation

# File input
args <- commandArgs(trailingOnly=TRUE)
gc= args[1]
gcVec=scan(gc,what=0,quiet=TRUE)
window_size= as.numeric(args[2])
len = as.numeric(args[3])
output = as.character(args[4])

## Script 
a <- 1:ceiling(len/window_size)
numPerBin=window_size/100
getGCMean <- function(x){    
	wind = gcVec[(((x-1)*(numPerBin))+1):(x*(numPerBin))]
    	mn = mean(wind, na.rm=TRUE)
    	if(is.nan(mn)){
      		return(NA)
    	} else{
      		return(mn)
    	}
}
gc_new = sapply(a,getGCMean)
write.table(as.data.frame(gc_new),file=paste(output,"_",window_size,".gc", sep=""),quote=F,col.names=F,row.names=F)
