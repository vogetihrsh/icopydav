
## do the mapability correction

# File input
args <- commandArgs(trailingOnly=TRUE)
map= args[1]
mapVec=scan(map)
window_size= as.numeric(args[2])
len = as.numeric(args[3])
output = as.character(args[4])

## Script 
a <- 1:ceiling(len/window_size)
numPerBin=window_size/100
getMapMean <- function(x){    
	wind = mapVec[(((x-1)*(numPerBin))+1):(x*(numPerBin))]
    	mn = mean(wind, na.rm=TRUE)
    	if(is.nan(mn)){
      		return(NA)
    	} else{
      		return(mn)
    	}
}
map_new = sapply(a,getMapMean)
write.table(as.data.frame(map_new),file=paste(output,"_",window_size,".map", sep=""),quote=F,col.names=F,row.names=F)
