outlierPercentage=0.01
gcCorrRes=0.001

# 
# function to compute average rd value for every gc value  
makegcBins <- function(bin,windSize,chr){  
  ##create vectors for gc and number of reads
  gcBin <- c()
  numReads <- c()
  numHits <- c()
  for(i in 1:((1/windSize)+1)){
    gcBin[[i]] <- ((i*windSize)-windSize)
    numReads[[i]] <- NA
    numHits[[i]] <- 0
  }
##fill the vectors
  for(i in 1:length(bin$gc)){
    if(!(is.na(bin$rd[[i]]))){
      val <- round(bin$gc[[i]]/windSize)+1
      if(!is.na(val)){
        ##initialize to zero the first time
        if(is.na(numReads[[val]])){
          numReads[[val]] <- 0
        }
        numReads[[val]] <- numReads[[val]] + bin$rd[[i]]
        numHits[[val]] <- numHits[[val]] + 1
      }
    }  
  }
  return(data.frame(gc=gcBin,avgreads=numReads/numHits,numreads=numReads))
}

# function to stripoutliers
stripOutliers <- function(nonNAbins,outlierPercentage){
	num=round(length(nonNAbins$avgreads)*outlierPercentage)
		top = sort(nonNAbins$avgreads)[length(nonNAbins$avgreads)-num]
		btm = sort(nonNAbins$avgreads)[1+num]
		rmlist = which(nonNAbins$avgreads > top)
		rmlist = append(rmlist,which(nonNAbins$avgreads < btm))
		for(i in 1:length(rmlist)){
			nonNAbins$avgreads[rmlist[i]] = NA
		}

	return(list(nonNAbins,rmlist))
}

# balanced center 
balancedCenter <- function(pos, numReads){
	center = 0
	netChange = sum((pos-center)*numReads)
	while(netChange > 0){
		center = center + 1
		netChange = sum((pos-center)*numReads)
	}
	margin = sum(numReads)*0.001
	adj = 0.5
	center = (center-1)+adj
	netChange = sum((pos-center)*numReads)
	while((abs(netChange) > margin) & (adj > 0)){
		adj = adj/2
		if(netChange > 0){
			center = center + adj
			}
		else{
			center = center - adj
			}
		netChange = sum((pos-center)*numReads)
	}
	return(center)
}

# returns stripped outliers 
returnOutliers <- function(results,rmlist,preOutliers){
  
  ##first, group adjacent windows that were removed
  counter = 1
  segs = vector("list")
  rmlist = sort(rmlist)
  segs[[1]] = rmlist[1]
    
  for(i in 2:length(rmlist)){    
    ##adjacent windows, merge
    if(rmlist[i]-1 == segs[[counter]][length(segs[[counter]])]){
      segs[[counter]] = append(segs[[counter]],rmlist[i])
    } else { #just add
      counter = counter + 1
      segs[[counter]] = rmlist[i]
    }
  }
  
  ##now, put them back in
  for(i in 1:length(segs)){
    ##edge case: beginning
    if(1 %in% segs[[i]]){
      
      ##can't avg, get the gc content of the closest non-outlier gc values (high)
      higc = preOutliers[(which(preOutliers$gc == preOutliers$gc[max(segs[[i]])])+1),]$gc
      ##get adjustment for that nearby window      
      val = results[which(results$gc==higc),]$adj
      
      for(j in 1:length(segs[[i]])){
        results = rbind(results,data.frame(adj=val,gc=preOutliers$gc[segs[[i]][j]]))
      }    
      
      ##other edge case: end
    } else if(max(segs[[i]]) > length(preOutliers$reads)){
      
      ##can't avg, get the gc content of the closest non-outlier gc values (low)
      logc = preOutliers[(which(preOutliers$gc == preOutliers$gc[min(segs[[i]])])-1),]$gc
      ##get adjustment for that nearby window
      val = results[which(results$gc==logc),]$adj
      
      for(j in 1:length(segs[[i]])){
        results = rbind(results,data.frame(adj=val,gc=preOutliers$gc[segs[[i]][j]]))
      }    
      
      ##in the middle, average
    } else{
      ##get the gc contents of the closest non-outlier gc values (high and low)
      logc = preOutliers[(which(preOutliers$gc == preOutliers$gc[min(segs[[i]])])-1),]$gc
      higc = preOutliers[(which(preOutliers$gc == preOutliers$gc[max(segs[[i]])])+1),]$gc
      
      ##now, average the adjustments for those two nearby windows
      val = mean(c(results[which(results$gc==higc),]$adj, results[which(results$gc==logc),]$adj))
      
      for(j in 1:length(segs[[i]])){
        results = rbind(results,data.frame(adj=val,gc=preOutliers$gc[segs[[i]][j]]))
      }    
    }
  }
  return(results)
}

# do gc correction based on the adjustment values obtained 
doCorrection <- function(bin, gcCorrRes, gcAdj, chr){
	adjustForGC <- function(x){
		if( (!(is.na(x["gc"]))) & (!(is.na(x["rd"])))){
	##don't correct windows with zero reads
			if(!(x["rd"]==0)){
				x["rd"] <- x["rd"] - gcAdj[round(x["gc"]/gcCorrRes)+1]
				if(x["rd"] < 0){
					x["rd"] = 1
				}
			}
		}
		return(x)
	}

## do the adjustment, then convert the matrix to a dataframe, and return
## the adjusted read depth in listified format
	z=data.frame(t(as.matrix(apply(bin,1,adjustForGC))))$rd
		return(z)
}
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
#--------------------------------------------------------------------------------------#
#construct gcBins datastructure as gcBins$gc gcBins$avg

# first step is to create gcBins 
# avgreads and numreads 
args <- commandArgs(trailingOnly=TRUE)
rdFile <- args[1]
gcFile <- args[2]
rdValues <- scan(rdFile,what=double())
gcValues <- scan(gcFile,what=double())
binInfo=data.frame(gc=gcValues,rd=rdValues);
gcBins = makegcBins(binInfo,gcCorrRes);

##strip out those without any reads
nonNAbins=na.omit(gcBins)
##remove outliers if specified
rmlist = NULL
preOutliers=NULL
if(outlierPercentage > 0){ 
	preOutliers=nonNAbins
	tmp = stripOutliers(nonNAbins, outlierPercentage)
	nonNAbins = tmp[[1]]
	rmlist = tmp[[2]]
}
nonNAbins = na.omit(nonNAbins)

# create loess object 
gc = nonNAbins$gc
reads = nonNAbins$avgreads
numreads = nonNAbins$numreads
reads.loess <- loess(reads ~ gc, span=0.75, data.frame(reads=reads, gc=gc))

# balance center
bmed=balancedCenter(reads.loess$fitted,numreads)  
reads.adj = reads.loess$fitted-bmed
reads.fix = reads-reads.adj
reads.resLoess = loess(reads.fix ~ gc, span=0.75, data.frame(reads.fix=reads.fix, gc=gc))
results = data.frame(adj=reads.adj,gc=gc)

# restore removed outliers
if(outlierPercentage > 0){
	results = returnOutliers(results,rmlist,preOutliers)
}

# restore NA bins we removed at the beginning
x = 1:length(gcBins$gc)
restoreNAs <- function(z){
	if(!(is.na(gcBins$avgreads[z]))){
		return(results[which(results$gc==gcBins$gc[z]),]$adj)
	} 
	else {
		return(NA)
	}
}
adjustments = sapply(x,restoreNAs)

# correct rd values
chr="chr1"; 
correctedRDValues =  doCorrection(binInfo, gcCorrRes, adjustments, chr)
			        
plotGC(rdValues,correctedRDValues,gcValues,"GC-plot")
# write the corrected values to the file 
write(correctedRDValues,file=rdFile,sep="\n");
