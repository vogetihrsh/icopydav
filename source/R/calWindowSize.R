

calWindParams <- function(numReads,fdr,genomeSize,oDisp, ploidy, minSize, ploidyPerc, medAdj, startDiv=1000,inputPloidy){#nullThis, nullBoth, startDiv=100){
	## answer has to be this close to the FDR rate (on the lower side)
	tolerance = 0.005
	#starting point for search
	windSize = genomeSize/startDiv
	divider = NULL

	## first, halve size until we get above FDR threshold
	found <- FALSE
	p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj, inputPloidy)
	
	if(p$fdr > fdr){
			stop("not enough reads to achieve the desired FDR rate")
		}

	while((found == FALSE) & (windSize > minSize)){
		windSize <- round(windSize / 2)
		p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj, inputPloidy)
		if(p$fdr > fdr){
			found = TRUE
			}
	}
	
	## zero in on a size that's within the desired parameters
	if(windSize > minSize){
		found = FALSE    
		adj = windSize/2
		while((found == FALSE) & (windSize > minSize)){    
			if(p$fdr > fdr){
				windSize <- round(windSize + adj)
				p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,inputPloidy=inputPloidy)
			}
			else if(p$fdr < (fdr - tolerance)){
				windSize <- round(windSize - adj)
				p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,inputPloidy=inputPloidy)

			} 
			else{
				found = TRUE
			}
				adj <- adj/2
			}
	}

	##if the window size is below the minimum size, have to
	##recalculate params
	if(windSize < minSize){
		windSize = minSize
	} else {
		windSize = floor(windSize/minSize)*minSize # rounding to multiplies of min size
	}
	p <- fdrRate(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,inputPloidy=inputPloidy)
	med=p$med
	div=p$div
	return(data.frame(binSize=windSize,div=div,med=med))
}


# fdrRate function calculates and returns the FDR value for a given window size
# use above normal and below normal 
fdrRate <- function(windSize,ploidy,genomeSize,numReads,oDisp,ploidyPerc, medAdj,inputPloidy){ #nullThis,nullBoth){
	numWinds <- genomeSize/windSize
	med <- (numReads/(((genomeSize*ploidyPerc$belowNormalPerc*(inputPloidy-1)/inputPloidy) +(genomeSize*ploidyPerc$aboveNormalPerc*(inputPloidy+1)/inputPloidy) +(genomeSize*ploidyPerc$normalPerc)) / windSize))*medAdj
	medAlt <- med*(ploidy/2);
	divFdr=NULL;
	divFdr = dividePeaks(med, medAlt, numWinds*ploidyPerc$normalPerc, numWinds*ploidyPerc$aboveNormalPerc, oDisp)
	return(data.frame(fdr=divFdr$fdr, div=divFdr$div, med=med))
}

# function that calculates number fdr given a window length
# ths function is called by FDR function
dividePeaks <- function(amed,bmed,aNum,bNum,oDisp){  
	thresholds =(amed+1):(bmed-1)
	mislabeledWinds <- function(thresh){
		low=(1-pnbinom(thresh,size=(amed/oDisp-1),mu=amed))*aNum
		high=(pnbinom(thresh,size=(bmed/oDisp-1),mu=bmed))*bNum
		return(low+high)
		}

	vals=sapply(thresholds,mislabeledWinds)
	lows = which(vals == min(vals))
	
	##choose the low point closest
	##to the halfway point
	diffFromMed = abs((lows+amed)-(amed+bmed/2))
	pos = which(diffFromMed == min(diffFromMed))
	if(length(pos) > 1){
		pos = pos[1]
		}

	return(data.frame(div=lows[pos]+amed,fdr=mislabeledWinds(lows[pos]+amed)/(aNum+bNum)))
}

ploidyPercentages <- function(effectiveGenomeSize,percCNGain,percCNLoss){
	  ##first, get the coverage that come from diploid chrs
	  normalPerc=1-percCNLoss - percCNGain;
	  return(data.frame(belowNormalPerc=percCNLoss,
					      normalPerc=normalPerc,
					      aboveNormalPerc=percCNGain))
}

									   

ploidy=3;
medAdj=1;
args <- commandArgs(trailingOnly = TRUE);
numReads=strtoi(args[1]);
genomeSize=strtoi(args[2]);
fdr= as.numeric(args[3]);
minSize= strtoi(args[4]);
percCNLoss= as.numeric(args[5]);
percCNGain= as.numeric(args[6]);
overDispersion= as.numeric(args[7]);
inputPloidy=strtoi(args[8]);
ploidyPerc=ploidyPercentages(genomeSize,percCNGain,percCNLoss);
pTrip <- calWindParams(numReads=numReads, fdr=fdr, genomeSize=genomeSize, oDisp=overDispersion,ploidy=ploidy,
		minSize=minSize, ploidyPerc=ploidyPerc,medAdj=medAdj,startDiv=1000,inputPloidy=inputPloidy );
pTrip$binSize;
