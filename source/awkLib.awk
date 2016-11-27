# currentState has start, end, gap and copynumber
# currentCNV has start, end and copynumber


# returns an empty object of currentState
function intializeCurrentState(currentState){
	currentState["start"]=currentState["end"]=-1;
	currentState["copynumber"] = currentState["gap"]=0;
	
}

# rounding off function 
function roundOff(i){
	ret=sprintf("%.0f",i);
	return ret;
}

# add CNVs to the current state
function updateCurrentState(currentState, currentCNV){
	
	# if currentState is empty assign currentCNV and return 
	if(currentState["start"]==-1){
		currentState["start"]=currentCNV["start"];
		currentState["end"] = currentCNV["end"];
		currentState["copynumber"] = currentCNV["copynumber"];
		currerntState["gap"] = 0;
		return ;
	}

	coverage = (currentState["end"] - currentState["start"] - currentState["gap"])*currentState["copynumber"];
	coverage += (currentCNV["end"] - currentCNV["start"])*currentCNV["copynumber"];
	regionSize = (currentState["end"]-currentState["start"]-currentState["gap"]) + (currentCNV["end"] - currentCNV["start"]);
	# update currentState variables
	currentState["gap"] = currentState["gap"] + currentCNV["start"] - currentState["end"];
 	currentState["end"] = currentCNV["end"];
	currentState["copynumber"] = coverage/regionSize ;

}
# checks if two CNs are same
function isCopyNumberSame(c1,c2){
	
	i=roundOff(c1);
	j=roundOff(c2);
	if(i==j || absolute(c1-c2)<0.4)
		return 1;
	else
		return 0;
}

# checks the merging ratio criteria
function checkCriteria(currentState, currentCNV){
	
	currentGap = currentCNV["start"]-currentState["end"];

	# special case for copy number zero
        if(currentState["copynumber"]< 0.6 && currentCNV["copynumber"] < 0.6 && currentGap<5000)
		return 1;
	totalGap=currentState["gap"]+ currentGap;
	totalCNVRegion = (currentState["end"] - currentState["start"]- currentState["gap"]) + (currentCNV["end"] - currentCNV["start"]);
	
	ratio = totalGap/(totalCNVRegion + totalGap);
	if(ratio < mergeFraction)
		return 1;
	return 0;
}

# merges cnvs from startIndex to stopIndex
function mergeCNV(mergedCNV,cnv,startIndex,stopIndex,size){
	if(stopIndex>size)
		stopIndex=size;

	mergedCNV[1]=cnv[startIndex,1];
	mergedCNV[2]=cnv[stopIndex,2];
	mergedCNV[3]=cnv[startIndex,3];
	
	x1=0;
	x2=0;
	i=startIndex;

	while(i<=stopIndex){
		csize=cnv[i,2]-cnv[i,1];
		x2+=csize;
		x1+=csize*cnv[i,4];
		i++;
	}
	mergedCNV[4]=x1/x2;
}

# returns curretnCNV object given CNV list and index
function getCurrentCNV(currentCNV,cnv,ind){
	currentCNV["start"] = cnv[ind,1];
	currentCNV["end"] = cnv[ind,2];
	currentCNV["copynumber"] = cnv[ind,4];

}

# this is the main fucntion
# should be used to understand the overflow of the program

function merge(finalList,cnv,size){
	fp=1;
	currentIndex=1;
	intializeCurrentState(currentState);
	while(currentIndex<=size){
		
		getCurrentCNV(currentCNV,cnv,currentIndex);
		maxIndex=startIndex=currentIndex;
		updateCurrentState(currentState, currentCNV);
		currentIndex++;
		getCurrentCNV(currentCNV,cnv,currentIndex);
		
		while(currentIndex<=size && isCopyNumberSame(currentState["copynumber"], currentCNV["copynumber"])==1){
			if(checkCriteria(currentState,currentCNV)==1)
				maxIndex=currentIndex;
			updateCurrentState(currentState, currentCNV);
			currentIndex++;
			getCurrentCNV(currentCNV,cnv,currentIndex);
		}
		
		# get merged CNV and append it to the list 
		mergeCNV(temp,cnv,startIndex,maxIndex,size);
		finalList[fp,1]=temp[1]
		finalList[fp,2]=temp[2]
		finalList[fp,3]=temp[3]
		finalList[fp,4]=temp[4]
		fp++;


		intializeCurrentState(currentState);
		currentIndex = maxIndex+1;
	}
}

# simple absolute function 
function absolute(n){
	if(n<0)
		return -n;
	else
		return n;
}

BEGIN{
		if(mergeFraction<0|| mergeFraction >1)
			mergeFraction=0.2;
		entryPointer=-1;
		prev=-1;
		cnvType=-1;
		curAvg=0;
		count=0;

		prevRD=-1;
		prevEndPoint=-1;
		prevStart=-1;

		PTR=0;
		diff=0;
	}

	{
		diff=absolute(prevRD-$3);
		if(entryPointer==1&&($3<=lc||$3>=hc)&&(diff<avg/2)&&prevEndPoint==$1)
		{
			curAvg=curAvg+$3;
			count=count+1;
		}

		else if(entryPointer==1)
		{
			entryPointer=-1;
			epoint=prevEndPoint;
			curAvg=curAvg/count;
			cn=(curAvg*2)/avg;
			if(count>1){
				PTR=PTR+1;
				CNVARR[PTR,1]=spoint;CNVARR[PTR,2]=epoint;CNVARR[PTR,3]=cnvType;CNVARR[PTR,4]=cn;
			}
			cnvType=-1;
			curAvg=0;
			count=0;
		}

		if(entryPointer==-1 && ($3<=lc || $3>=hc) )
		{
			entryPointer=1;
			spoint=$1;

			if($3<=lc)
				cnvType=0;
			else
				cnvType=1;
			count = count +1;
			curAvg = curAvg + $3;
		}
		prevEndPoint = $2;
		prevStart=$1;
		prevRD=$3;
	}

	END{
		if(entryPointer==1 && count > 1){
			epoint=prevEndPoint;
			curAvg=curAvg/count;
			cn=(curAvg*2)/avg;
			PTR=PTR+1;
			CNVARR[PTR,1]=spoint;CNVARR[PTR,2]=epoint;CNVARR[PTR,3]=cnvType;CNVARR[PTR,4]=cn;
		}
#	for(i=1;i<=PTR;i++)
#			print CNVARR[i,1],CNVARR[i,2],CNVARR[i,3],CNVARR[i,4];
		merge(CNVlist,CNVARR,PTR);
		size=length(CNVlist)/4;
		for(i=1;i<=size;i++)
			print  chrName,CNVlist[i,1],CNVlist[i,2], CNVlist[i,3], CNVlist[i,4] >> outfile
	}
