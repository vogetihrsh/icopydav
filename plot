#!/bin/bash

cFlag=0;
sFlag=0;
eFlag=0;

SOURCEDIR=$(dirname $(readlink -f $0));
inputFlag=0;

SHORTOPTS="i:do:"
LONGOPTS="hg18,hg19,genome:,help,start:,end:"
PROGNAME="plot_c.sh"

ARGS=$(getopt -s bash --options $SHORTOPTS  --longoptions $LONGOPTS --name $PROGNAME -- "$@" )

eval set -- "$ARGS"

while true; do
	case "$1" in
		-i) inputFile="$2"; inputFlag=1; shift 2;;
		-o) outputFile="$2"; outputFlag=1; shift 2;;
		-d) cFlag=1; shift ;;
		--start) sTart="$2"; sFlag=1; shift 2; ;;
		--end) end="$2"; eFlag=1; shift 2; ;;
		--hg18) genome="hg18"; genomeFlag=1; shift ;;
		--hg19) genome="hg19"; genomeFlag=1; shift ;;
		--genome) genome="$2"; genomeFlag=1; shift 2;;
		--help) echo "Usage: plot -i <input bed> -o <output prefix> <genome flag>" ;
			exit 1; ;;
		--) shift; break;;		
		*) echo "Internal error!";
	esac
done

if [[ $inputFlag -eq 0 || $genomeFlag -eq 0 || $outputFlag -eq 0 ]]; then
	echo "One or more required paramters is/are missing">&2;
		exit 1;
elif [[ ! -f $inputFile ]]; then
	echo "Error: $inputFile not found!" >&2;
		exit 1;
fi

if [[ $cFlag -eq 1 ]]; then
	if [[ $sFlag -eq 1 || $eFlag -eq 1 ]]; then
		sFlag=1;
		eFlag=1;
	else
		echo "Please enter start/end coordinates" >&2;
	        exit 1;
	fi
else
	echo "Generating plot for complete chromomsome, please enter coordinates with parameters '-d', '--start' and '--end' for specific genomic region visualization!";
	echo "-------------------------------------------------------------------";
fi

rm -rf $outputFile 
echo "Chromosome	Start	End	Call" >> $outputFile
awk -v OFS='\t' '{print $1,$2,$3,$5}' $inputFile >> $outputFile

Rscript --slave $SOURCEDIR"/source/R/plotCNV.R" $outputFile $genome $sTart $end;

rm -rf $outputFile
