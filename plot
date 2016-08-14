#!/bin/bash

SOURCEDIR=$(dirname $(readlink -f $0));

rm -rf $2
echo "Chromosome	Start	End	Call" >> $2
awk -v OFS='\t' '{print $1,$2,$3,$5}' $1 >> $2
Rscript --slave $SOURCEDIR"/source/R/plotCNV.R" $2;
rm -rf $2
