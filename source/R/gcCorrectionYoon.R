# GC-cORRECTION USING YOON ET AL METHOD
# CURRENTLY DIVIDES THE GC CONTENT DOMAIN INTO 1000 PARTS 
# gc correction code input: read-depth gc-content(in fractions)

args <- commandArgs(trailingOnly = TRUE)
x=scan(file=args[1],what=double()) # intial scanning 
n=length(x)								
n1=n/2;					# number of reads
y=matrix(x,ncol=2,byrow=T) # y[n1,2] accessing elements 2d matrix containing read depth and gc content (
rm(x)
for (i in 1:n1)
{
	z=y[i,2];
	z=1000*z;
	z=round(z);
	if(z==0)
	   z=1;
	y[i,2]=z;
	
}
med=median(y[,1]); # median of all the RD

a=vector(mode="list",length=1000) # temporary remove after calculation of medians
for( i in 1:n1)
{
	z=y[i,2];

	if (z>0)
	{
		a[[z]]=c(a[[z]],y[i,1]);
	}
}
mgc=numeric()	# contains medians for a given a gc content
for(i in 1:1000)
{
	mgc[i]=0
	if (length(a[[i]])>0)
	{
		mgc[i]=median(a[[i]])

	}
}
rm(a);		# delete  
for(i in 1:n1)
{
	z=y[i,1];
	j=y[i,2]
	if(j>0)
	{	
		z=(z*med)/mgc[j];
		y[i,1]=z;
	}
	
}
for(i in 1:n1)
{
	cat(y[i,1]);cat("\n");
}

	
