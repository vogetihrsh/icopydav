args <- commandArgs(trailingOnly=TRUE)
rdFile <- args[1]
gcFile <- args[2]

rdValues <- scan(rdFile,what=double())
gcValues <- scan(gcFile,what=double())
corrected=loess(rdValues~gcValues, span=0.75,data.frame(rdValues=rdValues,gcValues=gcValues))

write(corrected,file=rdFile)

