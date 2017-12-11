isColumnPresent <- function(columnName, listOfColumns) {
  return(columnName %in% listOfColumns)
}

isExacScorePositive <- function(df, i) {
  exacScores = unlist(strsplit(as.character(df[i, columnNameConstants$EXAC]), ","))
  for (exacScore in exacScores) {
    if (grepl("\\+", exacScore))
      return(TRUE)
  }
  return(FALSE)
}

hasAnyFunctionalOverlaps <- function(df, i, dfColumns) {
  #protein-coding regions/enhancers/miRNA target sites/lncRNA
  if ((isColumnPresent(columnNameConstants$GENE, dfColumns) &&
       df[i, columnNameConstants$GENE] != ".") ||
      (isColumnPresent(columnNameConstants$MIRNA_TARGET_SITES, dfColumns) &&
       df[i, columnNameConstants$MIRNA_TARGET_SITES] != ".") ||
      (isColumnPresent(columnNameConstants$ENHANCERS, dfColumns) &&
       df[i, columnNameConstants$ENHANCERS] != ".") ||
      (isColumnPresent(columnNameConstants$LNCRNA, dfColumns) &&
       df[i, columnNameConstants$LNCRNA] != ".")) {
    return(TRUE)
  }
  return(FALSE)
}

getDecipherScore <- function(df, i) {
  entries = unlist(strsplit(as.character(df[i, columnNameConstants$DECIPHER]), ","))
  minScore = 101
  for (entry in entries) {
    score = unlist(strsplit(entry, "\\|"))[3]
    score = as.numeric(substr(score, 1, nchar(score) - 1))
    if (score < minScore)
      score = minScore
  }
  return(minScore)
}

args = commandArgs(TRUE)
fileName = args[1];
sourceDir = args[2];
source(paste(sourceDir,"ColumnNames.R",sep=""),chdir=FALSE)
df = read.csv(file=fileName,sep="\t",header=FALSE)
df = df[-length(df)]
colnames(df) = unlist(df[1,])
df = df[-1,]
columnNameConstants = ColumnNames()
dfColumns = colnames(df)
priorites = list()
HIGH = "high"
MEDIUM = "medium"
LOW = "low"

for (i in 1:nrow(df)) {
  priorites[i] = "."
  # CHECK HIGH PROIRITY CONDITIONS
  # High priority: Overlap with OMIM, ClinVar, positive ExAC score (GENE+1.45), DECIPHER 0-25%
  if (isColumnPresent(columnNameConstants$OMIM_PHENOTYPE, dfColumns) &&
      df[i, columnNameConstants$OMIM_PHENOTYPE] != ".") {
    priorites[i] = HIGH
    
  }
  else if (isColumnPresent(columnNameConstants$CLINVAR_PATHOGENICITY, dfColumns) &&
           df[i, columnNameConstants$CLINVAR_PATHOGENICITY] != ".") {
    priorites[i] = HIGH
    
  }
  else if (isColumnPresent(columnNameConstants$EXAC, dfColumns) && df[i, columnNameConstants$EXAC] != "." &&
           isExacScorePositive(df, i)) {
    priorites[i] = HIGH
    
  }
  else if (isColumnPresent(columnNameConstants$DECIPHER, dfColumns) && df[i,columnNameConstants$DECIPHER] != "." &&
           getDecipherScore(df, i) <= 25) {
    priorites[i] = HIGH
    
  }
  # CHECK MEDIUM PROIRITY CONDITIONS
  # Medium priority: CNVs overlap with protein-coding regions/enhancers/miRNA target sites/lncRNA alone but not clincal data mentioned above.
  else if (hasAnyFunctionalOverlaps(df, i, dfColumns)) {
    priorites[i] = MEDIUM
  }
  # CHECK LOW PROIRITY CONDITIONS
  #  Low priority: CNVs have overlap with DGV annotations alone.
  #CNVs with negative ExAC values e.g., GENE-1.33 and CNVs with DECIPHER percentage > 25%
  else if (isColumnPresent(columnNameConstants$DGV_ACCESSIONS, dfColumns) &&
           df[i, columnNameConstants$DGV_ACCESSIONS] != ".") {
    priorites[i] = LOW
    
  }
  else if (isColumnPresent(columnNameConstants$EXAC, dfColumns) && df[i, columnNameConstants$EXAC] != "." &&
           !isExacScorePositive(df, i, dfColumns)) {
    priorites[i] = LOW
    
  }
  else if (isColumnPresent(columnNameConstants$DECIPHER, dfColumns) && df[i,columnNameConstants$DECIPHER] != "." &&
           getDecipherScore(df, i) > 25) {
    priorites[i] = LOW
    
  }
}
df$Priority = unlist(priorites)
write.table(df,
            file = fileName,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
)
