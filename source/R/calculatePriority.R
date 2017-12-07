isColumnPresent <- function(columnName, listOfColumns) {
  return(columnName %in% listOfColumns)
}

isExacScorePostive <- function(df, i) {
  if (df[i, columnNameConstants$EXAC] == ".")
    return(FALSE)
  exacScores = unlist(strsplit(df[i, columnNameConstants$EXAC], ","))
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
       df[i, columnNameConstants$LNCRNA] != ".")){
    return(TRUE)
  }
  return(FALSE)
}

getDecipherScore <- function(df,i){
  entries = unlist(strsplit(df[i,columnNameConstants$DECIPHER],","))
  minScore=101
  for(entry in entries){
    score = unlist(strsplit(entry,"\\|"))[3]
    score = as.numeric(substr(score,1,nchar(score)-1))
    if(score<minScore)
      score = minScore
  }
  return(minScore)
}

source("ColumnNames.R")
headers = "c1	c2	c3	c4	c5	c6\tGENE	LOCATION	Gene_Structural_Elements	DGV_accessions	Enhancers	miRNA_target_sites	Segmental_Duplications	Interspersed_repeats	Tandem_repeats	ClinVar_pathogenicity	ClinVar_Phenotype	OMIM_PHENOTYPE	DECIPHER	EXAC"
headers = matrix(unlist(strsplit(headers, "\t")))
df = read.delim2(
  file = 'temp.bed',
  sep = '\t',
  skip = 1,
  header = FALSE
)
colnames(df) = headers
columnNameConstants = ColumnNames()
dfColumns = colnames(df)
priorites = list()
HIGH = "high"
MEDIUM = "medium"
LOW = "low"

for (i in nrow(df)) {
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
  else if (isColumnPresent(columnNameConstants$EXAC, dfColumns) &&
           isExacScorePositive(df, i)) {
    priorites[i] = HIGH
    
  }
  else if (isColumnPresent(columnNameConstants$DECIPHER, dfColumns) &&
           getDecipherScore(df, i) <= 25) {
    priorites[i] = HIGH
    
  }
  # CHECK MEDIUM PROIRITY CONDITIONS
  # Medium priority: CNVs overlap with protein-coding regions/enhancers/miRNA target sites/lncRNA alone but not clincal data mentioned above.
  else if (hasAnyFunctionalOverlaps(df, i)) {
    priorites[i] = MEDIUM
  }
  # CHECK LOW PROIRITY CONDITIONS
  #  Low priority: CNVs have overlap with DGV annotations alone.
  #CNVs with negative ExAC values e.g., GENE-1.33 and CNVs with DECIPHER percentage > 25%
  else if (isColumnPresent(columnNameConstants$DGV_ACCESSIONS, dfColumns) &&
           df[i, columnNameConstants$DGV_ACCESSIONS] != ".") {
    priorites[i] = LOW
    
  }
  else if (isColumnPresent(columnNameConstants$EXAC, dfColumns) &&
           !isExacScorePositive(df, i, dfColumns)) {
    priorites[i] = LOW
    
  }
  else if (isColumnPresent(columnNameConstants$DECIPHER, dfColumns) &&
           getDecipherScore(df, i) >25) {
    priorites[i] = LOW
    
  }
}
