column_names = list(
  GENE = "GENE",
  LOCATION = "LOCATION",
  GENE_STRUCTURAL_ELEMENTS = "Gene_Structural_Elements",
  DGV_ACCESSIONS = "DGV_accessions",
  ENHANCERS = "Enhancers",
  MIRNA_TARGET_SITES = "miRNA_target_sites",
  SEGMENTAL_DUPLICATION = "Segmental_Duplications",
  INTERSPERSED_REPEATS = "Interspersed_repeats",
  TANDEM_REPEATS = "Tandem_repeats",
  CLINVAR_PATHOGENICITY = "ClinVar_pathogenicity",
  CLINVAR_PHENOTYPE = "ClinVar_Phenotype",
  OMIM_PHENOTYPE = "OMIM_PHENOTYPE",
  DECIPHER = "DECIPHER",
  EXAC = "EXAC",
  LNCRNA = "lncRNA"
)
class(column_names) = "ColumnNames"
ColumnNames <- function(){
  column_names
}