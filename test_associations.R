### NIH project - Placental miRNA
### Association testing 
### by F. White

# Usage: Rscript test_associations.R <outdir>

#======================================
### I/O
#======================================

file_T1 = "data/seq.clean.wide.tsv"
file_T2 = "data/dpcr.clean.wide.tsv"
file_hemolysis = "data/dpcr.hemolysis.tsv"
file_phenotypes = "data/phenotypes.tsv"

args = commandArgs(trailingOnly=TRUE)
outdir = args[1]
if (!dir.exists(outdir)) dir.create(outdir)

#======================================
### Global variables
#======================================

# import functions   
if(!exists("format_dPCR_data", mode="function")) source("functions.R")

outcomes = c("sex", "GDM", "matsuda", "birthweight", "placental_weight",
	  "insulin_delivery", "c_peptide_delivery")
dichotomic_outcomes = c("GDM", "sex")

#======================================
# Open dataset
#======================================

full_T1 = read.csv(file_T1, sep="\t")
T2 = read.csv(file_T2, sep="\t")

# T1 with overlapping samples
sample_T1 = intersect(full_T1$ID, T2$ID)
T1 = subset(full_T1, ID %in% sample_T1)

cat(" * number of samples at T1 (full): ", nrow(full_T1),"\n")
cat(" * number of samples at T1 (overlap): ", nrow(T1),"\n")
cat(" * number of samples at T2: ", nrow(T2),"\n")

#======================================
# Prepare phenotypes
#======================================

# Phenotypes used for stratification: sex and GDM 
# Phenotypes used as adjustment variables: 
# - BMI_V1 (proxi for BMI before pregnancy)
# - GA_V1 (gestational age for T1) 
# - hemolysis (for T2)
# Phenotypes tested as outcomes:
# - sex
# - GDM
# - matsuda (measured at T2)
# - birthweight
# - placental_weight
# - insulin_delivery (measured in cordblood)
# - c_peptide_delivery (measured in cordblood)

phenotypes = read.csv(file_phenotypes, sep="\t")
hemo = read.csv(file_hemolysis, sep="\t")
hemo$hemolysis = as.numeric(factor(hemo$hemolysis, levels = hemolysis_list))
phenotypes = merge(hemo, phenotypes, by="ID", all=TRUE)

# normal distribution
phenotypes$BMI_V1 = log2(as.numeric(phenotypes$BMI_V1))
phenotypes$matsuda = scale(log2(as.numeric(phenotypes$matsuda)))
phenotypes$birthweight = scale(as.numeric(phenotypes$birthweight)^2)
phenotypes$placental_weight = scale(sqrt(as.numeric(phenotypes$placental_weight)))
phenotypes$insulin_delivery = scale(log2(as.numeric(phenotypes$insulin_delivery)))
phenotypes$c_peptide_delivery = scale(sqrt(as.numeric(phenotypes$c_peptide_delivery)))

# remove outliers (|z-score| > 4)
phenotypes$insulin_delivery = ifelse(abs(phenotypes$insulin_delivery) > 4, NA, phenotypes$insulin_delivery)
phenotypes$c_peptide_delivery = ifelse(abs(phenotypes$c_peptide_delivery) > 4, NA, phenotypes$c_peptide_delivery)

# merge phenotypes with miRNA
full_T1 = prepare_for_association(full_T1, phenotypes)
T1 = prepare_for_association(T1, phenotypes)
T2 = prepare_for_association(T2, phenotypes)

#======================================
# Test T1
#======================================

T1_models = c("+BMI_V1+GA_V1","*sex+BMI_V1+GA_V1", "*GDM+BMI_V1+GA_V1")

cat("\n=== Testing full T1 ===\n")
full_T1_report = run_association_tests(full_T1, T1_models, "full_T1")

cat("\n=== Testing overlapping T1 ===\n")
T1_report = run_association_tests(T1, T1_models, "T1")

#======================================
### Tests T2
#======================================

T2_models = c("+hemolysis+BMI_V1","*sex+hemolysis+BMI_V1", "*GDM+hemolysis+BMI_V1")

cat("\n=== Testing T2 ===\n")
T2_report = run_association_tests(T2, T2_models, "T2")

