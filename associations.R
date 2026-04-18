### Matsuda index and Placental miRNA 
# 		measured in plasma at first (T1) 
#			and second (T2) trimester of pregnancy  
### by F. White

# Usage: Rscript associations.R <outdir> <indir>

library(data.table)
library(ggplot2)
library(ggpubr)
#======================================
### I/O
#======================================
args = commandArgs(trailingOnly=TRUE)
outdir = args[1]
indir = args[2]

if (is.na(indir)) indir = "data"
if (!dir.exists(outdir)) dir.create(outdir)

file_phenotypes = "data/phenotypes.tsv"
file_T1 = "data/seq.clean.wide.tsv"
file_hemolysis = paste0(indir,"/dpcr.hemolysis.tsv")
file_T2 = paste0(indir,"/dpcr.clean.wide.tsv")

#======================================
### Global variables
#======================================

# import functions   
if(!exists("test_associations", mode="function")) source("Functions.R")

mir_list = c("miR.517a.3p","miR.524.5p","miR.1283","let.7b.3p")

continuous_var_list = c("BMI_T1","GA_T1","matsuda","insulin_cordblood",
	"cpeptide_cordblood","birthweight","placental_weight")

factor_var_list = c("hemolysis","sex","GDM")

subset_list = c("complete","Female","Male","NGT","GDM")
subset_order = c("GDM","NGT","Male","Female","complete")

timepoint_labels = c(fullT1="T1 (complete)", T1="T1 (overlap)", T2="T2")

outcome_labels = c(GDM="GDM", matsuda="Matsuda index",
	insulin_cordblood = "Insulin (cord blood)", cpeptide_cordblood="C-peptide (cord blood)",
	birthweight="Birthweight", placental_weight="Placental weight")

#======================================
# Open dataset
#======================================

T2 = fread(file_T2)
included_cols = c("ID", mir_list) 
T2 = T2[, ..included_cols]
setnames(T2, c("ID", paste0(mir_list, "_T2")))
main_dataset_id = T2$ID

full_T1 = fread(file_T1)
full_T1 = full_T1[, ..included_cols]
setnames(full_T1, c("ID", paste0(mir_list, "_T1")))
fullT1_dataset_id = full_T1$ID
T1_dataset_id = intersect(main_dataset_id, fullT1_dataset_id)

phenotypes = fread(file_phenotypes)
hemo = fread(file_hemolysis)
phenotypes = merge(phenotypes, hemo, by="ID", all.x=T)
included_cols = c("ID", factor_var_list, continuous_var_list) 
phenotypes = phenotypes[, ..included_cols]

data = merge(T2, full_T1, by="ID", all=T)
data = merge(data, phenotypes, by="ID", all=T)


#======================================
# Prepare variables
#======================================

# transformation
list_log2 = c("BMI_T1","matsuda","insulin_cordblood")
list_sqrt = c("placental_weight","cpeptide_cordblood")
for (var in continuous_var_list){
	if (var %in% list_log2) data[, (var) := log2(as.numeric(get(var)))]
	if (var %in% list_sqrt) data[, (var) := sqrt(as.numeric(get(var)))]
}
for( var in factor_var_list) data[, (var) := factor(get(var))]
write.table(data, file=paste0(outdir, "/data.tsv"), sep="\t", row.names=F, quote=F)

#======================================
### Tests T2
#======================================
cat("\n=== Testing T2 ===\n")
T2_models = c("+hemolysis+BMI_T1","*sex+hemolysis+BMI_T1", "*GDM+hemolysis+BMI_T1")
T2_report = test_associations(data, T2_models, "T2", main_dataset_id)

#======================================
# Test T1
#======================================
cat("\n=== Testing full T1 ===\n")
T1_models = c("+BMI_T1+GA_T1","*sex+BMI_T1+GA_T1", "*GDM+BMI_T1+GA_T1")
full_T1_report = test_associations(data, T1_models, "fullT1", fullT1_dataset_id)

cat("\n=== Testing overlapping T1 ===\n")
T1_report = test_associations(data, T1_models, "T1", T1_dataset_id)

#======================================
### Forest plots
#======================================
cat("\n=== Producing forest plots ===\n")

all_results = rbindlist(list(T2_report, T1_report, full_T1_report))

cols_to_num = c("beta", "lowerCI", "upperCI", "se", "pvalue")
all_results[, (cols_to_num) := lapply(.SD, as.numeric), .SDcols = cols_to_num]

all_results[, miR := gsub("_T1$|_T2$", "", miR)]
all_results[, `:=`(miR = factor(miR, levels = mir_list),
  				timepoint = factor(timepoint, levels = names(timepoint_labels)),
  				subset = factor(subset, levels = subset_order))]
all_results = all_results[model %in% c("+BMI_T1+GA_T1", "+hemolysis+BMI_T1")]
all_results[, signif := fcase(
  pvalue < 0.05,  "*",
  default = "NS"
)]
all_results[, signif := factor(signif, levels = c("NS", "*"))]

pdf(paste0(outdir,"/forest_plots.pdf"), height=6, width=8.5)
for (outcome_name in names(outcome_labels)){
	df_sub = all_results[outcome == outcome_name]	
	if (nrow(df_sub) == 0) next
	plots = create_forest_plot(df_sub, outcome_name)

	combined_plot = ggarrange(plotlist=plots, nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
	print(annotate_figure(combined_plot, 
		top=text_grob(outcome_labels[outcome_name], face="bold", size=16)))
}
dev.off()

