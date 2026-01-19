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
subset_order = c("GDM","NGT","Male","Female","complete")
timepoint_list = c("T1 (complete)", "T1 (overlap)", "T2")
#======================================
# Open dataset
#======================================
# run analysis only if output does not already exist
if(!file.exists(paste0(outdir,"/T2_associations.tsv"))) {

	full_T1 = read.csv(file_T1, sep="\t")
	T2 = read.csv(file_T2, sep="\t")

	# T1 samples overlapping T2 samples (using hemolysis list because table T2 is missing the outliers)
	hemo = read.csv(file_hemolysis, sep="\t")
	sample_T1 = intersect(full_T1$ID, hemo$ID)
	T1 = subset(full_T1, ID %in% sample_T1)

	cat(" * number of samples at T1 (full): ", nrow(full_T1),"\n")
	cat(" * number of samples at T1 (overlap): ", nrow(T1),"\n")
	cat(" * number of samples at T2: ", nrow(T2),"\n")

	#======================================
	# Prepare phenotypes
	#======================================
	phenotypes = read.csv(file_phenotypes, sep="\t")

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

	# compute z-score
	# merge phenotypes with miRNA
	full_T1 = prepare_for_association(full_T1, phenotypes, "full_T1")
	T1 = prepare_for_association(T1, phenotypes, "T1")
	T2 = prepare_for_association(T2, phenotypes, "T2")
	
	hemo$hemolysis = as.numeric(factor(hemo$hemolysis, levels = hemolysis_list))
	T2 = merge(hemo, T2, by="ID", all=TRUE)
	
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

} else {
	cat("\n=== Associations already tested ===\n")
	
	full_T1 = read.csv(paste0(outdir, "/full_T1_data.tsv"), sep="\t")
	T1 = read.csv(paste0(outdir, "/T1_data.tsv"), sep="\t")
	T2  = read.csv(paste0(outdir, "/T2_data.tsv"), sep="\t")

	full_T1_report = read.csv(paste0(outdir,"/full_T1_associations.tsv"), sep="\t")	
	T1_report = read.csv(paste0(outdir,"/T1_associations.tsv"), sep="\t")	
	T2_report = read.csv(paste0(outdir,"/T2_associations.tsv"), sep="\t")	
}


#======================================
### Forest plots
#======================================
cat("\n=== Producing forest plots ===\n")

outcome_labels = c(sex = "Sex",
	GDM = "GDM",
	matsuda = "Matsuda index",
	birthweight = "Birthweight",
	placental_weight = "Placental weight",
	insulin_delivery = "Insulin (delivery)",
	c_peptide_delivery = "C-peptide (delivery)")

full_T1_report$timepoint = "T1 (complete)"
T1_report$timepoint = "T1 (overlap)"
T2_report$timepoint = "T2"
all_data = rbind(full_T1_report, T1_report, T2_report)

all_data$miR = factor(all_data$miR, levels=mir_list)
all_data$timepoint = factor(all_data$timepoint, levels=timepoint_list)
all_data$subset = factor(all_data$subset, levels = subset_order)
all_data = subset(all_data, model %in% c("+BMI_V1+GA_V1", "+hemolysis+BMI_V1"))

all_data$signif = ifelse(all_data$pvalue < 0.05, "*", "")
all_data$signif = ifelse(all_data$pvalue < 0.01, "**", all_data$signif)
all_data$signif = ifelse(all_data$pvalue < 0.001, "***", all_data$signif)
all_data$signif = factor(all_data$signif, levels = c("", "*", "**", "***"))

pdf(paste0(outdir,"/forest_plots.pdf"), height=7, width=9)
for (outcome_name in outcomes){

	df = subset(all_data, outcome == outcome_name)	
	plots = create_forest_plot(df, outcome_name)
	combined_plot = ggarrange(plotlist=plots, nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
	annotated_plot = annotate_figure(combined_plot,
		top=text_grob(outcome_labels[outcome_name], 
			face="bold", size=16))
	print(annotated_plot)
}
dev.off()

#======================================
### Boxplots 
#======================================
cat("\n=== Producing boxplots ===\n")

full_T1$timepoint = "T1 (complete)"
T1$timepoint = "T1 (overlap)"
T2$timepoint = "T2"
T2$hemolysis = NULL

all_quantification = unique(rbind(full_T1, T1, T2))
all_quantification$sex_factor = ifelse(all_quantification$sex==0, "M","F")
all_quantification$GDM_factor = ifelse(all_quantification$GDM==0, "NGT","GDM")
all_quantification$group = paste0(all_quantification$GDM_factor,"_",all_quantification$sex_factor)
all_quantification$group = factor(all_quantification$group, levels=c("NGT_M", "GDM_M", "NGT_F", "GDM_F"))

pdf(paste0(outdir,"/boxplots.pdf"), height=4, width=9)
for (timepoint_name in timepoint_list){
	
	plots = lapply(mir_list, function(mir){
		create_boxplot(subset(all_quantification, timepoint==timepoint_name), mir)})
	
	combined_plot = ggarrange(plotlist=plots, nrow=1, ncol=4, common.legend=TRUE, legend="bottom")
	print(annotate_figure(combined_plot, top=text_grob(timepoint_name, face="bold", size=16)))
}
dev.off()


#======================================
### Correlations over time points
#======================================
wide_quant = merge(T1[, c("ID", mir_list)], T2[, c("ID", mir_list, "sex", "GDM")], by="ID")
colnames(wide_quant) = c("ID", paste0("T1_", mir_list), paste0("T2_", mir_list),"sex", "GDM")
wide_quant$ID = NULL

NGT_df = subset(wide_quant, GDM == 0)
GDM_df = subset(wide_quant, GDM == 1)
male_df = subset(wide_quant, sex == 0)
female_df = subset(wide_quant, sex == 1)

#correlations = list()
subset_list = list(complete=wide_quant, NGT=NGT_df, GDM=GDM_df, Male=male_df, Female=female_df)
correlations = lapply(subset_list, function(df){
	M = cor(as.matrix(df[, grep("T1|T2", colnames(df))]))
	diag(M[grep("T2",rownames(M)), grep("T1", colnames(M))])
	})

correlations = data.frame(correlations)
colnames(correlations) = names(subset_list)
rownames(correlations) = mir_list

library(corrplot)
#Greys = c('#FFFFFF', '#F0F0F0', '#D9D9D9', '#BDBDBD', '#969696', '#737373', '#525252', '#252525', '#000000')
pdf(paste0(outdir, "/corrplot.pdf"), width=3.5, height=3.5)
corrplot(abs(as.matrix(correlations)), tl.col="black", addCoef.col=T, method="color", col=COL1(Greys, 10), col.lim = c(0, 0.5), is.corr=F)
#corrplot(abs(as.matrix(correlations)), tl.col="black", addCoef.col=T, method="color", col=Greys, col.lim = c(0, 0.5), is.corr=F)
dev.off()
