### NIH project - Placental miRNA
### miRNA-seq data processing and extraction 
### by F. White

#   miRNA-sequencing was conducted on plasma samples at first trimester. 
#   As a lookup with the dPCR analysis, the 4 miRNA of interest need to 
#   be extracted from the dataset.

# Usage: Rscript process_seq_data.R <outdir>
# Description:
# - library size normalization
# - CPM 
# - Extraction of 4 miRNA of interest

library(edgeR)
library(reshape2)

#======================================
### I/O
#======================================
args = commandArgs(trailingOnly=TRUE)
outdir = args[1]
if (is.na(outdir)) outdir = "data"
if (!dir.exists(outdir)) dir.create(outdir)

count_df = read.csv("/project/rrg-jacquesp-ab/Gen3G_share/data/freeze/plasma/mirnaseq_v1/exceRpt_miRNA_ReadCounts.mirbase21.clean.gct", sep="\t", check.names=F, row.names=1)

#======================================
### Global variables
#======================================
threshold_z_score = NA

# import functions   
if(!exists("format_dPCR_data", mode="function")) source("functions.R")

#======================================
### Data processing 
#======================================

# Normalisation TMM and CPM
dataObject = edgeR::DGEList(counts=count_df)
dataObject = edgeR::calcNormFactors(dataObject, method="TMM")
data = edgeR::cpm(dataObject, log=TRUE, prior.count=0.25) ### default

# extract miR of interest
data = data[grep(paste0(mir_list, collapse="|"), rownames(data)), ]
data = data.frame(t(data))

# clean column names
colnames(data) = gsub("hsa.", "", colnames(data))
colnames(data) = gsub("miR.517a.3p.miR.517b.3p", "miR.517a.3p", colnames(data))

# save data
data_wide = cbind(ID=rownames(data), data)
write.table(data_wide, file=paste0(outdir,"/seq.clean.wide.tsv"), sep="\t",row.names=F, quote=F)

data_long = reshape2::melt(data_wide, id.vars=c("ID"), variable.name="miR", value.name="normalized")
write.table(data_long, file=paste0(outdir,"/seq.clean.long.tsv"), sep="\t", row.names=F, quote=F)

outliers = plot_raw_data(data_long, "seq","outliers", threshold_z_score)
if (!is.na(threshold_z_score)) {
	cat(" * Exclusion de", length(outliers), "outliers (|z|>",threshold_z_score,"):\n")
	cat("  ", paste(outliers, collapse=", "), "\n")	
	#data_long = subset(data_long, !ID %in% outliers)
} 
