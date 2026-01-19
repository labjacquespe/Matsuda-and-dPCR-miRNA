### NIH project - Placental miRNA
### Digital PCR data processing 
### by F. White

#   Digital PCR was conducted on candidate miRNAs and 
#   the miRNA control cel-miR-39-3p in plasma samples
#   collected during OGTT.

# Usage: Rscript process_dPCR_data.R 
# Description:
# - The script process_dPCR_data.R is a script to compute CI in dPCR data.
# - The script also combine the samples quantified two times and keep
#   the replicate with lower CI.
# - Output file are written in long (miRNA stacked) and wide format

library(reshape2)

#======================================
### I/O
#======================================

input_file_original = "data/Donnees_Cohorte_Gen3G.ORIGINAL_data.tsv"
input_file_rescued = "data/Donnees_Cohorte_Gen3G.REQUANTIFICATION_data.tsv"

#======================================
### Global variables
#======================================

threshold_z_score = 4

# import functions   
if(!exists("format_dPCR_data", mode="function")) source("functions.R")

#=============================================================================
### Format data
#=============================================================================

initial_data = format_dPCR_data(input_file_original)
rescued_data = format_dPCR_data(input_file_rescued)

# save hemolysis info
write.table(initial_data[, c("ID", "hemolysis")],
	file="data/dpcr.hemolysis.tsv", sep="\t", row.names=F, quote=F)

# convert to long format
stacked_initial = stack_data(initial_data)
stacked_rescued = stack_data(rescued_data)

stacked_data = merge(stacked_initial,
	stacked_rescued[, c("uid","normalized","CI_lower", "CI_upper","range")], by="uid", all=T)

#=============================================================================
### Keep best batch
#=============================================================================

stacked_data$normalized = keep_lower_range(stacked_data, "normalized")
stacked_data$CI_lower = keep_lower_range(stacked_data, "CI_lower")
stacked_data$CI_upper = keep_lower_range(stacked_data, "CI_upper")

data_long = stacked_data[, c("ID", "miR", "normalized", "CI_lower", "CI_upper")]

#=============================================================================
### Visualize data
#=============================================================================

plot_rescued_data(initial_data, rescued_data)

outliers = plot_raw_data(data_long, "dpcr","outliers", threshold_z_score)
cat(" * Exclusion of", length(outliers), "outliers (|z|>",threshold_z_score,"):\n")

cat("  ", paste(outliers, collapse=", "), "\n")
data_long = subset(data_long, !ID %in% outliers)
cat(" * including ",length(unique(data_long$ID)), "samples\n")

o = plot_raw_data(data_long, "dpcr", "clean")

#=============================================================================
### Write outputs
#=============================================================================

cat(" * applying square root transformation to get normal distribution\n")

data_long$normalized = sqrt(data_long$normalized)
data_long$CI_lower = sqrt(data_long$CI_lower)
data_long$CI_upper = sqrt(data_long$CI_upper)

data_wide = reshape2::dcast(data_long, ID~miR, value.var="normalized")
data_wide = data_wide[, c("ID", mir_list)]

write.table(data_long, file="data/dpcr.clean.long.tsv",
	sep="\t", row.names=F, quote=F)

write.table(data_wide, file="data/dpcr.clean.wide.tsv",
	sep="\t", row.names=F, quote=F)

cat(" * Données sauvegardées:\n")
cat("   - data/dpcr.clean.long.tsv\n")
cat("   - data/dpcr.clean.wide.tsv\n")

