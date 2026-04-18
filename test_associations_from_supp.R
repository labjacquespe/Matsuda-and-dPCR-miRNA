### Matsuda index and Placental miRNA 
# 		measured in plasma at first (T1) 
#			and second (T2) trimester of pregnancy  
### by F. White

# Usage: Rscript test_associations_from_supp.R <outdir> 

# This script is using the published file from supplementary material
# 	to replicate the study results.

library(data.table)
library(ggplot2)
library(ggpubr)
#======================================
### I/O
#======================================
args = commandArgs(trailingOnly=TRUE)
outdir = args[1]
if (!dir.exists(outdir)) dir.create(outdir)

input_file = "data/supp_table_raw_data.tsv"

data = fread(input_file)[,-2]
colnames(data)[1] = "ID"

#======================================
### Global variables
#======================================
if(!exists("test_associations", mode="function")) source("Functions.R")

mir_list = c("miR.517a.3p","miR.524.5p","miR.1283","let.7b.3p")

continuous_var_list = c("BMI_T1","GA_T1","matsuda","insulin_cordblood",
	"cpeptide_cordblood","birthweight","placental_weight")
factor_var_list = c("hemolysis","sex","GDM")

# trimester suffix
T1_suffix = "seq"
T2_suffix = "normalized"
timepoint_labels = c(T1="T1", T2="T2")

subset_list = c("complete","Female","Male","NGT","GDM")
subset_order = c("GDM","NGT","Male","Female","complete")

outcome_labels = c(GDM="GDM", matsuda="Matsuda index",
	insulin_cordblood = "Insulin (cord blood)", cpeptide_cordblood="C-peptide (cord blood)",
	birthweight="Birthweight", placental_weight="Placental weight")

#======================================
### Functions
#======================================
test_associations2 = function(dt, models, timepoint_name){
	# create test grid
	test_list = CJ(model=models, mirna=paste0(mir_list, "_", timepoint_name), 
		outcome=names(outcome_labels), subset=subset_list)

	# Exclusion
	test_list = subset(test_list, !(grepl("\\*sex", model) & subset %in% c("Male","Female")))
	test_list = subset(test_list, !(outcome=="GDM" & subset %in% c("GDM","NGT")))
	test_list = subset(test_list, !(grepl("\\*GDM", model) & subset %in% c("GDM","NGT")))
	test_list = subset(test_list, !(grepl("\\*GDM", model) & outcome == "GDM"))
	cat(" * Test to perform:", nrow(test_list), "\n")
	
	# create empty report
	report = data.table(miR = as.character(test_list$mirna),outcome = as.character(test_list$outcome),
		model = as.character(test_list$model),subset = as.character(test_list$subset),
		nb = NA_real_, beta = NA_real_, se = NA_real_, 
		pvalue = NA_real_, lowerCI = NA_real_, upperCI = NA_real_)

	# exclude samples without sequencing data at T1
	if(timepoint_name == "T1") dt = dt[!is.na(get(paste0(mir_list[1],"_T1")))]
	cat(" * n=",nrow(dt), " at ",timepoint_name, "\n")
	# z-score all continuous variables
	for (var in c(continuous_var_list, paste0(mir_list, "_", timepoint_name))){
		dt[, (var) := as.numeric(scale(as.numeric(get(var))))]
	}

	for (i in 1:nrow(test_list)){
		cur_test = test_list[i]

		sub_df = switch(cur_test$subset,
		  "Male"   = dt[sex == 0],
		  "Female" = dt[sex == 1],
		  "NGT"    = dt[GDM == 0],
		  "GDM"    = dt[GDM == 1],
		  dt) # Default : pas de filtre (All)
		
		# run regression
		stats = run_reg(sub_df, cur_test$outcome, cur_test$model, cur_test$mirna)
		set(report, i, 5:10, as.list(stats))
	}
	report = report[!is.na(pvalue)]
	write.table(report, file=paste0(outdir,"/", timepoint_name, "_associations.tsv"), 
		sep="\t", row.names=FALSE, quote=FALSE)

	report$timepoint = timepoint_name
	return(report)
}

run_reg2 = function(df, outcome, model_str, mir) {
	formula = as.formula(paste(outcome, "~", mir, model_str))

	fam = if(outcome %in% factor_var_list) "binomial" else "gaussian"
	fit = glm(formula, data = df, family = fam)

	sum_fit = summary(fit)$coefficients

	idx = if(grepl("*", model_str, fixed=TRUE)) grep(":", rownames(sum_fit)) else 2 
	if(length(idx) == 0) return(rep(NA, 6))

	res = c(n = length(fit$fitted.values), beta = sum_fit[idx, 1], 
		  se = sum_fit[idx, 2], p = sum_fit[idx, 4])
	ci = confint.default(fit)[idx, ] # .default est beaucoup plus rapide pour data.table

	final_stats = c(res[1], res[2:3], res[4], ci)
	if(outcome %in% factor_var_list) {
		final_stats[c(2,3,5,6)] = exp(final_stats[c(2,3,5,6)])
	}

	return(c(final_stats[1], round(final_stats[2:3], 3), 
		   format(signif(final_stats[4], 3), scientific = T), 
		   round(final_stats[5:6], 3)))
}

create_forest_plot2 = function(df, outcome_name) {
	is_dichotomic = outcome_name %in% factor_var_list

	ref_line = if(is_dichotomic) 1 else 0
	ymin     = if(is_dichotomic) 0 else -1
	ymax     = if(is_dichotomic) 2 else 1
	ylab     = if(is_dichotomic) "OR" else "Beta"

	if (is_dichotomic) {
		if (outcome_name == "GDM") {
		  shape_values = c("Male" = 2, "Female" = 1, "complete" = 15)
		} else {
		  shape_values = c("GDM" = 4, "NGT" = 5, "complete" = 15)
		}
	} else {
		shape_values = setNames(c(4, 5, 2, 1, 15), subset_order)
	}
	plots = list()

	for (timepoint_name in timepoint_list) {

		subdf = df[timepoint == timepoint_name]
		plots[[timepoint_name]] = ggplot(data=subdf, aes(x=subset, y=as.numeric(beta))) + 
			geom_point(size=2, aes(shape=subset, color=signif)) +
			facet_wrap(~miR, strip.position="left", nrow=9) + 
			geom_errorbar(aes(ymin=as.numeric(lowerCI), ymax=as.numeric(upperCI), color=signif), width=0.5, cex=0.5) + 
			geom_hline(yintercept=ref_line, linetype=2) + 
			ylim(ymin, ymax) + 
			ylab(ylab) + 
			xlab("")+
			theme_light() +
			theme(legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.ticks.x=element_blank(), axis.text.x=element_text(face="bold"), axis.title=element_text(size=12,face="bold"), strip.text.y = element_text(hjust=0, vjust=1, angle=180, face="bold"),  panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank()) + 
			coord_flip() +
			geom_vline(xintercept=c(2.5, 4.5), linetype=3) +
			scale_colour_manual(values=c("grey75","black")) + 
			scale_shape_manual(values=shape_values) + 
			ggtitle(timepoint_name)	
	}
	return(plots)
}

create_boxplot2 = function(df, mir_name, timepoint_name){
	mir_name = paste0(mir_name, "_", timepoint_name)
	colour_values = c("#0088aa", "#3caaff", "#ff5555", "#ff8080")
	names(colour_values) = subgroup_list
	p = ggplot(df, aes(x=group, y=.data[[mir_name]], group=group, color=group)) +
		 geom_boxplot(outlier.shape=NA) +
		 geom_jitter(aes(shape=sex_factor), width=0.2) + 
		 stat_summary(fun=mean, geom="point", shape=3, size=3, color="black") + 
		 xlab("") + ylab(gsub("\\.", "-", mir_name)) + 
		 theme_light() +
		 theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
		 ylim(-3, 3) +
		 scale_colour_manual(values=colour_values)
	return(p)
}

#======================================
### Preprocessing
#======================================
included_cols = c(
  "ID",
  paste0(T2_suffix, "_", mir_list),
  paste0(T1_suffix, "_", mir_list),
  factor_var_list,
  continuous_var_list)

data = data[, ..included_cols]

setnames(data, c(
	"ID", 
	paste0(mir_list, "_T2"),
	paste0(mir_list, "_T1"),
  factor_var_list,
  continuous_var_list))


# transformation
list_log2 = c("BMI_T1","matsuda","insulin_cordblood")
list_sqrt = c("placental_weight","cpeptide_cordblood", paste0(mir_list, "_T2"))
for (var in continuous_var_list){
	if (var %in% list_log2) data[, (var) := log2(as.numeric(get(var)))]
	if (var %in% list_sqrt) data[, (var) := sqrt(as.numeric(get(var)))]
}

for( var in factor_var_list) data[, (var) := factor(get(var))]
write.table(data, file=paste0(outdir, "/data.tsv"), sep="\t", row.names=F, quote=F)

#======================================
# Test associations
#======================================
cat("\n=== Testing T2 ===\n")
T2_models = c("+hemolysis+BMI_T1","*sex+hemolysis+BMI_T1", "*GDM+hemolysis+BMI_T1")
T2_report = test_associations(data, T2_models, "T2", data$ID)
write.table(T2_report, file=paste0(outdir,"/", "T2", "_associations.tsv"), 
		sep="\t", row.names=FALSE, quote=FALSE)


cat("\n=== Testing T1 ===\n")
T1_models = c("+BMI_T1+GA_T1","*sex+BMI_T1+GA_T1", "*GDM+BMI_T1+GA_T1")
overlap_id = data[!is.na(get(paste0(mir_list[1], "_T1")))]$ID
T1_report = test_associations(data, T1_models, "T1", overlap_id)

#======================================
### Forest plots
#======================================
cat("\n=== Producing forest plots ===\n")

all_results = rbindlist(list(T1_report, T2_report))
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

pdf(paste0(outdir, "/forest_plots.pdf"), height = 6, width = 5)
for (outcome_name in names(outcome_labels)) {
	df_sub = all_results[outcome == outcome_name]
	if (nrow(df_sub) == 0) next
	plots = create_forest_plot(df_sub, outcome_name)

	combined_plot = ggarrange(plotlist=plots, nrow=1, ncol=2, common.legend=TRUE, legend="bottom")

	print(annotate_figure(combined_plot, 
		top=text_grob(outcome_labels[outcome_name], face="bold", size=16)))
}
dev.off()



