### Matsuda index and Placental miRNA 
# 		measured in plasma at first (T1) 
#			and second (T2) trimester of pregnancy  

### by F. White

# Usage: Rscript test_associations_from_supp <outdir> 

# The script use the published file from supplementary
# 	to replicate the study results.

library(data.table)
library(ggplot2)
library(ggpubr)
#======================================
### I/O
#======================================
args = commandArgs(trailingOnly=TRUE)
outdir = args[1]

input_file = "data/NIH miRNA project - Raw data.tsv"

data = fread(input_file, skip = 1)
data = data[, c(1:5, 14:17, 18:27)]

#======================================
### Functions
#======================================
run_association_tests_V2 = function(dt, models, timepoint_name){
	# create test grid
	test_list = CJ(model=models, mirna=paste0(mir_list, "_", timepoint_name), 
		outcome=outcome_list, subset=subset_list)

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
	return(report[!is.na(pvalue)])
}

run_reg = function(df, outcome, model_str, mir) {
	formula = as.formula(paste(outcome, "~", mir, model_str))

	fam = if(outcome %in% dichotomic_outcomes) "binomial" else "gaussian"
	fit = glm(formula, data = df, family = fam)

	sum_fit = summary(fit)$coefficients

	idx = if(grepl("*", model_str, fixed=TRUE)) grep(":", rownames(sum_fit)) else 2 
	if(length(idx) == 0) return(rep(NA, 6))

	res = c(n = length(fit$fitted.values), beta = sum_fit[idx, 1], 
		  se = sum_fit[idx, 2], p = sum_fit[idx, 4])
	ci = confint.default(fit)[idx, ] # .default est beaucoup plus rapide pour data.table

	final_stats = c(res[1], res[2:3], res[4], ci)
	if(outcome %in% dichotomic_outcomes) {
	final_stats[c(2,3,5,6)] = exp(final_stats[c(2,3,5,6)])
	}

	return(c(final_stats[1], round(final_stats[2:3], 3), 
		   format(signif(final_stats[4], 3), scientific = T), 
		   round(final_stats[5:6], 3)))
}

create_forest_plot = function(df, outcome_name) {
	is_dichotomic = outcome_name %in% dichotomic_outcomes

	ref_line = if(is_dichotomic) 1 else 0
	ymin     = if(is_dichotomic) 0 else -1
	ymax     = if(is_dichotomic) 2 else 1
	ylab     = if(is_dichotomic) "OR" else "Beta"

	if (is_dichotomic) {
		if (outcome_name == "sex") {
		  shape_values = c("GDM" = 4, "NGT" = 5, "complete" = 15)
		} else {
		  shape_values = c("Male" = 2, "Female" = 1, "complete" = 15)
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


#======================================
### Global variables
#======================================

mir_list = c("miR.517a.3p","miR.524.5p","miR.1283","let.7b.3p")

outcome_list = c("GDM", "matsuda", "birthweight", "placental_weight", "insulin_delivery", "c_peptide_delivery")
dichotomic_outcomes = c("GDM", "sex")

subset_list = c("complete","Female","Male","NGT","GDM")
subset_order = c("GDM","NGT","Male","Female","complete")

timepoint_list = c("T1", "T2")

#======================================
### Preprocessing
#======================================
setnames(data, c("GEO", paste0(mir_list, "_T2"), paste0(mir_list, "_T1"),
				 "hemolysis", "GA_V1", "BMI_V1", "sex", "GDM",
				 "matsuda", "birthweight", "placental_weight",
				 "insulin_delivery", "c_peptide_delivery"))

list_log2 = c("BMI_V1","matsuda","insulin_delivery")
list_sqrt = c("placental_weight","c_peptide_delivery")

for (var in outcome_list[-1]){
	# transformation
	if (var %in% list_log2) data[, (var) := log2(as.numeric(get(var)))]
	if (var %in% list_sqrt) data[, (var) := sqrt(as.numeric(get(var)))]

	# z-score
	data[, (var) := as.numeric(scale(get(var)))]

	# remove outliers
	data[abs(get(var)) > 4, (var) := NA]
}

for (mir in mir_list){

	# T1 : z-score
	m_t1 = paste0(mir, "_T1")
	data[, (m_t1) := as.numeric(scale(as.numeric(get(m_t1))))]

	# T2 : z-score + remove outliers
	m_t2 = paste0(mir, "_T2")
	data[, (m_t2) := as.numeric(scale(sqrt(as.numeric(get(m_t2)))))]
	data[abs(get(m_t2)) > 4, (m_t2) := NA]
}
data = data[!is.na(matsuda)]
write.table(data, file=paste0(outdir, "/processed_input.tsv"), sep="\t", row.names=F, quote=F)

#======================================
# Test associations
#======================================

cat("\n=== Testing T1 ===\n")
T1_models = c("+BMI_V1+GA_V1","*sex+BMI_V1+GA_V1", "*GDM+BMI_V1+GA_V1")
T1_report = run_association_tests_V2(data, T1_models, "T1")
write.table(T1_report, file=paste0(outdir,"/", "T1", "_associations.tsv"), 
		sep="\t", row.names=FALSE, quote=FALSE)

cat("\n=== Testing T2 ===\n")
T2_models = c("+hemolysis+BMI_V1","*sex+hemolysis+BMI_V1", "*GDM+hemolysis+BMI_V1")
T2_report = run_association_tests_V2(data, T2_models, "T2")


#======================================
### Forest plots
#======================================
cat("\n=== Producing forest plots ===\n")

outcome_labels = c(sex="Sex", GDM="GDM", matsuda="Matsuda index",
	birthweight="Birthweight", placental_weight="Placental weight",
	insulin_delivery = "Insulin (cord blood)", c_peptide_delivery="C-peptide (cord blood)")

T1_report$timepoint = "T1"
T2_report$timepoint = "T2"
all_data = rbindlist(list(T1_report, T2_report))
cols_to_num = c("beta", "lowerCI", "upperCI", "se", "pvalue")
all_data[, (cols_to_num) := lapply(.SD, as.numeric), .SDcols = cols_to_num]

all_data[, miR := gsub("_T1$|_T2$", "", miR)]
all_data[, `:=`(miR = factor(miR, levels = mir_list),
  				timepoint = factor(timepoint, levels = timepoint_list),
  				subset = factor(subset, levels = subset_order))]

all_data = all_data[model %in% c("+BMI_V1+GA_V1", "+hemolysis+BMI_V1")]

all_data[, signif := fcase(
  pvalue < 0.05,  "*",
  default = "NS"
)]
all_data[, signif := factor(signif, levels = c("NS", "*"))]

pdf(paste0(outdir, "/forest_plots.pdf"), height = 6, width = 5)
for (outcome_name in outcome_list) {
	df_sub = all_data[outcome == outcome_name]
	if (nrow(df_sub) == 0) next

	plots = create_forest_plot(df_sub, outcome_name)

	combined_plot = ggarrange(plotlist=plots, nrow=1, ncol=2, common.legend=TRUE, legend="bottom")

	print(annotate_figure(combined_plot, 
		top=text_grob(outcome_labels[outcome_name], face="bold", size=16)))
}
dev.off()
