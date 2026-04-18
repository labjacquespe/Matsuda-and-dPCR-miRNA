
#======================================
### Global variables
#======================================

library(ggplot2)
library(ggpubr)
library(ggrepel)
pd = position_dodge(0.9)

mir_list = c("miR.517a.3p","miR.524.5p","miR.1283","let.7b.3p")
hemolysis_list = c("N","1+","2+","3+","4+","5+")

#======================================
### Utility functions
#======================================

# Convert percent to intervals (divided by 2)
format_pct = function(list){
	list = as.numeric(sub(x=list, pattern="%", replacement=""))/100
	return(list/2)
}

# stack data into long format
stack_data = function(df){
	new_df = NULL
	for (mir in mir_list){
		mir_df = df[, c("ID", 
			paste0("normalized_",mir),
	 		paste0("CI_lower_",mir),
	 		paste0("CI_upper_",mir))]
		mir_df$miR = mir
		colnames(mir_df) = c("ID", "normalized", "CI_lower", "CI_upper", "miR")
		mir_df = mir_df[, c("ID", "miR","normalized", "CI_lower", "CI_upper")]
		new_df = rbind(new_df, mir_df)
	}
	# remove failed quantification (some readings show 0)
	new_df$normalized = ifelse(new_df$normalized==0, NA, new_df$normalized)

	new_df$uid = paste0(new_df$ID, "_", new_df$miR) 
	new_df$range = new_df$CI_upper - new_df$CI_lower 
	return(new_df)
}


# Keep value with lower absolute CI
keep_lower_range = function(df, value){
	result = df[, paste0(value,".x")]  # Default : batch 1
	# Replace by batch 2 if lower CI or data missing
	idx_batch2 = !is.na(df$range.y) & (df$range.y < df$range.x | is.na(df$range.x))
	result[idx_batch2] = df[idx_batch2, paste0(value,".y")]
	return(result)
}

#======================================
### dPCR processing functions
#======================================

format_dPCR_data <- function(input_file) {

	# open file and keep only data at fasting glucose (T0) 
	data = read.csv(input_file, sep="\t")
	data = subset(data, Time == "T0")
	data = data[, c( "ID", "hemolysis", mir_list, "cel.miR.39.3p",
			 paste0("normalized_", mir_list), 
			 paste0("CI_", mir_list), "CI_cel.miR.39.3p")]

	# convert percent CI to absolute CI 
	for (mir in c(mir_list, "cel.miR.39.3p")){
		data[, paste0("CI_", mir)] = format_pct(data[, paste0("CI_", mir)])
	}

	# compute CI for normalized values
	CI_lower_control = data[,"cel.miR.39.3p"] - data[,"CI_cel.miR.39.3p"]
	CI_upper_control = data[,"cel.miR.39.3p"] + data[,"CI_cel.miR.39.3p"]
	for (mir in c(mir_list)){
		CI_abs = data[, mir] * data[, paste0("CI_", mir)]
		min_mir = data[, mir] - CI_abs
		max_mir = data[, mir] + CI_abs
		data[, paste0("CI_lower_", mir)] = round(1000 * min_mir / CI_upper_control, 2)
		data[, paste0("CI_upper_", mir)] = round(1000 * max_mir / CI_lower_control, 2)
	}

	# remove useless columns
	data = data[, c("ID", "hemolysis",
		paste0("normalized_", mir_list), 
		paste0("CI_lower_", mir_list), 
		paste0("CI_upper_", mir_list)
		)]
	
	data$hemolysis = as.numeric(factor(data$hemolysis, levels=hemolysis_list)) - 1
	
	return(data)
}


# compare normalized levels and CI between batch for samples quantified twice
plot_rescued_data = function(df1, df2){

	sample_list = c(df1$ID, df2$ID)
	quantified_twice = sample_list[duplicated(sample_list)]

	# wide format
	df = rbind(cbind(df1, batch="1"), cbind(df2, batch="2"))
	df$uid = paste0(df$ID, "_", df$batch)
	df = subset(df, ID %in% quantified_twice)
	
	plots = lapply(mir_list, function(mir){
		ggplot(df, aes(x=uid, y=.data[[paste0("normalized_", mir)]], group=ID, color=batch)) + 
			geom_errorbar(aes(ymin=.data[[paste0("CI_lower_", mir)]], ymax=.data[[paste0("CI_upper_", mir)]]), position=pd) + 
			geom_point(position=pd) + 
			theme(legend.position="none",axis.text.x=element_text(angle=45, vjust=0.5)) + 
			labs(x="Rescued samples", y=paste0("Normalized levels ", mir))})
	
	pdf("plot/process_dpcr.rescued.pdf", width=5, height=10)
	print(ggarrange(plotlist=plots, nrow=4, ncol=1))
	dev.off()
}

plot_raw_data = function(df, dataset, name, threshold_z_score=NA){
	# long format
	histograms = lapply(mir_list, function(mir){
		ggplot(subset(df, miR==mir), aes(x=normalized)) + 
			geom_histogram(bins=30) + 
	 		xlab(paste0("Normalized levels ", mir)) + ylab("Samples")})


	scatters = lapply(mir_list, function(mir){
		subdf = subset(df, miR==mir)
		CI_low = ifelse(dataset=="dpcr", subdf$CI_lower, subdf$normalized)
		CI_up = ifelse(dataset=="dpcr", subdf$CI_upper, subdf$normalized)

		ggplot(subdf, aes(x=as.character(ID), y=as.numeric(normalized), label=as.character(ID))) + 
			geom_errorbar(aes(ymin=CI_low, ymax=CI_up), position=pd) + 
			geom_point(position=pd) + 
	  		geom_text_repel(data=subset(subdf, scale(normalized)>threshold_z_score)) + 
			theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
			labs(x="Samples", y=paste0("Normalized levels ", mir))})
	
	outliers = lapply(mir_list, function(mir){
		subdf = subset(df, miR==mir)
		subset(subdf, scale(normalized)>threshold_z_score)[, c("ID", "miR")]
		})

	pdf(paste0("plot/process_",dataset,".",name,".pdf"), width=5, height=10)
	print(ggarrange(plotlist=histograms, nrow=4, ncol=1))
	print(ggarrange(plotlist=scatters, nrow=4, ncol=1))
	dev.off()
	return(outliers)
}

#======================================
### Prepare for associations
#======================================
prepare_for_association = function(df, phenotypes, timepoint_name){

 	# combine datafame with phenotypes
	df = merge(df, phenotypes, by="ID", all.x=TRUE)

 	# z-score miRNA levels (scale function scales the columns)
	df[, mir_list] = scale(df[, mir_list])

 	# remove samples without main outcome (matsuda)
	complete_samples = na.omit(df[, c("ID", "matsuda")])$ID
	df = subset(df, ID %in% complete_samples)
	cat(" * samples included in ",timepoint_name, ": ", length(complete_samples),"\n")
	
	write.table(df, file=paste0(outdir,"/", timepoint_name, "_data.tsv"), 
		sep="\t", row.names=FALSE, quote=FALSE)

	return(df)
}

#======================================
### Regression
#======================================

test_associations = function(dt, models, timepoint_name, sample_list){

	if(timepoint_name=="fullT1"){
		time_suffix = "T1"
	} else {
		time_suffix = timepoint_name
	}

	### prepare tests
	# create test grid
	test_list = CJ(model=models, mirna=paste0(mir_list, "_", time_suffix), 
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

	### prepare dataset
	# included samples from tested dataset
	dt = dt[ID %in% sample_list]
	dt = dt[!is.na(matsuda)]
	cat(" * n=", nrow(dt), " at ",timepoint_name, "\n")
	
	# z-score all continuous variables
	for (var in c(continuous_var_list, paste0(mir_list, "_", time_suffix))){
		dt[, (var) := as.numeric(scale(as.numeric(get(var))))]
		dt[abs(get(var)) > 4, (var) := NA]
	}
	### run tests
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
	write.table(report, file=paste0(outdir,"/", timepoint_name, "_associations.tsv"), 
		sep="\t", row.names=FALSE, quote=FALSE)
	

	pdf(paste0(outdir,"/boxplots.",timepoint_name,".pdf"), height=3.5, width=7)
	plots = lapply(mir_list, function(mir){
				create_boxplot(dt, mir, timepoint_name)})
	combined_plot = ggarrange(plotlist=plots, nrow=1, ncol=4, common.legend=TRUE, legend="bottom")
	print(annotate_figure(combined_plot, 
			top=text_grob(timepoint_labels[timepoint_name], face="bold", size=16)))
	dev.off()

	report$timepoint = timepoint_name
	report = report[!is.na(pvalue)]
	return(report)
}

run_reg = function(df, outcome, model_str, mir) {
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

create_forest_plot = function(df, outcome_name) {
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

	for (timepoint_name in names(timepoint_labels)) {
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
			ggtitle(timepoint_labels[timepoint_name])	
	}
	return(plots)
}

create_boxplot = function(dt, mir_name, timepoint_name){
	
	subgroup_list = c("NGT_M", "GDM_M","NGT_F", "GDM_F")
	dt[, `:=`(sex_factor = ifelse(sex==0, "M", "F"),
			GDM_factor = ifelse(GDM==0,"NGT","GDM"))]
	dt[, group := factor(paste0(GDM_factor, "_", sex_factor),levels=subgroup_list)]

	if(timepoint_name=="fullT1"){
		time_suffix = "T1"
	} else {
		time_suffix = timepoint_name
	}

	mir_fullname = paste0(mir_name, "_", time_suffix)
	colour_values = c("#0088aa", "#3caaff", "#ff5555", "#ff8080")
	names(colour_values) = subgroup_list
	p = ggplot(dt, aes(x=group, y=.data[[mir_fullname]], group=group, color=group)) +
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

