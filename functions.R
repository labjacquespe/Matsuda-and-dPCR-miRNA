
#======================================
### Global variables
#======================================

mir_list = c("miR.517a.3p","miR.524.5p","miR.1283","let.7b.3p")
hemolysis_list = c("N","1+","2+","3+")

library(ggplot2)
library(ggpubr)
library(ggrepel)
pd = position_dodge(0.9)


#======================================
### Global functions
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
	
	data$hemolysis = factor(data$hemolysis, levels=hemolysis_list)


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
		CI_low = ifelse(dataset=="dpcr", subdf$CI_lower, NA)
		CI_up = ifelse(dataset=="dpcr", subdf$CI_upper, NA)

		ggplot(subdf, aes(x=as.character(ID), y=normalized, label=as.character(ID))) + 
			geom_errorbar(aes(ymin=CI_low, ymax=CI_up), position=pd) + 
			geom_point(position=pd) + 
	  		geom_text_repel(data=subset(subdf, scale(normalized)>threshold_z_score)) + 
			theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
			labs(x="Samples", y=paste0("Normalized levels ", mir))})
	
	outliers = lapply(mir_list, function(mir){
		subdf = subset(df, miR==mir)
		subset(subdf, scale(normalized)>threshold_z_score)$ID
		})

	pdf(paste0("plot/process_",dataset,".",name,".pdf"), width=5, height=10)
	print(ggarrange(plotlist=histograms, nrow=4, ncol=1))
	print(ggarrange(plotlist=scatters, nrow=4, ncol=1))
	dev.off()
	return(unique(unlist(outliers)))
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
	complete_samples = na.omit(df[, c("ID", mir_list, "matsuda")])$ID
	df = subset(df, ID %in% complete_samples)
	
	write.table(df, file=paste0(outdir,"/", timepoint_name, "_data.tsv"), 
		sep="\t", row.names=FALSE, quote=FALSE)

	return(df)
}

# Filter dataframe according to the provided subset
get_subset = function(df, subset_name, outcome, model){
	can_be_tested = TRUE
	# Exclusion
	if (outcome == "GDM" & grepl("\\*GDM", model)) can_be_tested = FALSE
	if (outcome == "sex" & grepl("\\*sex", model)) can_be_tested = FALSE

	if(subset_name == "Male"){
		df = subset(df, sex==0)
		if (grepl("\\*sex", model)) can_be_tested = FALSE
	} else if (subset_name == "Female"){
		df = subset(df, sex==1)
		if (grepl("\\*sex", model)) can_be_tested = FALSE
	} else if (subset_name == "NGT"){
		df = subset(df, GDM==0)
		if (grepl("\\*GDM", model)) can_be_tested = FALSE
	} else if (subset_name == "GDM"){
		df = subset(df, GDM==1)
		if (grepl("\\*GDM", model)) can_be_tested = FALSE
	}

	if (subset_name %in% c("Male","Female") & outcome == "sex") can_be_tested = FALSE
	if (subset_name %in% c("NGT","GDM") & outcome == "GDM") can_be_tested = FALSE
	
	return(list(df, can_be_tested))	
}

#======================================
### Regression
#======================================


run_association_tests = function(df, models, timepoint_name){
	# create test grid
	test_list = expand.grid(
		model = models, 
		mirna = mir_list, 
		outcome = outcomes, 
		subset = c("complete","Male","Female","NGT","GDM"),
		stringsAsFactors = FALSE
	)

	# Exclusion
	test_list = subset(test_list, !(outcome=="sex" & subset %in% c("Male","Female")))
	test_list = subset(test_list, !(grepl("\\*sex", model) & subset %in% c("Male","Female")))
	test_list = subset(test_list, !(outcome=="GDM" & subset %in% c("GDM","NGT")))
	test_list = subset(test_list, !(grepl("\\*GDM", model) & subset %in% c("GDM","NGT")))
	test_list = subset(test_list, !(grepl("\\*GDM", model) & outcome == "GDM"))
	test_list = subset(test_list, !(grepl("\\*sex", model) & outcome == "sex"))
	cat(" * Test to perform:", nrow(test_list), "\n")
	
	# Initialise rapport
	report = as.data.frame(matrix(NA, nrow=nrow(test_list), ncol=10))
	colnames(report) = c("miR","outcome","model","subset","nb","beta","se","pvalue","lowerCI","upperCI")
	
	# for loop tests
	for (i in 1:nrow(test_list)){
		mir = test_list[i, "mirna"]
		model = test_list[i, "model"]
		outcome = test_list[i, "outcome"]
		subset_name = test_list[i, "subset"]
		
		# filter by strata
		result = get_subset(df, subset_name, outcome, model)
		data = result[[1]]
		can_be_tested = result[[2]]
		
		if (!can_be_tested) next
		
		# run regression
		if (outcome %in% dichotomic_outcomes){
			stats = test_logistic_reg(data, outcome, model, mir)
		} else {
			stats = test_linear_reg(data, outcome, model, mir)
		}
		
		report[i,] = c(mir, outcome, model, subset_name, stats)
	}
	
	# clean and write
	report = na.omit(report)
	write.table(report, file=paste0(outdir,"/", timepoint_name, "_associations.tsv"), 
		sep="\t", row.names=FALSE, quote=FALSE)
	
	return(report)
}

# run logistic regression
test_logistic_reg = function(df, outcome, model, mir){
	fit = glm(as.formula(paste0('df[, outcome]', "~",'df[, mir] ', model)), data=df, family="binomial")

	if(grepl("\\*", model)){
		return(extract_reg_interaction(fit, outcome))
	} else {
		return(extract_reg(fit, outcome))
	}
}

# run linear regression
test_linear_reg = function(df, outcome, model, mir){
	fit = lm(as.formula(paste0('df[, outcome]', "~",'df[, mir] ', model)), data=df)

	if(grepl("\\*", model)){
		return(extract_reg_interaction(fit, outcome))
	} else {
		return(extract_reg(fit, outcome))
	}
}

# Extract stats from model without interaction
extract_reg = function(fit, outcome){
	nb = length(fit$fitted.values)
	complete_fit = summary(fit)

	beta = complete_fit$coefficients[2,1]
	se = complete_fit$coefficients[2,2]
	pvalue = format(signif(complete_fit$coefficients[2,4], 3), scientific = T)
	lower = confint(fit)[2, 1]
	upper = confint(fit)[2, 2]

	# Exp for dichotomic outcomes
	if (outcome %in% dichotomic_outcomes){
		return(c(nb, 
			round(exp(beta), 3), 
			round(exp(se), 3),
			pvalue,
			round(exp(lower), 3),
			round(exp(upper), 3)))
	} else {
		return(c(nb,
		    round(beta, 3), 
		    round(se, 3), 
		    pvalue,
		    round(lower, 3), 
		    round(upper, 3)))
	}
}

# Extract stats from model with interaction
extract_reg_interaction = function(fit, outcome){
	nb = length(fit$fitted.values)
	complete_fit = summary(fit)
	
	results = data.frame(complete_fit$coefficients)
	idx_row = grep(":", rownames(results)) # interaction term	

	beta = complete_fit$coefficients[idx_row, 1]
	se = complete_fit$coefficients[idx_row, 2]
	pvalue = format(signif(complete_fit$coefficients[idx_row, 4], 3), scientific = T)
	lower = confint(fit)[idx_row, 1]
	upper = confint(fit)[idx_row, 2]
	
	# Exp for dichotomic outcomes
	if (outcome %in% dichotomic_outcomes){
		return(c(nb, 
			round(exp(beta), 3), 
			round(exp(se), 3),
			pvalue,
			round(exp(lower), 3),
			round(exp(upper), 3)))
	} else {
		return(c(nb,
		    round(beta, 3), 
		    round(se, 3), 
		    pvalue,
		    round(lower, 3), 
		    round(upper, 3)))
	}
}



create_forest_plot = function(df, outcome_name, dataset){
	ref_line = ifelse(outcome_name %in% dichotomic_outcomes, 1, 0)
	ymin = ifelse(outcome_name %in% dichotomic_outcomes, 0, -1)
	ymax = ifelse(outcome_name %in% dichotomic_outcomes, 2, 1)
	ylab = ifelse(outcome_name %in% dichotomic_outcomes, "OR", "Beta")

	if (outcome_name %in% dichotomic_outcomes){
		shape_values = if(outcome_name == "sex") {
			c(4, 5, 15)
		} else {
			c(2, 1, 15)
		}
		names(shape_values) = if(outcome_name == "sex") {
			c("GDM", "NGT","complete")
		} else {
			c("Male", "Female","complete")
		}
	} else {
		shape_values = c(4,5,2,1,15)
		names(shape_values) = subset_order
	}


	plots = list()
	for (timepoint_name in timepoint_list){
		subdf = subset(df, timepoint == timepoint_name)
		plots[[timepoint_name]] = ggplot(data=subdf, aes(x=subset, y=as.numeric(beta))) + 
			geom_point(size=2, aes(shape=subset, color=signif)) +
			facet_wrap(~miR, strip.position="left", nrow=9) + 
			geom_errorbar(aes(ymin=as.numeric(lowerCI), ymax=as.numeric(upperCI), color=signif), width=0.5, cex=0.5) + 
			geom_hline(yintercept=ref_line, linetype=3) + 
			ylim(ymin, ymax) + 
			ylab(ylab) + 
			theme_minimal() +
			theme(legend.position="none", axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.text.x=element_text(face="bold"), axis.title=element_text(size=12,face="bold"), strip.text.y = element_text(hjust=0, vjust=1, angle=180, face="bold")) + 
			coord_flip() +
			scale_colour_manual(values=c("grey75","grey30","grey20","black")) + 
			scale_shape_manual(values=shape_values) + 
			ggtitle(timepoint_name)
	}
	return(plots)
}


create_boxplot = function(df, mir_name){
	p = ggplot(df, aes(x=group, y=.data[[mir_name]], group=group, color=group)) +
		 geom_boxplot(outlier.shape=NA) +
		 geom_jitter(aes(shape=sex_factor), width=0.2) + 
		 stat_summary(fun=mean, geom="point", shape=3, size=3, color="black") + 
		 xlab("") + ylab(gsub("\\.", "-", mir_name)) + 
		 theme(legend.position="none", axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
		 ylim(-3, 3)
	return(p)
}

