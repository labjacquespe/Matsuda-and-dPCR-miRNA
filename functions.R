
#======================================
### Global variables
#======================================

mir_list = c("miR.517a.3p","miR.524.5p","miR.1283","let.7b.3p")
mir_colors = c("#f8766d","#00b0f6","#00bf7d","#e76bf3")
hemolysis_list = c("N","1+","2+","3+")

library(ggplot2)
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
		#if(mir == mir_list[1]){
		#	new_df = mir_df
		#} else{
		#	new_df = rbind(new_df, mir_df)
		#}
	}
	# remove failed quantification
	#new_df = subset(new_df, normalized!=0)
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

plot_raw_data = function(df, dataset, name){
	# long format
	histograms = lapply(mir_list, function(mir){
		ggplot(subset(df, miR==mir), aes(x=normalized)) + 
			geom_histogram(bins=30) + 
	 		xlab(paste0("Normalized levels ", mir)) + ylab("Samples")})

	scatters = lapply(mir_list, function(mir){
		subdf = subset(df, miR==mir)
		ggplot(subdf, aes(x=as.character(ID), y=normalized, label=as.character(ID))) + 
			geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper), position=pd) + 
			geom_point(position=pd) + 
	  		geom_text_repel(data=subset(subdf, scale(normalized)>3)) + 
			theme(legend.position="none", axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
			labs(x="Samples", y=paste0("Normalized levels ", mir))})
	
	outliers = lapply(mir_list, function(mir){
		subdf = subset(df, miR==mir)
		subset(subdf, scale(normalized)>3)$ID
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
prepare_for_association = function(df, phenotypes){

 	# combine datafame with phenotypes
	df = merge(df, phenotypes, by="ID", all.x=TRUE) 

 	# z-score miRNA levels (scale function scales the columns)
	df[, mir_list] = scale(df[, mir_list])

 	# remove samples without main outcome (matsuda)
	complete_samples = na.omit(df[, c("ID", mir_list, "matsuda")])$ID
	df = subset(df, ID %in% complete_samples)
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
	report = as.data.frame(matrix(NA, nrow=nrow(test_list), ncol=9))
	colnames(report) = c("miR","model","subset","nb","beta","se","pvalue","lowerCI","upperCI")
	
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
		
		report[i,] = c(mir, paste0(outcome,"~",model), subset_name, stats)
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
	complete_fit = summary(fit)

	beta = round(complete_fit$coefficients[2,1], 3)
	se = round(complete_fit$coefficients[2,2], 3)
	pvalue = format(signif(complete_fit$coefficients[2,4], 3), scientific = T)
	lower = round(confint(fit)[2, 1], 3)
	upper = round(confint(fit)[2, 2], 3)

	# Exp for dichotomic outcomes
	if (outcome %in% dichotomic_outcomes){
		return(c(length(fit$fitted.values), exp(beta), exp(se), pvalue, exp(lower), exp(upper)))
	} else {
		return(c(length(fit$fitted.values), beta, se, pvalue, lower, upper))
	}
}

# Extract stats from model with interaction
extract_reg_interaction = function(fit, outcome){
	complete_fit = summary(fit)
	
	results = data.frame(complete_fit$coefficients)
	idx_row = grep(":", rownames(results)) # interaction term	

	beta = round(complete_fit$coefficients[idx_row,1], 3)
	se = round(complete_fit$coefficients[idx_row,2], 3)
	pvalue = format(signif(complete_fit$coefficients[idx_row, 4], 3), scientific = T)
	lower = round(confint(fit)[idx_row, 1], 3)
	upper = round(confint(fit)[idx_row, 2], 3)
	
	# Exp for dichotomic outcomes
	if (outcome %in% dichotomic_outcomes){
		return(c(length(fit$fitted.values), exp(beta), exp(se), pvalue, exp(lower), exp(upper)))
	} else {
		return(c(length(fit$fitted.values), beta, se, pvalue, lower, upper))
	}
}