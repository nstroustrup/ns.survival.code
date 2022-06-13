ns_survfit_truncate_trailing_censoring = function(formula){
	options(datadist=NULL)
	f = terms(formula);
	r = eval.parent(attr(f,"variables"))
	p = attr(f,"factors")
	group = r[[which(p==1)]]
	deaths = r[[which(p!=1)]]
	for (g in unique(group)){
		#largest death time
		max_death_time = max(deaths[group == g & deaths[,2]==1,1])
		deaths[group == g & deaths[,1]>max_death_time,1] = max_death_time;
	}
	return(survfit(deaths~group));
}
ns_get_rms_contrasts = function(model,reference_label,bj_group_levels){
	l = list();
	l[["bj_group"]] = as.character(reference_label)
	l2 = list();
	l2[["bj_group"]] = as.character(bj_group_levels)
	return(rms::contrast(model,l,l2,usebootcoef=TRUE))
}
ns_get_rms_contrasts_2 = function(model,reference_label,bj_group_levels){
	l = list();
	l[["bj_group_2"]] = as.character(reference_label)
	l2 = list();
	l2[["bj_group_2"]] = as.character(bj_group_levels)
	return(rms::contrast(model,l,l2,usebootcoef=TRUE))
}
ns_process_bj_result = function(ret_prop,return_coefficients,bj_group,reference_label,dbg=F){
	bj_residual = residuals(ret_prop,type="censored");
	
	if (return_coefficients){
		#bj produces incorrect coefficient labels
		#so we need to make them from scratch
		nn=levels(bj_group);
		nn = nn[nn!=as.character(reference_label)]
		#strip any 
		pos = regexpr("=",nn[2])[1]
		if (!is.na(pos))
		nn = c(as.character(reference_label),substring(nn,pos+1));

		#organize everything and return it.
		N = length(ret_prop$coefficients);
		
		contr = ns_get_rms_contrasts(ret_prop,nn,reference_label);
		coefficients = data.frame(bj_group=nn,
			      coefficient= contr$Contrast,
			      se = contr$SE,
			      lower=contr$Lower,
			      upper=contr$Upper,
			      p=contr$Pvalue);
		   if (dbg)browser()
		
		rownames(coefficients) = NULL;
		return(list(bj_residual = exp(bj_residual[,1]),
		    intercept =  coef(ret_prop)[1],
		    intercept_sd = sqrt(ret_prop$var[1,1]),
		    coefficients = coefficients,
		    model=ret_prop));		
	}else{
		return (exp(bj_residual[,1]));
	}
}
#' @export
ns_single_bj_group = function(formula,reference_group=NA,return_coefficients=T,dbg=F){
	options(datadist=NULL)
	f = terms(formula);
	r = eval.parent(attr(f,"variables"))
	p = attr(f,"factors")
	bj_group = r[[which(p==1)]]
	deaths = r[[which(p!=1)]]
	if (class(deaths) != "Surv")
		stop("Must provide a Surv object");
	bj_group = as.factor(as.character(bj_group))
	
	if (is.na(reference_group))
		reference_group = unique(bj_group)[1]
	reference_factor_level = which(as.character(levels(bj_group)) == as.character(reference_group))
	num_groups = length(unique(bj_group))
	#contrasts(bj_group) <- contr.treatment(num_groups,base=reference_factor_level);
	#now we run the regression
	ret_prop  <- bj(deaths~bj_group, x=TRUE, y=TRUE,time.inc=100000, control=list(iter.max=450,max.cycle=120,eps=10^-10))
	
	#bootcov(ret_prop)
#	if (dbg)browser()
	return(ns_process_bj_result(ret_prop,return_coefficients,bj_group,reference_group,dbg));
}
#' @export
ns_single_bj = function(deaths,
			reference_label=NA,
			bj_group_column="bj_group",
			time_column="Age.at.Death..d..Raw",
			censored_column="Censored",
			return_coefficients=T){
	options(datadist=NULL)		
	#first we set up the categorical variables for our regression
	if (!(bj_group_column %in% names(deaths)))
		stop("The parameter bj_group_column  must be specified")
	bj_group = as.factor(deaths[,bj_group_column])
	if (is.na(reference_label)){
		reference_label = unique(bj_group)[1]
	}
	reference_factor_level = which(as.character(levels(bj_group)) == as.character(reference_label))
	
	num_groups = length(unique(bj_group))
	#contrasts(bj_group) <- contr.treatment(num_groups,base=reference_factor_level);
	contrasts(bj_group) = NULL;
	#now we run the regression
	ret_prop  <- bj(Surv(deaths[,time_column],1-deaths[,censored_column])~
				bj_group, x=TRUE, y=TRUE,time.inc=.0001,
				 control=list(iter.max=500,max.cycle=30,eps=10^-5))
		 
	return(ns_process_bj_result(ret_prop,return_coefficients,bj_group,reference_label));
	    			      
}

ns_test_bj_effect_on_KS = function(deaths,bj_result,group_column,number_of_replicates=1000){
	groups = unique(deaths[,group_column]);
	if (length(groups) != 2)
		stop("Only two values should exist in the group column");
	group_1 = deaths[,group_column]==groups[[1]]
	Y = rep(NA,number_of_replicates)
	links_to_coeff = match(deaths$bj_group,bj_result$coefficients$bj_group)
	for (i in 1:number_of_replicates){
		tbj=  exp(bj_result$intercept+bj_result$coefficients$c + rnorm(dim(bj_result$coefficients)[1],0,1)*bj_result$coefficients$sd)
		td = deaths$Age.at.Death..d..Raw/tbj[links_to_coeff]
		s1 = Surv(td[group_1],1-deaths$Censored[group_1])
		s2 = Surv(td[!group_1],1-deaths$Censored[!group_1])
		dist = ns_ks_test(s1,s2,return_timeseries=F);
		Y[i] = dist$Y;	
	}
	return (Y);
}
#' @export
ns_multiple_bj = function(deaths,
			bj_group_column="bj_group",
			bj_group_2_column="Plate.Name",
			time_column="Age.at.Death..d..Raw",
			censored_column="Censored",
			reference_group=NA,
			reference_group_2=NA){
			
	options(datadist=NULL)

	#first we set up the categorical variables for our regression
	d = deaths;
	
	
	if (!(bj_group_column %in% names(d)))
		stop(paste("The bj group column,",bj_group_column,", must be specified in deaths data frame"))
	if (!(bj_group_2_column %in% names(d)))
		stop(paste("The bj group column 2, ",bj_group_2_column,", must be specified in deaths data frame"))
		
	bj_group = as.factor(as.character(d[,bj_group_column]))
	bj_group_2 = 	as.factor(as.character(d[,bj_group_2_column]))
	Censored = 1-d[,censored_column];
	
	t = d[,time_column]
	
	if (length(levels(bj_group_2)) < 2 && length(levels(bj_group)) < 2){
		#if there is only one device and one subject type, there is no need for a regression.
		#we just return the original death times.
		d$bj_residual = t
		d$device_corrected_death_time = t

		return(list(model = NA,
		    intercept = NA,
		    group_coefficients = NA,
		    device_coefficients = NA,
			deaths = d))
	}
	
	#set up the reference category for the regression
	if (!is.na(reference_group)){
		multi_aft_reference_group_label = as.character(reference_group);
	}else{
		multi_aft_reference_group_label	= as.character(levels(bj_group)[1]);
	}	
	multi_aft_reference_factor_level = which(as.character(levels(bj_group)) == multi_aft_reference_group_label)
	if (length(multi_aft_reference_factor_level) == 0)
			stop(paste("The specified reference group",multi_aft_reference_group_label,"is not present in the data"))
			
	if (length(levels(bj_group_2)) < 2){
		
		#if there is only one device, we can do a single regression
		#contrasts(bj_group) <- contr.treatment(length(levels(bj_group)),base=multi_aft_reference_factor_level);
		contrasts(bj_group) = NULL;
		d_surv = Surv(t,Censored);
		regression_data =data.frame(t=t,Censored = Censored,bj_group=bj_group)
		data_dist <<- datadist(regression_data);
		options(datadist="data_dist")
		scaling_device_bj_model <- bj(d_surv~bj_group, data=regression_data,time.inc=.0001, control=list(iter.max=60), x=TRUE, y=TRUE)
		group_coefficients = ns_get_rms_contrasts(scaling_device_bj_model,reference_group,levels(bj_group))
		group_coefficients_levels = levels(bj_group)
		
		r_prop <- resid(scaling_device_bj_model, type="censored")
		bj_residual = exp(r_prop[,1])
		#group_coefficients$N = 0;
		#for (i in 1:dim(group_coefficients)[1])
		#	group_coefficients$N[i] = sum(d$Event.Frequency[bj_group == group_coefficients$group[i] & Censored == 0])
		
		intercept = scaling_device_bj_model$coefficients[["Intercept"]]
		
		return(list(model = scaling_device_bj_model,
		    data_dist = data_dist,
		    regression_data = regression_data,
		    intercept = intercept,
		    group_coefficients = group_coefficients,
		    group_coefficients_levels=group_coefficients_levels,
		    device_coefficients = NA,
		    bj_residual = bj_residual))
	
	}
	#set up the reference device for the regression
	
	if (!is.na(reference_group_2)){
			multi_aft_reference_device_label = as.character(reference_group_2);
	}else multi_aft_reference_device_label = as.character(levels(bj_group_2)[1])
	multi_aft_reference_device_level = which(as.character(levels(bj_group_2)) == multi_aft_reference_device_label)
	if (length(multi_aft_reference_device_level) == 0){
		stop(paste("Requested reference device ",multi_aft_reference_device_label, "could not be found among levels:",paste(levels(bj_group_2))))	
	}
	#stop(multi_aft_reference_device_level)
	#contrasts(bj_group_2) <- contr.treatment(length(levels(bj_group_2)),base=multi_aft_reference_device_level);	
	contrasts(bj_group_2) = NULL;
	
	#if there is only one specified experimental covariate group,
	
	
	#we perform a single-covariate regression with each device as the covariate
	if (length(levels(bj_group)) < 2){
	
	
		regression_data = data.frame(t=t,Censored=Censored,bj_group=bj_group,bj_group_2=bj_group_2)
		data_dist <<- datadist(regression_data);
		options(datadist="data_dist")
	 	scaling_device_bj_model <- bj(Surv(t,Censored)~bj_group_2, time.inc=.0001, control=list(iter.max=60), x=TRUE, y=TRUE)
		group_2_coefficients = scaling_device_bj_model$coefficients
		group_2_n = levels(bj_group_2);
		group_n = levels(bj_group);
		regression_data = data.frame();
	}
	else{
		#if there are more than one covariate group, we 
		#perform a multivariate regression on each specified covariate and each devices.
		
		
		#contrasts(bj_group) <- contr.treatment(length(levels(bj_group)),base=multi_aft_reference_factor_level);
		contrasts(bj_group) = NULL;

		#set up ranges of variables for the RMS package to allow visualization of bj results
		regression_data = data.frame(t=t,Censored=Censored,bj_group=bj_group,bj_group_2=bj_group_2)
		data_dist <<- datadist(regression_data);
		options(datadist="data_dist")
		scaling_device_bj_model <- bj(Surv(t,Censored)~bj_group + bj_group_2, data=regression_data,time.inc=.0001, control=list(iter.max=60), x=TRUE, y=TRUE)
	
	}
	
	#if there was only one covariate
	#the coefficient for that covariate is 0,
	#and the intercept exists in the device coefficient array.
	if (length(levels(bj_group)) < 2){
		
		group_coefficients_levels = levels(bj_group)
		
		#we get a "no-effect" RMS contrast via subterfuge
		group_coefficients = ns_get_rms_contrasts_2(scaling_device_bj_model,multi_aft_reference_device_label,multi_aft_reference_device_label)
		
		intercept = scaling_device_bj_model$coefficients[["Intercept"]]
		group_2_coefficients = ns_get_rms_contrasts_2(scaling_device_bj_model,multi_aft_reference_device_label,levels(bj_group_2));
		
	}else{
		#if there were multiple covariates
		#the intercept exists in the group coefficient array.
		group_coefficients = ns_get_rms_contrasts(scaling_device_bj_model,multi_aft_reference_group_label,levels(bj_group))
		group_coefficients_levels = levels(bj_group)
		#browser()
		group_2_coefficients = ns_get_rms_contrasts_2(scaling_device_bj_model,multi_aft_reference_device_label,levels(bj_group_2))
		#print("AUTO")
		#browser()
		intercept = scaling_device_bj_model$coefficients[["Intercept"]]
	}
	group_coefficients$N = 0;
	for (i in 1:length(group_coefficients[["Contrast"]]))
		group_coefficients$N[i] = sum(d$Event.Frequency[bj_group == group_coefficients$group[i] & Censored == 0])
	
	rownames(group_2_coefficients) = NULL
	#calculate the residual of the single or multiple regression model
	r_prop <- resid(scaling_device_bj_model, type="censored")
	bj_residual = exp(r_prop[,1]);
	bj_corrected_death_time = rep(NA,length(bj_residual));
	bj_corrected_death_time_2 = rep(NA,length(bj_residual));
	#we can produce "device-normalized" death times
	#by centering the data around the intercept and 
	#adding back in the effect of each non-device category.
	
	all_bj_groups = unique(bj_group);
	for (i in 1:length(all_bj_groups)){
		indx = bj_group == all_bj_groups[[i]]
		cc = group_coefficients$Contrast[as.numeric(names(group_coefficients$Contrast)) == as.numeric(all_bj_groups[[i]]) ]
		#browser()
		bj_corrected_death_time[indx] = exp(r_prop[indx,1]+ cc + intercept);
	}
	
	all_bj_groups_2 = unique(bj_group_2);
	for (i in 1:length(all_bj_groups_2)){
		indx = bj_group_2 == all_bj_groups_2[[i]]
		cc = group_2_coefficients$Contrast[as.numeric(names(group_2_coefficients$Contrast)) == as.numeric(all_bj_groups_2[[i]])]
		bj_corrected_death_time_2[indx] = exp(r_prop[indx,1]+ cc + intercept);
	}

	
	#rm(data_dist)
	return(list(model = scaling_device_bj_model,
		    data_dist = data_dist,
		    regression_data = regression_data,
		    intercept = intercept,
		    group_coefficients = group_coefficients,
		    group_coefficients_levels = group_coefficients_levels,
		    group_2_coefficients = group_2_coefficients,
		    bj_residual = bj_residual,
		    bj_corrected_death_time=bj_corrected_death_time,
		    bj_corrected_death_time_2=bj_corrected_death_time_2))
}