
ns_get_percentile_from_sfit = function(sfit,percentile){
	res = percentile;
	if (length(which(percentile > 1)) != 0 || length(which(percentile < 0)) != 0)
		stop("Invalid percentile");
	for (i in 1:length(percentile)){
		if (percentile[[i]] == 1){
			res[i] = max(sfit$time)
		}
		else if (percentile[[i]] == 0){
			res[i] = min(sfit$time)
		}
		else res[i] = quantile.survfit(sfit,percentile[[i]],conf.int=F)[1];
	}
	return (res);
	
}