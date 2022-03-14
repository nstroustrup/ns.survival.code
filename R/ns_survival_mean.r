ns_parse_survival_output = function(x){
	x$call = "";
	xx <- x;
	width = getOption("width");
	options(width = 10000);
	if (!is.null(x$strata))
		names(x$strata) = gsub(" ","&$%",names(x$strata))
	
	ox <- capture.output(print(x,print.rmean=TRUE,rmean="individual"));
	options(width = width);
	start_row = 3;
	colnames = strsplit(ox[start_row], split=' ')[[1]]
	missing_data = 0;
	if(colnames[9] == "missingness"){	#missing data error message is sometimes written first
		missing_data = as.numeric(colnames[4])
		start_row = 4;
		colnames = strsplit(ox[start_row], split=' ')[[1]]
		#browser()
	}
	colnames = colnames[colnames!=""]
	colnames = c("group",colnames);
	
	tmp <- strsplit(unlist(ox[(start_row+1):(length(ox)-1)[1]]), split=' ')
	#only one group
	if (length(ox) == 5)
		tmp[[1]][1]="group";
	for (i in 1:length(tmp))
		tmp[[i]] = tmp[[i]][tmp[[i]] != ""]
	tmp = matrix(unlist(tmp),length(tmp),length(tmp[[1]]),byrow=T)
	tmp[,1] = gsub("\\&\\$\\%"," ",tmp[,1])
	colnames(tmp) = colnames;
	tmp <- as.data.frame(tmp,stringsAsFactors = F);
	tmp$group = as.factor(tmp$group);
	tmp[["*rmean"]] = as.double(tmp[["*rmean"]]);
	tmp[["*se(rmean)"]] = as.double(tmp[["*se(rmean)"]]);
	tmp[["median"]] = as.double(tmp[["median"]]);
	tmp[["events"]] = as.double(tmp[["events"]]);
	tmp[["missing_data"]] = missing_data;
	v = strsplit(as.character(tmp[,1]),split="=");
	v = matrix(unlist(v),length(v),length(v[[1]]),byrow=T)
	tmp$group_short = v[,dim(v)[2]];
	#browser()
	return(tmp);
}
ns_survival_median <- function(x){
       return(ns_parse_survival_output(x)[["median"]]);
}
ns_survival_mean <- function(x){
       return(ns_parse_survival_output(x)[["*rmean"]]);
}

ns_survival_stderr <- function(x) {
	return(ns_parse_survival_output(x)[["*se(rmean)"]])
}

ns_survival_statistcs <- function(x) {
	return(ns_parse_survival_output(x))
}

ns_survival_median <- function(x){
       return(ns_parse_survival_output(x)[["median"]]);
}
