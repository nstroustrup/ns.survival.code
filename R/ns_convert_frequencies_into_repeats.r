#' @export
ns_convert_frequencies_into_repeats<- function(raw_data,frequency_column="Event.Frequency"){
#	print(  raw_data[,frequency_column])
	temp = raw_data[rep(row.names(raw_data), raw_data[,frequency_column]), 1:ncol(raw_data)];
	temp[,frequency_column][temp[,frequency_column]>1] = 1;
	return(temp)
}