ns_lighten<-function(x,val=.75){
	x = x+val;
	if (x>1)
	return(1);
	return (x);
}
ns_lighten_color = function(x,val=.75,opacity=1){
	x = col2rgb(x);
	err_bar_col = rgb(ns_lighten(x[1]/255),
			ns_lighten(x[2]/255),
			ns_lighten(x[3]/255),opacity);
}

#some code taken from muhaz
#' @export
ns_kp_haz = compiler::cmpfun(function (time, status, strata, q = 1, min_dt=1, method = "nelson",age_adjustment=data.frame(age=c(1),adjustment=c(1))) 
{
    if (missing(time)) 
        stop("Argument \"time\" is missing, with no default")
    if (any(is.na(time))) 
        stop("Time values can not be NA")
    if (any(is.nan(time))) 
        stop("Time values can not be Infinite")
    if (any(time < 0)) 
        stop("Time values must be >= 0")
    if (any(!is.numeric(time))) 
        stop("Time must be a numeric vector")
    if (missing(status)) 
        stop("Argument \"status\" is missing, with no default")
    if (any(is.na(status))) 
        stop("Status values can not be NA")
    if (any(!is.numeric(status))) 
        stop("Status must be a numeric vector")
    if (length(status) != length(time)) 
        stop("No. of observations in \"time\" and \"status\" must match")
    status <- as.integer(status)
    status[status != 0] <- 1
    if (all(status == 0)) 
        stop("No events occur in this data set")
    if (missing(strata)) 
        qstrata <- FALSE
    else {
        if (length(strata) != length(time)) 
            stop("\"Strata\" vector is the wrong length")
        qstrata <- TRUE
    }
    if (!is.numeric(q)) 
        stop("Agument \"q\" must be a numberic value")
    if (is.na(q)) 
        stop("q may not be NA")
    if (is.nan(q)) 
        stop("q may not be Infinite")
    q <- as.integer(q)
    if (q < 1) 
        stop("q must be positive")
    imethod <- pmatch(method, c("nelson", "product-limit"))
    if (is.na(imethod)) 
        stop("method must be one of \"nelson\" or \"product-limit\"")
    if (!qstrata) 
        strata <- rep(1, length(time))
    ind <- order(strata)
    strata <- strata[ind]
    time <- time[ind]
    status <- status[ind]
    ustrata <- unique(strata)
    dtime <- vector()
    haz <- vector()
    var <- vector()
    dstrata <- vector()
    
    for (j in 1:length(ustrata)) {
    	#set up variables containing data only for this strata
        cur.strata <- ustrata[j]
        cur.time <- time[strata == cur.strata]
        cur.status <- status[strata == cur.strata]
        ind <- order(cur.time)
        cur.time <- cur.time[ind]
        cur.n <- length(cur.time)
        cur.status <- cur.status[ind]
        cur.dtime <- unique(cur.time[cur.status != 0])
        cur.nd <- length(cur.dtime)
      
        cur.adj = 1/age_adjustment$adjustment[match(cur.time,age_adjustment$age,nomatch=1)]
       
        cur.weighted_n = sum(cur.adj);
        #first calculate the cumulative hazard function
        if (cur.nd > q) {
            H <- vector()
            VH <- vector()
            for (i in 1:cur.nd) {
      		temp.time = cur.time[cur.time <= cur.dtime[i]]
        	adj = cur.adj[cur.time <= cur.dtime[i]];
                temp.status <- (cur.status[cur.time <= cur.dtime[i]])*adj
                temp.n <- (1:length(temp.status));
                temp.n = cumsum(1/age_adjustment$adjustment[match(temp.time,age_adjustment$age,nomatch=1)])
  
                if (imethod == 1) 
                  H[i] <- sum(temp.status/(cur.weighted_n - temp.n + 1))
                else if (imethod == 2) 
                  H[i] <- sum(-log(1 - (temp.status/(cur.weighted_n - 
                    temp.n + 1))))
                if (imethod == 1) 
                  VH[i] <- sum(temp.status/((cur.weighted_n - temp.n + 
                    1) * (cur.weighted_n - temp.n + 1)))
                else VH[i] <- sum(temp.status/((cur.weighted_n - temp.n + 
                  1) * (cur.weighted_n - temp.n)))
            }
            #numeric differentiation
            for (i in seq(1, cur.nd - q)) {
            	#if (i%%100 == 0){
               #		print(paste(i,cur.nd - q))
                #	flush.console()
                #}
                dt = cur.dtime[(q+i):length(cur.dtime)]-cur.dtime[i];
                dQ = which(dt >= min_dt);
                if (length(dQ) == 0){
                	break;
                }else{
                  dQ = min(dQ)+q-1
                }
		dtime <- c(dtime, ((cur.dtime[dQ + i] + cur.dtime[i])/2))
		ttt <- (cur.dtime[dQ + i] - cur.dtime[i])
           
                haz <- c(haz, (H[dQ + i] - H[i])/ttt)
     
                var <- c(var, (VH[dQ + i] - VH[i])/(ttt^2))
                dstrata <- c(dstrata, cur.strata)
            }
          #  browser()
        }
    }
    time <- dtime
    strata <- dstrata
    if (qstrata) 
        return(list(time = time, haz = haz, var = var, strata = strata))
    else return(list(time = time, haz = haz, var = var))
})
#some code taken from muhaz
#' @export
ns_kp_cause_specific_haz = function (time, status, strata, q = 1, min_dt=1, method = "nelson",causes,age_adjustment=data.frame(age=c(1),adjustment=c(1))) 
{
    if (missing(time)) 
        stop("Argument \"time\" is missing, with no default")
    if (any(is.na(time))) 
        stop("Time values can not be NA")
    if (any(is.nan(time))) 
        stop("Time values can not be Infinite")
    if (any(time < 0)) 
        stop("Time values must be >= 0")
    if (any(!is.numeric(time))) 
        stop("Time must be a numeric vector")
    if (missing(status)) 
        stop("Argument \"status\" is missing, with no default")
    if (any(is.na(status))) 
        stop("Status values can not be NA")
    if (any(!is.numeric(status))) 
        stop("Status must be a numeric vector")
    if (length(status) != length(time)) 
        stop("No. of observations in \"time\" and \"status\" must match")
    status <- as.integer(status)
    status[status != 0] <- 1
    if (all(status == 0)) 
        stop("No events occur in this data set")
    if (missing(strata)) 
        qstrata <- FALSE
    else {
        if (length(strata) != length(time)) 
            stop("\"Strata\" vector is the wrong length")
        qstrata <- TRUE
    }
    if (!is.numeric(q)) 
        stop("Agument \"q\" must be a numberic value")
    if (is.na(q)) 
        stop("q may not be NA")
    if (is.nan(q)) 
        stop("q may not be Infinite")
    q <- as.integer(q)
    if (q < 1) 
        stop("q must be positive")
    imethod <- pmatch(method, c("nelson", "product-limit"))
    if (is.na(imethod)) 
        stop("method must be one of \"nelson\" or \"product-limit\"")
    if (!qstrata) 
        strata <- rep(1, length(time))
    ind <- order(strata)
    strata <- strata[ind]
    time <- time[ind]
    status <- status[ind]
    ustrata <- unique(strata)
    dtime <- vector()
    haz <- vector()
    var <- vector()
    dstrata <- vector()
    dcause <- vector()
    
    cause.levels = levels(causes)
    if (is.null(cause.levels))
    	stop("causes must be a factor");
    for (j in 1:length(ustrata)) {
    	#set up variables containing data only for this strata
        cur.strata <- ustrata[j]
        cur.time <- time[strata == cur.strata]
        cur.status <- status[strata == cur.strata]
        cur.causes = causes[strata == cur.strata]
        ind <- order(cur.time)
        cur.time <- cur.time[ind]
        cur.n <- length(cur.time)
        cur.status <- cur.status[ind]
        cur.causes = cur.causes[ind]
        cur.dtime <- unique(cur.time[cur.status != 0])
       
        cur.nd <- length(cur.dtime)
        
        
	cur.adj = 1/age_adjustment$adjustment[match(cur.time,age_adjustment$age,nomatch=1)]
	       
        cur.weighted_n = sum(cur.adj);
        
        #first calculate the cumulative hazard function
        if (cur.nd > q) {
            H <- matrix(0,nrow=cur.nd,ncol=length(cause.levels))
            VH <- matrix(0,nrow=cur.nd,ncol=length(cause.levels))
          
            for (cc in 1:length(cause.levels)){
        
		    cur.cause.status = cur.status;
		    cur.cause.status[cur.causes!=cause.levels[cc]]=0 #censor other causes
		    for (i in 1:cur.nd) {
		    
		    		temp.time = cur.time[cur.time <= cur.dtime[i]]
		        	adj = cur.adj[cur.time <= cur.dtime[i]];                

				temp.status <-   cur.cause.status[cur.time <= cur.dtime[i]]*adj
				temp.n <- 1:length(temp.status)
                		temp.n = cumsum(1/age_adjustment$adjustment[match(temp.time,age_adjustment$age,nomatch=1)])
                		
				if (imethod == 1) 
				  H[i,cc] <- sum(temp.status/(cur.weighted_n - temp.n + 1))
				else if (imethod == 2) 
				  H[i,cc] <- sum(-log(1 - (temp.status/(cur.weighted_n - 
				    temp.n + 1))))
				if (imethod == 1) 
				  VH[i,cc] <- sum(temp.status/((cur.weighted_n - temp.n + 
				    1) * (cur.n - temp.n + 1)))
				else VH[i,cc] <- sum(temp.status/((cur.weighted_n - temp.n + 
				  1) * (cur.n - temp.n)))
				  
     
			}
            }
            #numeric differentiation
            cc = 1:length(cause.levels)
    	    for (i in seq(1, cur.nd - q)) {
    	    
	              	#if (i%%100 == 0){
	                 #		print(paste(i,cur.nd - q))
	                  #	flush.console()
	                  #}
	                  dt = cur.dtime[(q+i):length(cur.dtime)]-cur.dtime[i];
	                  dQ = which(dt >= min_dt);
	                  if (length(dQ) == 0){
	                  	break;
	                  }else{
	                    dQ = min(dQ)+q-1
	                  }
	  		dtime <- c(dtime, rep(((cur.dtime[dQ + i] + cur.dtime[i])/2),length(cc)))
	  		ttt <- (cur.dtime[dQ + i] - cur.dtime[i])
	             
	                  haz <- c(haz, (H[dQ + i,cc] - H[i,cc])/ttt)
	                  var <- c(var, (VH[dQ + i,cc] - VH[i,cc])/(ttt^2))
	                  dstrata <- c(dstrata, cur.strata)
	                  dcause <- c(dcause,cc);
            }
        }
    }
    
        
    if (length(dtime)==0)
    	stop("Population size is less than specified q value");
    dcause = cause.levels[dcause]
    if (qstrata) 
        return(list(time = dtime, haz = haz, var = var, strata = dstrata, cause=dcause))
    else return(list(time = dtime, haz = haz, var = var,cause=dcause))
}



#' @export
ns_plot_hazard_by_groups <- function(deaths,group_column_name,
					level_colors,style_colors=NA,style_lty=NA,style_line_thickness=NA,
					style=c("km"),
					number_of_periods=12,
					minimum_period_duration_in_days=.02,
					parametric_fits=c(),
					draw_lines=c(T),draw_points=c(F),skip_N_less_than=0,
					logx=F,logy=T,xlim=c(.01,NA),add=FALSE,draw_legend=T,thick_center_half=F,
					curve_labels=NA,
					ylim=c(NA,NA),draw_grid_lines=T,
					plot_percentile_limits=c(.05,.95),point_size=1,muhaz_smoothness_factor=1,parameterization_legend=F,
					highlight_parametric_fit_region=F,legend_position="bottomright",
					quantile_label_style=c(NA),
					quantiles_to_label=c(),
					quantile_point_pch=20,
				    	#specificy how many death times should be used to calculate
					#each step of the cumulative hazard?
					#(fewer == smoother / less correct)
					khaz_span=50,
					#specificy minimum time step used to calculate the cumulative hazard
					#as an absolute time
					#(larger == smoother / less correct)
					khaz_min_dt = NA,
					#specificy minimum time step used to calculate the cumulative hazard
					#as a fraction of the total time during which deaths were observed
					#(fewer == smoother / less correct)
					khaz_number_of_steps = 20,
					#specify how many knots should be used to calculate the splines used
					#to calculate the error bars?
					#(fewer == smoother / less correct)
					km_spline_knots=60,
					#what's the lambda smoothness factor used to calculate the splines used
					#to calculate the error bars?
					#(fewer == smoother / less correct)
					km_spline_lambda=10,
					#what's the alpha smoothness factor used to calculate the splines used
					#to calculate the error bars?
					#(fewer == smoother / less correct),
					km_spline_a=.01,
				     plot_ci=T,
				     xlabel = "time (days)",ylabel="Hazard Rate (1/Days)",
				     level_order=NA,
				     text_size=par("cex"),
				     plot_spacing=c(2,4,.2,.2),
				     x_axis_tick=NA,
				     y_axis_tick=NA,
				     x_shift = 0){
	
	if(length(level_order) == 1 && is.na(level_order)){
		levels = unique(as.character(deaths[,group_column_name]));
		levels = levels[order(levels)]
	}else{
		if (length(level_order) != length(levels))
			stop(paste("Invalid level order length:",length(level_order),"supplied,",length(levels),"needed"))
		levels = as.character(level_order)
		if (length(which((levels %in% level_order) == FALSE))>0 ||
			length(which((level_order %in% levels) == FALSE))>0){
			#browser()
			#print(levels[which((levels %in% level_order) == FALSE)])
			#print(level_order)
			stop("Invalid level orders")
		}
	}
	
	if (length(khaz_span) != length(levels)){
		if (length(khaz_span) == 1){
			khaz_span = rep(khaz_span,length(levels));
		}else stop("khaz_span must have a length of 1 or be equal to the number of groups being plotted")
	}
	
	if (length(km_spline_knots) != length(levels)){
			if (length(km_spline_knots) == 1){
				km_spline_knots = rep(km_spline_knots,length(levels));
			}else stop("km_spline_knots must have a length of 1 or be equal to the number of groups being plotted")
	}
	
	if (length(levels) != length(level_colors))
		stop(paste("Provided ",length(level_colors), " colors for " , length(levels), " different values found in column ", group_column_name));;
	
	style_colors_specified = length(style_colors) > 1 || (length(style_colors) == 1 && !is.na(style_colors) )
	quantile_label_style_specified = length(quantile_label_style) > 1 || (length(quantile_label_style) == 1 && !is.na(quantile_label_style) )
	style_lty_specified = length(style_lty) > 1 || (length(style_lty) == 1 && !is.na(style_lty) )
	style_line_thickness_specified = length(style_line_thickness) > 1 || (length(style_line_thickness) == 1 && !is.na(style_line_thickness) )
	
	point_type = 21
	
	if (style_colors_specified && length(style_colors) != length(style))
		stop("Style colors must have the same length as styles")
	if (style_lty_specified && length(style_lty) != length(style))
		stop("Style lty must have the same length as styles")
	if (style_line_thickness_specified && length(style_line_thickness) != length(style))
		stop("Style line thickness must have the same length as styles")
	if (quantile_label_style_specified && length(quantile_label_style) != length(style))
		stop("Quantile label style must have the same length as styles")
	
	non_parametric_plot_types = c("km","period","mu");
	parametric_fit_style_i =  !(style %in% non_parametric_plot_types);
	parameterization = style[ parametric_fit_style_i]
	parametric_colors=c();
	parametric_lty=c();
	if (style_colors_specified)
		parametric_colors = style_colors[parametric_fit_style_i];
	if (style_lty_specified)
		parametric_lty = style_lty[parametric_fit_style_i];
	if (style_line_thickness_specified)
		parametric_line_thickness = style_line_thickness[parametric_fit_style_i]
	
	#stop(number_of_periods)
	if (length(number_of_periods) == 1){
		number_of_periods = rep(number_of_periods,length(levels))
	}
	else{
		if (length(number_of_periods) != length(levels)){
			stop("number_of_periods must either be of length one or correspond to the number of distinct curves being plotted");
		}
	
	}
	#default_colors = c("red","blue","green","orange","purple");
	par(las=1)
	if (length(parametric_fits) != 0){
		#identify and remove lifespan data for which we don't have fits
		
		if (!(group_column_name %in% colnames(parametric_fits))){
			stop(paste("Column",group_column_name, "could not be found in parametric_fits"));
		}
	}

	only_deaths = deaths$Age.at.Death..d..Raw[deaths$Censored==0];
	#plot(deaths$Age.at.Death..d..Raw[deaths$Censored == 0])
	xl = c(quantile(only_deaths,0.005),
	       quantile(only_deaths,0.995));
	yl= c(.001,50);
	if (length(ylim)==2 && !is.na(ylim[1]))
		yl[1] = ylim[1];
	if (length(ylim)==2 && !is.na(ylim[2]))
		yl[2] = ylim[2];
	
	if (length(ylim)==2 && !is.na(xlim[1]))
		xl[1] = xlim[1]
	if (length(ylim)==2 && !is.na(xlim[2]))
		xl[2] = xlim[2]
	
	#print(paste("Plotting between ", xl[1], " and ", xl[2]))
	mt_1 = 0;

	
	if (length(muhaz_smoothness_factor) == 1){
		muhaz_smoothness_factor = rep(muhaz_smoothness_factor,length(levels));
	}
	plotted = F
	if (add)
		plotted = T;
	print(paste("Plotting ",length(levels),"hazard functions"))
	hazard = list();
	#levels = levels[order(as.character(levels))]
	hazard[["groups"]] = levels;
	for (i in 1:length(levels)){
	
		#We can draw quantile labels according to each hazard estimation method
		#to do that, we need to first collect the full hazard curve
		#so that we can subsequently plot the correct data
		hazard_values_to_use_for_labels = list();
		
		
		level_color = level_colors[i];
		level_lty = 1;
		level_line_thickness = 1;
		
		print(paste("Plotting group",levels[[i]]))
		current_deaths = subset(deaths,deaths[,group_column_name] == levels[[i]]);
		#only_deaths_sub = current_deaths$Age.at.Death..d..Raw[current_deaths$Censored==0];
		#dl = c(min(only_deaths_sub),max(only_deaths_sub));
		#print(plot_percentile_limits)
		s = Surv(current_deaths$Age.at.Death..d..Raw,1-current_deaths$Censored);
		sfit = survfit(s~1);
		time_bounds = ns_get_percentile_from_sfit(sfit,plot_percentile_limits);
		print(time_bounds)
		if (is.na(time_bounds[2]))
			time_bounds[2] = max(current_deaths$Age.at.Death..d..Raw)
		print(time_bounds)
		#time_bounds = c(quantile(only_deaths_sub,plot_percentile_limits[1])[[1]],quantile(only_deaths_sub,plot_percentile_limits[2])[[1]]);
		#browser()
		#majority_time_bounds = c(quantile(only_deaths_sub,.1)[[1]],quantile(only_deaths_sub,.9)[[1]]);
		majority_time_bounds = ns_get_percentile_from_sfit(sfit,c(.1,.9));
		
		p_width = (time_bounds[2]-time_bounds[1])/number_of_periods[i]
		
		#print(paste("COUNT",dim(d),"MIN ", dl[1],"MAX",dl[2],"WIDTH:",p_width));
		#print(paste("UNIQUE",unique(d$Age.at.Death..d..Raw[d$Censored==0])));

		if (!plotted){
			plotted = T
		#	print(i)
			dcex = par("cex")
			par(cex=text_size)
			ll = "";
			if (logx)
				ll = "x";
			if (logy)
				ll=paste(ll,"y",sep="")
			#stop(ll)
			#print(paste("XL:",xl))

			par(mar=plot_spacing)  
			plot(x=c(1),y=c(1),type="n",log=ll, ylim=yl,xlim=xl,xlab="",ylab="",axes=F)
			
			if (length(grep("x",ll))==0){
					xlm = seq(xl[1],xl[2],(xl[2]-xl[1])/4)
			}else{
				#print(paste("xl",paste(xl)))
				
				xl1 = floor(4*log10(xl[1]))/4;
				xl2 = ceiling(4*log10(xl[2]))/4;
				
				if (xl1 > 1) xl1 = log10(floor(10^xl1))
				if (xl2 > 1) xl2 = log10(ceiling(10^xl2))
				dx = (10^xl2-10^xl1)/8
				xlm = seq(10^xl1,10^xl2,dx)
			}
			if(!is.na(x_axis_tick))
				xlm = x_axis_tick
		
			if (length(grep("y",ll))==0){
				ylm = seq(yl[1],yl[2],(yl[2]-yl[1])/4)
			} else{
				yl1 = floor(log10(yl[1]));
				yl2 = ceiling(log10(yl[2]));
				ylm = 10^(seq(yl1,yl2,1))
			}
			if(!is.na(y_axis_tick))
				ylm = y_axis_tick

			mgp.axis(1,xlm,axistitle=xlabel,mgp=c(2,.5,0))
			mgp.axis(2,ylm,axistitle=ylabel,mgp=c(6,1,0))
			par(cex=dcex)
			#mgp.axis(1,xlabel)
			#mgp.axis(2,ylabel)
			if (0&draw_grid_lines){
				for (n in -4:3){
					abline(h=10^n,col="#EEEEEE");
				}
				ab1 = ceiling(xl[1]*2)/2;
				ab2 = floor(xl[2]*2)/2;
				if (ab1 < ab2)
					for (n in seq(ceiling(xl[1]*2)/2,floor(xl[2]*2)/2,.5))
						abline(v=n,col="#EEEEEE");
			}
		}

		if (length(current_deaths$Age.at.Death..d..Raw) < skip_N_less_than){
			warning(paste("Skipping curve contaning only ", length(current_deaths$Age.at.Death..d..Raw), " animals as it smaller than specified minimum:", skip_N_less_than ))
			next;
		}
		#browser()
		if (p_width == 0 || dim(current_deaths)[1] == 0)
			next
			
			
		#if (!color_by_groups)
		#	level_colors[[i]] = "#000000"
		ccex = .5*point_size;
			
		if ("km" %in% style){
		
			style_i = which(style=="km");
			if (style_colors_specified && !is.na(style_colors[style_i]))
				level_color = style_colors[style_i]
			if (style_lty_specified && !is.na(style_lty[style_i]))
				level_lty = style_lty[style_i]
			if (style_line_thickness_specified && !is.na(style_line_thickness[style_i]))
				level_line_thickness = style_line_thickness[style_i]
			
			draw_l = draw_lines[which(style=="km")]
			draw_p <- draw_points[which(style=="km")]
			
			
			
			if (!is.na(khaz_number_of_steps) && !is.na(khaz_min_dt))
				stop("khaz_number_of_steps and khaz_min_dt cannot both be specified");
			if (!is.na(khaz_number_of_steps)){
				tr = c(min(current_deaths$Age.at.Death..d..Raw),
					max(current_deaths$Age.at.Death..d..Raw))
				min_dt = diff(tr)/khaz_number_of_steps;
			}else if (!is.na(khaz_min_dt)){
				min_dt = khaz_min_dt;
			} else stop("Either khaz_number_of_steps or khaz_min_dt must be specified");
			#browser();
			res <- ns_kp_haz(current_deaths$Age.at.Death..d..Raw+x_shift,min_dt=min_dt,
				status=1-current_deaths$Censored,q=khaz_span[[i]],method="product-limit")
		
			t= res[["time"]]-x_shift;
			h = res[["haz"]];
			v_h = sqrt(res[["var"]])
			h_u = h+sqrt(v_h)
			h_l = h-sqrt(v_h)
			sub = is.infinite(h) | is.na(h) | is.na(h_u) | is.na(h_l)# | h_l <= 0;
			t = t[!sub];
			h = h[!sub];
			h_l = h_l[!sub]
			h_u = h_u[!sub]
				
			
			time_sub_index = t >= time_bounds[1] & t < time_bounds[2];
			t = t[time_sub_index];
			h = h[time_sub_index];
			h_u = h_u[time_sub_index];
			h_l = h_l[time_sub_index];
			#browser()
			fit<-bshazard(Surv(current_deaths$Age.at.Death..d..Raw+x_shift, 1-current_deaths$Censored) ~ 1,
				degree=4,lambda=km_spline_lambda,alpha=km_spline_a,nk=km_spline_knots[[i]])
			fitt = summary(fit);
			fit_h = fitt$HazardEstimates
			poly_haz = data.frame(t=fitt$HazardEstimates[,1]-x_shift,
					   h=fitt$HazardEstimates[,2],
					   l = fitt$HazardEstimates[,3],
					   u = fitt$HazardEstimates[,4])
			poly_haz = subset(poly_haz, poly_haz$t >= time_bounds[1] & 
					           poly_haz$t < time_bounds[2]);   
			
			#plot(t,h,log='xy',type='n')
			#lines(t,h,col=level_color,lty=level_lty)
			#lines(poly_haz$t,poly_haz$h,col="red")
			#browser()
			hazard_values_to_use_for_labels[[style_i]] = data.frame(t=t,h=h);

			majority_time_bounds_index =  t >= majority_time_bounds[1] & t <= majority_time_bounds[2]

			t_offset = 0;
			if (!is.na(parametric_fits) && !is.null(parametric_fits$time_offset))
				t_offset = parametric_fits[parametric_fits[,group_column_name] == levels[[i]],]$time_offset;

			if (plot_ci){
				polygon(c(poly_haz$t+t_offset,rev(poly_haz$t+t_offset)),
					 c(poly_haz$u,rev(poly_haz$l)),
					 col=ns_lighten_color(level_color,1,.5),lty=0)
				#polygon(c(t,rev(t)),c(h_l,rev(h_u)),col=ns_lighten_color(level_color,1,.5),lty=0)
				#lines(t,c(h_l),col=ns_lighten_color(level_color,1,.8),lty=level_lty)
				#lines(t,h_u,col=ns_lighten_color(level_color,1,.8),lty=level_lty)
			
			}
			if (draw_p)
				points(t+t_offset,h,col=level_color,lty=level_lty,pch=point_type,bg=level_color,cex=ccex);
			if(thick_center_half){
				if (draw_p)
				points(t[majority_time_bounds_index]+t_offset,h[majority_time_bounds_index],col=level_color,bg=level_color,lty=level_lty,pch=point_type,cex=2*ccex,lwd=level_line_thickness);
				if (draw_l)
				lines(t[majority_time_bounds_index]+t_offset,h[majority_time_bounds_index],col=level_color,bg=level_color,lty=level_lty,pch=point_type,lwd=level_line_thickness);
			}
			else{
			if (draw_l)
				lines(t+t_offset,h,pch=point_type,bg=level_color,cex=ccex,col=level_color,lty=level_lty,lwd=level_line_thickness);
			}

		}

		hazard[["period"]] = list();

		#hazard_values = data.frame(t=x,h=values);
		if ("period" %in% style){
			style_i = which(style=="period");
			if (style_colors_specified && !is.na(style_colors[style_i]))
				level_color = style_colors[style_i]
			if (style_lty_specified && !is.na(style_lty[style_i]))
				level_lty = style_lty[style_i]
			if (style_line_thickness_specified && !is.na(style_line_thickness[style_i]))
				level_line_thickness = style_line_thickness[style_i]
			
			draw_l = draw_lines[which(style=="period")]
			draw_p <- draw_points[which(style=="period")]
			
			if (0){
			#print(draw_lines)
			#stop()
			capture.output(
				p_a1 <- pehaz(current_deaths$Age.at.Death..d..Raw+x_shift, 1-current_deaths$Censored,
				width=p_width,
				min.time=min(current_deaths$Age.at.Death..d..Raw+x_shift),
				max.time=max(current_deaths$Age.at.Death..d..Raw+x_shift)
				)
			);
			
			#if (length(parameterization) != 0) ccex=.8*point_size;
			x = p_a1$Cuts[1:(length(p_a1$Cuts)-1)]-x_shift;
			#x= rollmean(p_a1$Cuts,2)[indx]-x_shift;
			vals = p_a1$Hazard; #p_a1$Events/p_a1$At.Risk/p_a1$Width;
			}


			if (1){
				pop_size = dim(current_deaths)[1];
				#we do *not* include censored times in our percentile calculations
				#as we want to have an equal number of deaths in each group
				s = survfit(Surv(current_deaths$Age.at.Death..d..Raw+x_shift)~1)
				percentile_width = 1/number_of_periods[i]
				percentile_cuts = seq(0,1,percentile_width)
				percentile_times = ns_get_percentile_from_sfit(s,percentile_cuts)
				#remove dulicate times (which occur if the survival curves jumps
				#between consecutive percentiles in one step
				#percentile_times = percentile_times[which(c(1,diff(percentile_times)) > .02)]
				p2 = rep(NA,length(percentile_times));
				nn = 1;
				p2[1] = percentile_times[1];
				#browser()
				for (m in 2:length(percentile_times)){
					if (percentile_times[m]- p2[nn] > minimum_period_duration_in_days){
						nn = nn + 1;
						p2[nn] = percentile_times[m];
					}
				}
				percentile_times = p2[!is.na(p2)]
				number_at_risk = matrix(nrow=length(percentile_times)-1,ncol=1)
				num_deaths = matrix(nrow=length(percentile_times)-1,ncol=1)
				duration = matrix(nrow=length(percentile_times)-1,ncol=1)
				for (k in 2:length(percentile_times)){
					number_at_risk[k-1] = sum(current_deaths$Age.at.Death..d..Raw >= percentile_times[k-1])
					num_deaths[k-1] = sum(current_deaths$Age.at.Death..d..Raw >= percentile_times[k-1] &
							    current_deaths$Age.at.Death..d..Raw < percentile_times[k] &
							    current_deaths$Censored == 0)	
					duration[k-1] = percentile_times[k]-percentile_times[k-1]
				
				}
				#browser()
				x1 = rollmean(percentile_times,2)#[1:(length(percentile_times)-1)]
				#plot(x1,num_deaths,type='l',ylab="num_deaths");
				#plot(x1,number_at_risk/number_at_risk[1],type='l',ylab="at_risk");
				#plot(x1,duration,type='l',ylab="duration",log="y")
				#plot(x=c(1),y=c(1),type="n",log=ll, ylim=yl,xlim=xl,xlab="",ylab="hazard")
				#browser()
				vals1 = num_deaths/number_at_risk/duration;
				#if (min(current_deaths$Age.at.Death..d..Raw) < 1)
				#	browser()
				#plot(x1,vals1,log="y",type='l',col="red")
				#lines(x,vals,col="black")
				#browser()
				x = x1;
				vals = vals1;
			}
			
			
			
			time_sub_index = x >= time_bounds[1] & x < time_bounds[2];
			
			x = x[time_sub_index];
			vals = vals[time_sub_index];
			
			hazard_values_to_use_for_labels[[style_i]] = data.frame(t=x,h=vals);
			
			
			#stop("WHA")
			
			majority_time_bounds_index =  x >= majority_time_bounds[1] & x <= majority_time_bounds[2]
			
			t_offset = 0;
			if (!is.na(parametric_fits) && !is.null(parametric_fits$time_offset))
				t_offset = parametric_fits[parametric_fits[,group_column_name] == levels[[i]],]$time_offset;

			if (draw_l){
			
				#points(x,vals,col=level_colors[[i]],lty=line_styles[[i]],pch=point_type,bg=level_color,cex=ccex);
				if(thick_center_half){
					lines(x[majority_time_bounds_index]+t_offset,vals[majority_time_bounds_index],col=level_color,lty=level_lty,pch=point_type,bg=level_color,lwd=level_line_thickness);
				}
				else{
					lines(x+t_offset,vals,col=level_color,lty=level_lty,pch=point_type,bg=level_color,cex=ccex,lwd=level_line_thickness);
				}
			}
			if (draw_p){
			#	print(paste(length(x),length(vals)))
				points(x+t_offset,vals,col=level_color,lty=level_lty,pch=point_type,bg=level_color,cex=ccex);
				if (thick_center_half)
					points(x[majority_time_bounds_index]+t_offset,vals[majority_time_bounds_index],col=level_color,lty=level_lty,pch=point_type,bg=level_color,type=ttype,cex=2*ccex,);
					
			}
			
			if (length(curve_labels) > 1 || length(curve_labels)==1 && !is.na(curve_labels)){
				
				text(x[length(x)]*1.25+t_offset,vals[length(vals)]*1.25,as.double(curve_labels[[i]]),col=level_color,lty=level_lty);
			}
			hazard$period[[levels[[i]]]] = data.frame(t=x,h=vals);
		}
		

		if ("mu" %in% style){
		
			style_i = which(style=="mu");
			if (style_colors_specified && !is.na(style_colors[style_i]))
				level_color = style_colors[style_i]
			if (style_lty_specified && !is.na(style_lty[style_i]))
				level_lty = style_lty[style_i]
			if (style_line_thickness_specified && !is.na(style_line_thickness[style_i]))
				level_line_thickness = style_line_thickness[style_i]


			tmax =max(current_deaths$Age.at.Death..d..Raw)
			tmin =  min(current_deaths$Age.at.Death..d..Raw);
			f_a1 <- muhaz(current_deaths$Age.at.Death..d..Raw+x_shift, 1-current_deaths$Censored,
			b.cor="b",
			bw.method="l",
			bw.pilot= (tmax-tmin)/(8*.5*muhaz_smoothness_factor[[i]]*sum(1-current_deaths$Censored)^.02),
			min.time=tmin,max.time=tmax
			)#,bw.grid=gr)
			#stop();
			
			x= f_a1$est.grid - x_shift;#rollmean(f_a1$est.grid,2);
			
			majority_time_bounds_index =  x >= majority_time_bounds[1] & x <= majority_time_bounds[2]
			time_bounds_index =  x >= time_bounds[1] & x <= time_bounds[2]
			
			t_offset = 0;
			if (!is.na(parametric_fits) && !is.null(parametric_fits$time_offset))
				t_offset = parametric_fits[parametric_fits[,group_column_name] == levels[[i]],]$time_offset;

			lines(x[time_bounds_index]+t_offset,f_a1$haz.est[time_bounds_index],col=level_color,lty=level_lty,lwd=level_line_thickness)
			if (thick_center_half)
				lines(x[majority_time_bounds_index]+t_offset,f_a1$haz.est[majority_time_bounds_index],col=level_color,lty=level_lty,lwd=2*level_line_thickness)
			#lines(f_a1,col=level_colors[[i]],lwd=level_line_thickness)
			
			hazard[[levels[[i]]]] = data.frame(t=x[time_bounds_index],h=f_a1$haz.est[time_bounds_index]);
			
			
			hazard_values_to_use_for_labels[[style_i]] = data.frame(t=x[time_bounds_index],h=f_a1$haz.est[time_bounds_index]);
							
		}
		if (length(parametric_fits) > 0 & length(parameterization) != 0){
			cp <- parametric_fits[parametric_fits[,group_column_name] == as.character(levels[[i]]),]
			
			if (dim(cp)[1] == 0){
				browser();
				stop(paste("Could not find group",levels[[i]],"in parametric fits"));
			}
			#dl = c(min(current_deaths$Age.at.Death..d..Raw[current_deaths$Censored==0]),max(current_deaths$Age.at.Death..d..Raw[d$Censored==0]));
			x = seq(time_bounds[1],time_bounds[2],length=50);
			for (k in 1:length(parameterization)){
				if (parameterization[[k]]=="weibul"){
					m = log(cp$weibul_scale);
					s = 1/cp$weibul_shape
					haz <- dsurvreg(x, distribution='weibull',mean=m,scale=s)/(1-psurvreg(x, mean=m,scale=s))
				}
				else if (parameterization[[k]]=="cropped_weibul"){
					m = log(cp$cropped_weibul_scale);
					s = 1/cp$cropped_weibul_shape
					s = 1/cp$cropped_weibul_shape_95_conf_l
					s = 1/cp$cropped_weibul_shape_95_conf_h
					haz <- dsurvreg(x, distribution='weibull',mean=m,scale=s)/(1-psurvreg(x, mean=m,scale=s))
				}
				if (parameterization[[k]]=="inverse_gaussian"){
					haz = hinvgauss(x,mu=cp$inverse_gaussian_mu,
							  sigma2=(cp$inverse_gaussian_sigma)^2);
				}
				else if (parameterization[[k]]=="cropped_inverse_gaussian"){
					haz = hinvgauss(x,mu=cp$cropped_inverse_gaussian_mu,
							  sigma2=(cp$cropped_inverse_gaussian_sigma)^2);
				}
				else if (parameterization[[k]]=="lognormal"){
					m = cp$lognormal_intercept;
					s = cp$lognormal_scale;
				#	haz <- dlnorm(x, meanlog=m,sdlog=s)/(1-plnorm(x, meanlog=m,sdlog=s))

					haz <- dsurvreg(x, distribution='lognormal',mean=m,scale=s)/(1-psurvreg(x,distribution='lognormal', mean=m,scale=s))
					#plot(x,haz)
					#stop(haz)
				}
				else if (parameterization[[k]]=="gompertz"){
					shape = cp$gompertz_shape;
					rate = cp$gompertz_rate;
					haz = flexsurv::dgompertz(x=x,shape=shape,rate=rate)/(1-flexsurv::pgompertz(x,shape=shape,rate=rate));
				
				}else if (parameterization[[k]]=="gompertz_gamma_frailty"){
					shape = cp$gompertz_gamma_frailty_shape;
					rate = cp$gompertz_gamma_frailty_rate;
					theta = cp$gompertz_gamma_frailty_theta;
					haz = dmgompertz(x=x,shape=shape,rate=rate,theta=theta)/(1-pmgompertz(x,shape=shape,rate=rate,theta=theta));
				}else if (parameterization[[k]]=="gompertz_alt_gamma_frailty"){
					
					initial = cp$gompertz_alt_gamma_frailty_initial;
					scale = cp$gompertz_alt_gamma_frailty_scale;
					theta = cp$gompertz_alt_gamma_frailty_theta;
					#browser()
					haz = dmgompertz2(x=x,initial=initial,scale=scale,theta=theta)/(1-pmgompertz2(x,initial=initial,scale=scale,theta=theta));
				}
				else if (parameterization[[k]]=="cropped_gompertz_gamma_frailty"){
					shape = cp$cropped_gompertz_gamma_frailty_shape;
					rate = cp$cropped_gompertz_gamma_frailty_rate;
					theta = cp$cropped_gompertz_gamma_frailty_theta;
					haz = dmgompertz(x=x,shape=shape,rate=rate,theta=theta)/(1-pmgompertz(x,shape=shape,rate=rate,theta=theta));
				}else if (parameterization[[k]]=="gompertz_makeham"){
					shape = cp$gompertz_makeham_shape;
					rate = cp$gompertz_makeham_rate;
					constant = cp$gompertz_makeham_constant;
					haz = dns_makeham(x=x,shape=shape,rate=rate,makeham=constant)/(1-pns_makeham(x,shape=shape,rate=rate,makeham=constant));
				}else if (parameterization[[k]]=="cropped_gompertz_makeham"){
					shape = cp$cropped_gompertz_makeham_shape;
					rate = cp$cropped_gompertz_makeham_rate;
					constant = cp$cropped_gompertz_makeham_constant;
					haz = dns_makeham(x=x,shape=shape,rate=rate,makeham=constant)/(1-pns_makeham(x,shape=shape,rate=rate,makeham=constant));
				}
				else if (parameterization[[k]]=="cropped_gompertz"){
					shape = cp$cropped_gompertz_shape;
					rate = cp$cropped_gompertz_rate;
					haz = flexsurv::dgompertz(x=x,shape=shape,rate=rate)/(1-flexsurv::pgompertz(x,shape=shape,rate=rate));
				
				}
				else if (parameterization[[k]]=="loglogistic"){
					m = cp$loglogistic_scale;
					s = cp$loglogistic_shape;
					
					haz = hllogis(x=x, shape=s, scale = m)
				}
				else if (parameterization[[k]]=="cropped_loglogistic"){
				
					m = cp$cropped_loglogistic_scale;
					s = cp$cropped_loglogistic_shape;
					haz = hllogis(x=x, shape=s, scale = m)
				}
				else if (parameterization[[k]]=="generalized_f"){
				#print(paste("WHA",cp$generalized_f_P))
					haz = hgenf(x=x,
						mu=cp$generalized_f_mu,
						sigma=cp$generalized_f_sigma,
						Q=cp$generalized_f_Q,
						P=cp$generalized_f_P)
				}
				else if (parameterization[[k]]=="weibul_gamma_frailty"){
					m = cp$weibul_gamma_frailty_scale;
					s = cp$weibul_gamma_frailty_shape;
					h = cp$weibul_gamma_frailty_theta;
					
					#stop(paste(m,s,h))
					haz = dmllogis(x, shape=s, scale = m,theta=h)/(1-pmllogis(x, shape=s, scale = m,theta=h))
				}
				else if (parameterization[[k]]=="cropped_weibul_gamma_frailty"){
					m = cp$cropped_weibul_gamma_frailty_scale;
					s = cp$cropped_weibul_gamma_frailty_shape;
					h = cp$cropped_weibul_gamma_frailty_theta;
					haz = dmllogis(x, shape=s, scale = m,theta=h)/(1-pmllogis(x, shape=s, scale = m,theta=h))
				}
				else if (parameterization[[k]]==" "){
					warning("Skipping blank parameterization specification")
					haz=rep(NA,length(x))
				}
				else stop(paste("ns_plot_hazard() doesn't know how to handle the",parameterization[[k]],"parameterization"))
				
				hazard_values_to_use_for_labels[[k]] = data.frame(t=x,h=haz);
	
				if (style_colors_specified && !is.na(parametric_colors[[k]]))
					level_color = parametric_colors[[k]];
				if (style_lty_specified && !is.na(parametric_lty[[k]]))
					level_lty = parametric_lty[[k]];
					
				if (style_line_thickness_specified && !is.na(parametric_line_thickness))
					level_line_thickness = parametric_line_thickness[[k]]
				
				#browser()
				lwd_fit = level_line_thickness;
				if (highlight_parametric_fit_region){
					lwd_fit = 2*lwd_fit;
					lwd_unfit = lwd_fit;
				}else{
					lwd_fit = lwd_unfit=1
				}
				#browser()
				t_offset = 0;
				if (!is.null(cp$time_offset))
					t_offset = cp$time_offset;

				if(!highlight_parametric_fit_region){
					lines(x+t_offset,haz,col=level_color,lty=level_lty,lwd=lwd_fit);
				}
				else{
					if (length(grep("crop",parameterization[[k]]))==1){
						fitted_range = c(cp$cropped_parametric_fit_lower_bound,cp$cropped_parametric_fit_upper_bound);
					}
					else{

						fitted_range = c(cp$parametric_fit_lower_bound,cp$parametric_fit_upper_bound);
					}
					fit_index = x >= fitted_range[1] & x <= fitted_range[2];
					unfit_index_l = which(x <= fitted_range[1]);
					unfit_index_h = which(x >= fitted_range[2]);
					if (length(unfit_index_l) > 0 & length(fit_index)>0)
						unfit_index_l = c(unfit_index_l,which(fit_index)[1]);
					if (length(unfit_index_h) > 0 & length(fit_index)>0)
						unfit_index_h = c(which(fit_index)[length(which(fit_index))],unfit_index_h);
					
					#print(unfit_index_h)
					lines(x[fit_index]+t_offset,haz[fit_index],col=level_color,lty=level_lty,lwd=lwd_fit)
					lines(x[unfit_index_l]+t_offset,haz[unfit_index_l],col=level_color,lty=level_lty,lwd=lwd_unfit)
					lines(x[unfit_index_h]+t_offset,haz[unfit_index_h],col=level_color,lty=level_lty,lwd=lwd_unfit)
					
					
				}
			}
		}
	
		if (length(quantiles_to_label) != 0){
			for (i in 1:length(style)){
				if (is.na(quantile_label_style[[i]]))
					next;

				h = hazard_values_to_use_for_labels[[i]]

				ssfit = survfit(Surv(current_deaths$Age.at.Death..d..Raw,1-current_deaths$Censored)~1);
				quantile_times_to_plot = ns_get_percentile_from_sfit(ssfit,quantiles_to_label)
			#	browser()
				hazards_to_plot = rep(NA,length(quantile_times_to_plot))
				h <- h[order(h$t),];
				#find two calculated hazards in between which the requested hazard falls,
				#and find the hazard as a linear interpolation between them
				for (j in 1:length(quantile_times_to_plot)){
					hh_u = min(which(h$t >= quantile_times_to_plot[j]));
					if (hh_u > 1){
						hh_l = hh_u-1;
					}else hh_l = hh_u;

					t_frac = (quantile_times_to_plot[j] - h$t[hh_l])/(h$t[hh_u] - h$t[hh_l]);
					hazards_to_plot[[j]] = h$h[hh_l] + (h$h[hh_u]-h$h[hh_l])*t_frac;

				}
				#print(quantile_times_to_plot);
				#print(hazards_to_plot);
				#stop()
				faded_col = col2rgb(level_color)/255
				faded_col=rgb(faded_col[1],faded_col[2],faded_col[3],.5)
				if (quantile_label_style[[i]] == "lines"){
					for (j in 1:length(quantile_times_to_plot)){
						lines(c(quantile_times_to_plot[j],quantile_times_to_plot[j],xl[1])+t_offset,
						      c(yl[1],hazards_to_plot[j],hazards_to_plot[j]),col=faded_col);
					}
				}
				if (quantile_label_style[[i]] == "v_lines"){
									for (j in 1:length(quantile_times_to_plot)){
										lines(c(quantile_times_to_plot[j],quantile_times_to_plot[j])+t_offset,
										      c(yl[1],hazards_to_plot[j]),col=faded_col);
									}
				}
				if (is.double(quantile_label_style[[i]])){
					for (j in 1:length(quantile_times_to_plot)){
						points(quantile_times_to_plot[j],
						       hazards_to_plot[j],
						       pch=quantile_point_pch,
						       cex=quantile_label_style[[i]],
						       col=level_color);
					}
				}
			}
		}
	}
	if(parameterization_legend){
	 
		if (length(style_colors) != 0){
			colors = parametric_colors;
		}else colors = rep("#000000",length(parameterization))
		if (length(style_lty) != 0){
			lty = parametric_lty;
		}else lty = rep(1,length(parameterization))


		legend(legend_position,legend = parameterization,lty=lty,col=colors,cex=.75,bty='n');
	}
	return(hazard);
}
	

ns_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#from web somewhere
interleave <- function(v1,v2){
	ord1 <- 2*(1:length(v1))-1
	ord2 <- 2*(1:length(v2))
	return (c(v1,v2)[order(c(ord1,ord2))])
}

#' @export
ns_plot_colored_hazard_by_groups <- function(deaths,group_column_name,
					color_scheme="rainbow",
					color_column=NA,
					colors = NA,
					second_style_column=NA,
					shuffle_colors=F,...){
	
	group_levels = as.character(unique(deaths[,group_column_name]));
	color_levels = as.character(group_levels);
	#get unique number for each color and line
	if (!is.na(color_column)){
		color_levels = as.character(unique(deaths[,color_column]));
		color_id_lookup = list();

		for (i in 1:length(color_levels)){
			color_id_lookup[[as.character(color_levels[i])]] = i;
		}
	}
	if (!is.na(second_style_column)){
		second_style_levels = unique(deaths[,second_style_column]);
		second_style_lookup = list();
		
		for (i in 1:length(second_style_levels))
			second_style_lookup[[second_style_levels[i]]] = i;
	}
	#if rainbow, the second style is the line style
	if(!is.na(colors)){
		if (length(color_levels) != length(colors))
			stop(paste(length(colors),"colors specified for",length(color_levels),"groups"));
		colors_to_use = colors
		
		line_styles_to_use=matrix(nrow=length(group_levels),ncol=1);
				for (k in 1:length(group_levels)){
					indx = deaths[,group_column_name]==group_levels[[k]];
					if (!is.na(second_style_column)){
						dominant_second_style_id = second_style_lookup[[ns_mode(deaths[indx,second_style_column])]];
						line_styles_to_use[[k]] = (dominant_second_style_id-1)*1+1;
					}
					else line_styles_to_use[[k]] = 1;
		}
	}
	else if (color_scheme=="rainbow"){
		colors_to_use = rainbow(length(color_levels),s=1,v=.8)
		
		line_styles_to_use=matrix(nrow=length(group_levels),ncol=1);
		for (k in 1:length(group_levels)){
			indx = deaths[,group_column_name]==group_levels[[k]];
			if (!is.na(second_style_column)){
				dominant_second_style_id = second_style_lookup[[ns_mode(deaths[indx,second_style_column])]];
				line_styles_to_use[[k]] = (dominant_second_style_id-1)*1+1;
			}
			else line_styles_to_use[[k]] = 1;
		}
	}
	
	#if ramp, the second style is a different part of the color spectrum
	else if (color_scheme == "ramp"){
		gradient_functions = list();
		gradient_functions[[1]] = (colorRampPalette(c("black", "red","orange"), space = "rgb"));
		gradient_functions[[2]] = (colorRampPalette(c("gray", "blue","green"), space = "rgb"));
		gradient_functions[[3]] = (colorRampPalette(c("purple", "pink","white"), space = "rgb"));
		if (!is.na(second_style_column) && length(second_style_levels) > 3)
			stop(paste("Too many (", length(second_style_levels), ") secondary levels found in column ", second_style_column));
		
		line_styles_to_use=matrix(nrow=length(group_levels),ncol=1);
		colors_to_use=matrix(nrow=length(group_levels),ncol=1);
		for (k in 1:length(group_levels)){
			line_styles_to_use[[k]] = 1;
			indx = deaths[,group_column_name]==group_levels[[k]];
			dominant_color_id =     color_id_lookup[[as.character(ns_mode(deaths[indx,color_column]))]];
			
			#pull the right color from the correct gradient
			if (!is.na(second_style_column)){
				dominant_second_style_id = second_style_lookup[[ns_mode(deaths[indx,second_style_column])]];
			}
			else dominant_second_style_id = 1;
			
			colors_to_use[k] = (gradient_functions[[dominant_second_style_id]])(length(color_levels))[dominant_color_id];
		}
	}
	else stop(paste("Unknown color scheme: ", color_scheme));
	if (shuffle_colors)
	colors_to_use = colors_to_use[sample.int(length(colors_to_use),length(colors_to_use),replace=F)]

	hazard = ns_plot_hazard_by_groups(deaths=deaths,group_column_name=group_column_name,
				 level_colors=colors_to_use,...);
		 
	return(list(colors=colors_to_use,line_styles=line_styles_to_use,hazard = hazard));
}


#' @export
ns_plot_colored_hazard_by_groups_with_legend<- function(deaths,group_column_name,legend_position="topleft",legend_title="",...){
	color_and_styles = 
		ns_plot_colored_hazard_by_groups(deaths,group_column_name=group_column_name,...);
		
	ss = 1:length(color_and_styles[["line_styles"]])
	
	levels = unique(as.character(deaths[,group_column_name]));
	levels = levels[order(levels)]
	legend(legend_position,legend = levels,title=legend_title,lty= color_and_styles[["line_styles"]],col= color_and_styles[["colors"]],cex=.75,
		bg="#FFFFFF")
	return(color_and_styles);
}