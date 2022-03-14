
library("rms");
library("flexsurv");
library("STAR");
#detach("package:eha", unload=TRUE)
library("eha");
library("gsl")

ns_pweiner_diff= function(t,b){
	return( 
		( b/sqrt(2*(pi*t)^3) )*
			exp((b^2)/(2*t))*
	   gsl::gamma_inc(0,(b^2)/(2*t))
	)
}

ns_invguass_boundary_parameters = function(mu,sigma2,drift_rate){
		boundary = drift_rate*mu;
		new_sigma2 = sigma2*boundary^2;
	return(data.frame(
		drift_rate=drift_rate,
		sigma2 = new_sigma2,
		boundary = boundary
		)
	)
}

ns_loglogistic_alpha_weibull_equivalent <- function(t,shape,scale){
	return( shape -  ( shape*( (t/scale)^shape ) )/( 1+ (t/scale)^shape ) );
}

if (0){
	#from OO Aalen Frailty 1994
	dns_aalen_frailty_drift = function(x, mu , sigma2, theta){
		c = 1;
		t = x;
		return(  
			(c/sqrt(2*pi)) * (1/(t*sqrt(t*t*sigma2+t) )) * exp( -(c-mu*t)^2/ (2*(t*t*mu*mu+t)) )
		      )
	}
	pns_aalen_frailty_drift = function(q, mu , sigma2, theta){
		c = 1;
		t = q;
		return(  
			pnorm( (c-mu*t)/sqrt(t*t*sigma2+t)) - exp(2*c*mu+s*c*c*sigma2)*pnorm( (-c-2*c*t*sigma2-mu*t)/sqrt(t*t*sigma2+t) )
		      )
	}
	hns_aalen_frailty_drift = function(x,mu,sigma2,theta){
		return(
			dns_aalen_frailty_drift(x,mu,sigma2,theta)
				/
				pns_aalen_frailty_drift(x,mu,sigma2,theta)
		)
	}
}



#STAR::pinvgauss modified to correctly handle very small sigma2
ns_pinvgauss =cmpfun(function (q, mu = 1, sigma2 = 1, boundary = NULL, lower.tail = TRUE, log.p = FALSE) {
	#print(paste(mu,sigma2));
    if (any(q < 0)) 
        stop("q must contain positive values")
    if (any(mu <= 0)) 
   	 return(NaN);#stop("mu must be positive")
    if (all(!is.null(sigma2))) 
        if (any(sigma2 <= 0)) 
            return(NaN);#stop("sigma2 must be positive")
    if (all(!is.null(boundary))) 
        if (any(boundary <= 0)) 
            stop("boundary must be positive")
    if (all(!is.null(sigma2)) && all(!is.null(boundary))) 
        stop("One of sigma2 or boundary must be specified, not both")
    if (all(!is.null(boundary))) {
        sigma2 <- (1/boundary)^2
        mu <- boundary * mu
    }
    t <- q/mu
    v <- sqrt(q * sigma2)
    if (lower.tail & !log.p) 
        return(
		pnorm((t - 1)/v)+ 
		exp(   
			2/(mu * sigma2)+ pnorm(-(t + 1)/v,log.p=T)
		)
        )
    if (!lower.tail & !log.p) 
        return(1 - (
        	pnorm((t - 1)/v) + 
        	exp( 
        		2/(mu * sigma2) +  pnorm(-(t + 1)/v,log.p=T)
        	   )
        	   )
        	)
    if (lower.tail & log.p) 
        return(
        	log(
       			pnorm((t - 1)/v) +
       			exp( 
       				2/(mu * sigma2) + pnorm(-(t + 1)/v,log.p=T)
       			)
       		)
       		)
    if (!lower.tail & log.p) 
        return(
        	log(
        		1 - (
        		pnorm((t - 1)/v) + 
        		exp(
        			2/(mu * sigma2) + pnorm(-(t + 1)/v,log.p=T)
        		)
        		)
        	)
        )
})
pns_invgauss = ns_pinvgauss;

dns_invgauss = cmpfun(function (x, mu = 1, sigma2 = 1, boundary = NULL, log = FALSE) 
{
    if (any(x <= 0)) 
        stop("y must contain positive values")
    if (any(mu <= 0)) 
        return(NaN);#stop("mu must be positive")
    if (all(!is.null(sigma2))) 
        if (any(sigma2 <= 0)) 
            return(NaN);#stop("sigma2 must be positive")
    if (all(!is.null(boundary))) 
        if (any(boundary <= 0)) 
            return(NaN);#stop("boundary must be positive")
    if (all(!is.null(sigma2)) && all(!is.null(boundary))) 
        stop("One of sigma2 or boundary must be specified, not both")
    if (all(!is.null(boundary))) {
        sigma2 <- (1/boundary)^2
        mu <- boundary * mu
    }
    tmp <- -(x - mu)^2/(2 * x * sigma2 * mu^2) - (log(2 * pi * 
        sigma2) + 3 * log(x))/2
    if (!log) 
        tmp <- exp(tmp)
    tmp
})

qns_invgauss = cmpfun(function(p, mu = 1, sigma2 = 1, boundary = NULL) 
{
    if (any(p < 0 | p > 1)) 
        stop("p must lie between 0 and 1")
    if (any(mu <= 0)) 
          return(NaN);#stop("mu must be positive")
    if (all(!is.null(sigma2))) 
        if (any(sigma2 <= 0)) 
              return(NaN);#stop("sigma2 must be positive")
    if (all(!is.null(boundary))) 
        if (any(boundary <= 0)) 
              return(NaN);#stop("boundary must be positive")
    if (all(!is.null(sigma2)) && all(!is.null(boundary))) 
        stop("One of sigma2 or boundary must be specified, not both")
    if (all(!is.null(boundary))) {
        sigma2 <- (1/boundary)^2
        mu <- boundary * mu
    }
    len <- max(length(p), length(mu), length(sigma2))
    if (length(p) != len) {
        if (length(p) == 1) 
            p <- rep(p, len)
        else stop("length of p incorrect")
    }
    if (length(mu) != len) {
        if (length(mu) == 1) 
            mu <- rep(mu, len)
        else stop("length of m incorrect")
    }
    if (length(sigma2) != len) {
        if (length(sigma2) == 1) 
            sigma2 <- rep(sigma2, len)
        else stop("length of sigma2 incorrect")
    }
    theta <- 1/mu/sigma2
    approx <- mu * exp(qnorm(p) * sqrt(1/theta) - 0.5/theta)
    sapply(1:len, function(idx) {
        if (identical(p[idx], 0)) 
            return(0)
        if (identical(p[idx], 1)) 
            return(Inf)
        interval <- approx[idx] * c(0.95, 1.05)
        h <- function(q) pinvgauss(q, mu[idx], sigma2[idx]) - 
            p[idx]
        while (h(interval[1]) * h(interval[2]) > 0) interval <- interval * 
            c(0.9, 1.1)
        uniroot(h, interval)$root
    })
})

custom.ns_invgauss <- list(name="ns_invgauss",
	pars=c("mu","sigma2"),
	location="mu",
	transforms=c(base::log, base::log),
	inv.transforms=c(base::exp, base::exp),
	inits=function(t){ c(stats::median(t),.001) }
)


#modified from STAR::invgaussMLE

ns_invgaussMLE = function (yi, ni = numeric(length(yi)) + 1, si = numeric(length(yi)) + 1, parameterization = "sigma2",inits=NA) 
{
    if (inherits(yi, "spikeTrain")) 
        yi <- diff(yi)
    if (any(yi < 0)) 
        stop("yi elements must be non-negative")
    yi <- as.numeric(yi)
    if (!identical(length(yi), length(ni))) 
        stop("yi and ni should have the same length")
    if (any(ni < 0)) 
        stop("ni elements must be non-negative")
    if (!identical(class(ni), "integer") && !identical(ni, round(ni))) 
        stop("ni should be a vector of positive integers")
    if (!identical(length(yi), length(si))) 
        stop("yi and si should have the same length")
    if (any(si < 0)) 
        stop("si elements must be non-negative")
    if (!identical(si, round(si))) 
        stop("si should be a vector of positive integers")
    if (any(si > ni)) 
        stop("si elements should not be greater than ni elements")
   
    
    start.time <- Sys.time()
    s.dot <- sum(si)
    n.dot <- sum(ni)
    ci <- ni - si
    c.dot <- sum(ci)
    
     minusLogLik <- cmpfun(function(p) {
            if (missing(p)) {
                txt <- paste("This function argument should be a 2 component vector:\n", 
                    "  component 1 is the log of the mu parameter,\n", 
                    "  component 2 is the log of the sigma2 parameter,\n", 
                    "using the 'sigma2' parameterization of the inverse Gaussian.\n")
                cat(txt)
            }
            else {
                mu <- exp(p[1])
                sigma2 <- exp(p[2])
              #  browser()
                 -(ifelse(s.dot > 0, sum(dns_invgauss(yi[si > 0], mu, 
                    sigma2, log = TRUE) * si[si > 0]), 0) + ifelse(c.dot > 
                    0, sum(ns_pinvgauss(yi[ci > 0], mu, sigma2, lower.tail = FALSE, 
                    log.p = TRUE) * ci[ci > 0]), 0))
           #    print(paste(mu,sigma2,r));
              # if (abs(r)<.0000001)
              # 	browser()
               #return(r);
            }
    })
    
    
    if (s.dot == n.dot) {
    	if (is.na(inits)){
		mu.hat <- weighted.mean(yi, ni)
		inv.mean <- weighted.mean(1/yi, ni)
		sigma2.hat <- inv.mean - 1/mu.hat
	}else{
		mu.hat = inits[1];
		sigma2.hat = init[2];
	}
        estimate <- c(mu.hat, sigma2.hat)
        observedI <- matrix(c(n.dot/(sigma2.hat * mu.hat^3), 
            0, 0, n.dot/sigma2.hat^2/2), nrow = 2, byrow = TRUE)
        se <- sqrt(diag(solve(observedI)))
        l <- -minusLogLik(log(estimate))
    }
    else {
	    if (is.na(inits)){
		if (s.dot >= 10) {
		    mu.hat <- weighted.mean(yi, si)
		    inv.mean <- weighted.mean(1/yi, si)
		    sigma2.hat <- inv.mean - 1/mu.hat
		}
		else {
		    mu.hat <- weighted.mean(yi, ni)
		    inv.mean <- weighted.mean(1/yi, ni)
		    sigma2.hat <- inv.mean - 1/mu.hat
		}
	}else{
			mu.hat = inits[1];
			sigma2.hat = init[2];
	}
        mleFit <- optim(par = log(c(mu.hat, sigma2.hat)), fn = minusLogLik, 
            
            method = "SANN", hessian = TRUE,
            control=list(maxit=400,tmax=10,temp=1,trace=T,REPORT=10))
        estimate <- exp(mleFit$par)
        newVar <- (1/estimate) %o% (1/estimate)
        observedI <- mleFit$hessian * newVar
        se <- sqrt(diag(solve(observedI)))
        l <- -mleFit$value
    }
    if (parameterization == "sigma2") {
        names(estimate) <- c("mu", "sigma2")
        names(se) <- c("mu", "sigma2")
        rFct <- function(mu, sigma2) -minusLogLik(log(c(mu, sigma2))) - 
            l
    }
    else {
        boundary.hat <- (sigma2.hat)^(-0.5)
        mu.hat <- mu.hat/boundary.hat
        if (s.dot == n.dot) {
            estimate <- c(mu.hat, boundary.hat)
            observedI <- observedI * matrix(c(boundary.hat^2, 
                -2/boundary.hat^2, -2/boundary.hat^2, 4/boundary.hat^6), 
                nrow = 2, byrow = TRUE)
            se <- se * c(1/boundary.hat, boundary.hat^3/2)
        }
        else {
            estimate <- c(estimate[1] * sqrt(estimate[2]), 1/sqrt(estimate[2]))
            newVar <- newVar * (c(estimate[2], -2/estimate[2]^3) %o% 
                c(estimate[2], -2/estimate[2]^3))
            observedI <- mleFit$hessian * newVar
            se <- sqrt(diag(solve(observedI)))
        }
        names(estimate) <- c("mu", "boundary")
        names(se) <- c("mu", "boundary")
        rFct <- function(mu, boundary) -minusLogLik(log(c(mu * 
            boundary, 1/boundary^2))) - l
    }
    
     stop.time <- Sys.time()    
    result <- list(estimate = estimate, se = se, logLik = l, 
        r = rFct, mll = minusLogLik, call = match.call(),start.time=start.time,stop.time=stop.time)
    class(result) <- "durationFit"
    return(result)
}

ns_fit_invgauss <-function(surv_obj,inits=NA){
	#d = data.frame(t=surv_obj[,1],cc=surv_obj[,2]);
	
	#invg2 = invGauss::invGauss(formula.mu=Surv(t,cc)~1,data=d)
	invg = ns_invgaussMLE(yi=surv_obj[,1],
				  si=surv_obj[,2],
				  parameterization="sigma2",inits=NA);

	AIC = -2*invg$logLik + 2*2
	upper = invg$estimate+invg$se
	lower = invg$estimate-invg$se
	#browser()
	return(list( mu = invg$estimate[1],sigma2 = invg$estimate[2],AIC=AIC,upper=upper,lower=lower))
}

ns_fit_parameterization_group = function(formula,parametric_fits){
	ff = terms(formula);
	rr = eval.parent(attr(ff,"variables"))
	pp = attr(ff,"factors")
	grouping = rr[[which(pp==1)]]
	deaths = rr[[which(pp!=1)]]
	if (class(deaths) != "Surv")
		stop("Must provide a Surv object");
	
	fits = NULL;
	successful_fits = list()
	num_groups = length(unique(grouping))
	#contrasts(bj_group) <- contr.treatment(num_groups,base=reference_factor_level);
	for (g in unique(grouping)){
		dd = deaths[grouping==g,];
		par_fit = ns_fit_parameterization(dd,c("inverse_gaussian"),death_column=1,censored_column=2,rev_censored_group=F);
	#	browser()
		par_fit$fits$group = g;
		fits = rbind(par_fit$fits,fits);
		successful_fits[[as.character(g)]] = par_fit$successful_fits;
	}
	#browser()
	return(list(fits=fits,successful_fits=successful_fits));
}

ns_fit_parameterization <- function(deaths,parametric_fits=NA,death_column="Age.at.Death..d..Raw",censored_column="Censored",rev_censored_group=T){

	lognormal_scale <- NA;
	lognormal_intercept <- NA;
	lognormal_AIC <- NA;

	weibull_shape <- NA
	weibull_scale <- NA;
	weibull_shape_95_conf_l <- NA;
	weibull_shape_95_conf_h <- NA;
	weibull_scale_95_conf_l <- NA;
	weibull_scale_95_conf_h <- NA;
	weibull_AIC <- NA;

	weibull_gamma_frailty_shape <- NA;
	weibull_gamma_frailty_scale <- NA;
	weibull_gamma_frailty_theta <- NA;
	weibull_gamma_frailty_shape_95_conf_l <- NA;
	weibull_gamma_frailty_shape_95_conf_h <- NA;
	weibull_gamma_frailty_scale_95_conf_l <- NA;
	weibull_gamma_frailty_scale_95_conf_h <- NA;
	weibull_gamma_frailty_theta_95_conf_l <- NA;
	weibull_gamma_frailty_theta_95_conf_h <- NA;
	weibull_gamma_frailty_AIC <- NA;

	inverse_gaussian_mu <- NA;
	inverse_gaussian_sigma2 <- NA;
	inverse_gaussian_mu_95_conf_l <- NA;
	inverse_gaussian_mu_95_conf_h <- NA;
	inverse_gaussian_sigma_95_conf_l <- NA;
	inverse_gaussian_sigma_95_conf_h <- NA;
	inverse_gaussian_AIC = NA

	gompertz_gamma_frailty_shape <- NA;
	gompertz_gamma_frailty_rate <- NA;
	gompertz_gamma_frailty_theta <- NA;
	gompertz_gamma_frailty_shape_95_conf_l <- NA;
	gompertz_gamma_frailty_shape_95_conf_h <- NA;
	gompertz_gamma_frailty_rate_95_conf_l <- NA;
	gompertz_gamma_frailty_rate_95_conf_h <- NA;
	gompertz_gamma_frailty_theta_95_conf_l <- NA;
	gompertz_gamma_frailty_theta_95_conf_h <- NA;
	gompertz_gamma_frailty_AIC <- NA;

	gompertz_alt_gamma_frailty_initial <- NA;
	gompertz_alt_gamma_frailty_scale <- NA;
	gompertz_alt_gamma_frailty_theta <- NA;
	gompertz_alt_gamma_frailty_initial_95_conf_l <- NA;
	gompertz_alt_gamma_frailty_initial_95_conf_h <- NA;
	gompertz_alt_gamma_frailty_scale_95_conf_l <- NA;
	gompertz_alt_gamma_frailty_scale_95_conf_h <- NA;
	gompertz_alt_gamma_frailty_theta_95_conf_l <- NA;
	gompertz_alt_gamma_frailty_theta_95_conf_h <- NA;
	gompertz_alt_gamma_frailty_AIC <- NA;

	gompertz_makeham_shape <- NA;
	gompertz_makeham_rate <- NA;
	gompertz_makeham_constant <- NA;
	gompertz_makeham_shape_95_conf_l <- NA;
	gompertz_makeham_shape_95_conf_h <- NA;
	gompertz_makeham_rate_95_conf_l <- NA;
	gompertz_makeham_rate_95_conf_h <- NA;
	gompertz_makeham_constant_95_conf_l <- NA;
	gompertz_makeham_constant_95_conf_h <- NA;
	gompertz_makeham_AIC <- NA;

	loglogistic_shape <-NA;
	loglogistic_scale <-NA;
	loglogistic_shape_95_conf_l <- NA;
	loglogistic_shape_95_conf_h <- NA;
	loglogistic_scale_95_conf_l <- NA;
	loglogistic_scale_95_conf_h <- NA;
	loglogistic_AIC <- NA;

	loglogistic_t_slope_zero <-NA;
	loglogistic_weibull_equivalent_alpha_percentile_25<-NA;
	loglogistic_weibull_equivalent_alpha_percentile_50<-NA;
	loglogistic_weibull_equivalent_alpha_percentile_75<-NA;

	gompertz_shape <- NA;
	gompertz_rate <- NA;
	gompertz_shape_95_conf_l <- NA;
	gompertz_shape_95_conf_h <- NA;
	gompertz_rate_95_conf_l <- NA;
	gompertz_rate_95_conf_h <- NA;
	gompertz_AIC <-  NA;

	generalized_f_mu <- NA;
	generalized_f_sigma <- NA;
	generalized_f_Q <- NA;
	generalized_f_P <- NA;
	generalized_f_AIC <- NA;
	
	generalized_extreme_value_shape = NA;
	generalized_extreme_value_scale = NA;
	generalized_extreme_value_loc = NA;
	generalized_extreme_value_AIC = NA;

	if (!(is.numeric(death_column) && death_column <= dim(deaths)[1]) && ! (death_column %in% names(deaths)))
		stop(paste("Death-time column",death_column, "does not exist in data"));
	if (rev_censored_group){
		survival_object <<- Surv(deaths[,death_column],1-deaths[,censored_column]);
	}else{
		survival_object <<- Surv(deaths[,death_column],deaths[,censored_column]);
	}
	
	rrr = range(deaths[,death_column])
	parametric_fit_lower_bound = rrr[1];
	parametric_fit_upper_bound = rrr[2];
  
  	lognormal_regression_model =
	weibull_regression_model =
	weibull_gamma_frailty_regression_model =
	gompertz_gamma_frailty_regression_model =
	gompertz_alt_gamma_frailty_regression_model =
	gompertz_makeham_regression_model =
	gompertz_regression_model =
	generalized_extreme_value_regression_model =
	loglogistic_regression_model =NA;
	successful_fits = NULL;

	if ("lognormal" %in% parametric_fits){
	tryCatch(
		{	
			model_loglog <- flexsurvreg(survival_object~1,dist="lnorm");
			lognormal_regression_model <- model_loglog

			lognormal_intercept = model_loglog$res[1]
			lognormal_scale = model_loglog$res[2]
			lognormal_AIC = model_loglog[["AIC"]]
			successful_fits = c("lognormal",successful_fits)
		}
			,
			error = function(e){
				warning(e);}
		);
	}
			#browser()
	#print(length(lognormal_regression_model));print("<-ZO");
	if ("weibull" %in% parametric_fits){
		tryCatch(
		{
			
			model_weibull <- flexsurvreg(survival_object~1,dist="weibull");
			weibull_regression_model <- model_weibull;

			weibull_shape = model_weibull$res[1]
			weibull_scale = model_weibull$res[2]
			weibull_shape_95_conf_l= model_weibull$res[1,2]
			weibull_shape_95_conf_h = model_weibull$res[1,3]
			weibull_scale_95_conf_l = model_weibull$res[2,2]
			weibull_scale_95_conf_h = model_weibull$res[2,3]


			weibull_AIC = model_weibull[["AIC"]]
			successful_fits = c("weibull",successful_fits)
			#browser()
		},
		error = function(e) {print(paste("Data length:",dim(deaths)[1]));

		warning(paste("Could not complete Weibul Fit:",e$message))
		#plot(sfit);
		}
	);
	}
	if ("generalized_extreme_value" %in% parametric_fits){
	 library("evd")

	tryCatch(
		{
			if (any((deaths[,censored_column]==1) == rev_censored_group))
				stop("generalized extreme value model fitting not implemented for censored data");
			#browser();
			model_gev <- fgev(deaths[,death_column],std.err=F,control=list("maxit"=1000));
			#browser()
			generalized_extreme_value_regression_model <- model_gev;
			generalized_extreme_value_shape = model_gev$estimate["shape"]
			generalized_extreme_value_scale = model_gev$estimate["scale"]
			generalized_extreme_value_loc = model_gev$estimate["loc"]
			generalized_extreme_value_AIC = AIC(model_gev);
			successful_fits = c("generalized_extreme_value",successful_fits)
		#	browser()

		},
		error = function(e) {print(paste("Data length:",deaths[,death_column]));
		warning(paste("Could not complete Generalized Extreme value Fit:",e$message))}
	);
		
		
		
	}
	if ("weibull_gamma_frailty" %in% parametric_fits){
		tryCatch(
		{
			model_mllogis <- flexsurvreg(survival_object~1,dist=custom.mllogis);
			weibull_gamma_frailty_regression_model <- model_mllogis;
			weibull_gamma_frailty_shape = model_mllogis$res[1]
			weibull_gamma_frailty_scale = model_mllogis$res[2]
			weibull_gamma_frailty_theta = model_mllogis$res[3]
			weibull_gamma_frailty_shape_95_conf_l = model_mllogis$res[1,2]
			weibull_gamma_frailty_shape_95_conf_h = model_mllogis$res[1,3]
			weibull_gamma_frailty_scale_95_conf_l = model_mllogis$res[2,2]
			weibull_gamma_frailty_scale_95_conf_h = model_mllogis$res[2,3]
			weibull_gamma_frailty_theta_95_conf_l = model_mllogis$res[3,2]
			weibull_gamma_frailty_theta_95_conf_h = model_mllogis$res[3,3]
			weibull_gamma_frailty_AIC = model_mllogis[["AIC"]]
			successful_fits = c("weibull_gamma_frailty",successful_fits)

		},
		error = function(e) {print(paste("Data length:",dim(deaths)[1]));
		warning(paste("Could not complete Weibul with Gamma Frailty Fit:",e$message))}
		);
	}
	if ("inverse_gaussian" %in% parametric_fits){
		tryCatch(
		{	
			#browser()
			ig <- ns_fit_invgauss(survival_object)
			inverse_gaussian_mu = ig$mu
			inverse_gaussian_sigma2 = ig$sigma2;
			inverse_gaussian_mu_95_conf_l = ig$upper[1]
			inverse_gaussian_mu_95_conf_h = ig$lower[1]
			inverse_gaussian_sigma_95_conf_l = ig$upper[2]
			inverse_gaussian_sigma_95_conf_h = ig$lower[2]
			inverse_gaussian_AIC = ig$AIC
			successful_fits = c("inverse_gaussian",successful_fits)
		},
		error = function(e) {stop(e);
		warning(paste("Could not complete inverse gaussian fit :",e$message))}
		);
	}
	if (("gompertz_gamma_frailty" %in% parametric_fits)){
		tryCatch(
			{
			model_mgompertz <- flexsurvreg(survival_object~1,dist=custom.mgompertz);
			gompertz_gamma_frailty_regression_model <- model_mgompertz;

			gompertz_gamma_frailty_shape = model_mgompertz$res[1]
			gompertz_gamma_frailty_rate = model_mgompertz$res[2]
			gompertz_gamma_frailty_theta = model_mgompertz$res[3]
			gompertz_gamma_frailty_shape_95_conf_l = model_mgompertz$res[1,2]
			gompertz_gamma_frailty_shape_95_conf_h = model_mgompertz$res[1,3]
			gompertz_gamma_frailty_rate_95_conf_l = model_mgompertz$res[2,2]
			gompertz_gamma_frailty_rate_95_conf_h = model_mgompertz$res[2,3]
			gompertz_gamma_frailty_theta_95_conf_l = model_mgompertz$res[3,2]
			gompertz_gamma_frailty_theta_95_conf_h = model_mgompertz$res[3,3]
			gompertz_gamma_frailty_AIC = model_mgompertz[["AIC"]]
			successful_fits = c("gompertz_gamma_frailty",successful_fits)

			},
			error = function(e) {print(paste("Data length:",dim(deaths)[1]));
			warning(paste("Could not complete Gompertz with Gamma Frailty Fit:",e$message))}
		);
	}
	if(("gompertz_alt_gamma_frailty" %in% parametric_fits)){
		tryCatch(
			{
			#print("start");
			model_mgompertz_alt <- flexsurvreg(survival_object~1,dist=custom.mgompertz2#,
							#	method="Nelder-Mead"#,
							#	control=list(trace=99)
								);
						#		stop()
			gompertz_alt_gamma_frailty_regression_model <- model_mgompertz_alt;

			gompertz_alt_gamma_frailty_initial = model_mgompertz_alt$res[1]
			gompertz_alt_gamma_frailty_scale = model_mgompertz_alt$res[2]
			gompertz_alt_gamma_frailty_theta = model_mgompertz_alt$res[3]
			gompertz_alt_gamma_frailty_initial_95_conf_l = model_mgompertz_alt$res[1,2]
			gompertz_alt_gamma_frailty_initial_95_conf_h = model_mgompertz_alt$res[1,3]
			gompertz_alt_gamma_frailty_scale_95_conf_l = model_mgompertz_alt$res[2,2]
			gompertz_alt_gamma_frailty_scale_95_conf_h = model_mgompertz_alt$res[2,3]
			gompertz_alt_gamma_frailty_theta_95_conf_l = model_mgompertz_alt$res[3,2]
			gompertz_alt_gamma_frailty_theta_95_conf_h = model_mgompertz_alt$res[3,3]
			gompertz_alt_gamma_frailty_AIC = model_mgompertz_alt[["AIC"]]
			
			successful_fits = c("gompertz_alt_gamma_frailty",successful_fits)
			#browser()
			},
			error = function(e) {print(paste("Data length:",dim(deaths)[1]));
			warning(paste("Could not complete Gompertz Alternate Parameterization with Gamma Frailty Fit:",e$message))}
		);
	}
	if (("gompertz_makeham" %in% parametric_fits)){
		#tryCatch(
		#	{
			browser()
			model_makeham <- flexsurvreg(survival_object~1,dist=custom.makeham);
			gompertz_makeham_regression_model <- model_makeham;

			gompertz_makeham_shape = model_makeham$res[1]
			gompertz_makeham_rate = model_makeham$res[2]
			gompertz_makeham_constant = model_makeham$res[3]
			gompertz_makeham_shape_95_conf_l = model_makeham$res[1,2]
			gompertz_makeham_shape_95_conf_h = model_makeham$res[1,3]
			gompertz_makeham_rate_95_conf_l = model_makeham$res[2,2]
			gompertz_makeham_rate_95_conf_h = model_makeham$res[2,3]
			gompertz_makeham_constant_95_conf_l = model_makeham$res[3,2]
			gompertz_makeham_constant_95_conf_h = model_makeham$res[3,3]
			gompertz_makeham_AIC = model_makeham[["AIC"]]
			
			successful_fits = c("gompertz_makeham",successful_fits)

		#	},
		#	error = function(e) {print(paste("Data length:",dim(deaths)[1]));
		#	warning(paste("Could not complete Gompertz Makeham:",e$message))}
		#);
	}
	#generalized_f_distribution
	if (("generalized_f" %in% parametric_fits)){
		tryCatch(
			{

			model_f <- flexsurvreg(survival_object~1,dist="genf");
			generalized_f_regression_model <- model_f;
			generalized_f_mu = model_f$res[1]
			generalized_f_sigma = model_f$res[2]
			generalized_f_Q = model_f$res[3]
			generalized_f_P = model_f$res[4]
			generalized_f_AIC = model_f[["AIC"]]
			successful_fits = c("generalized_f",successful_fits)

			},
			error = function(e) {
				warning(paste("Could not complete Generalized F Fit:",e$message));
			}
		);
	}
	#loglogistic
	if (("loglogistic" %in% parametric_fits)){
		#tryCatch(
		#	{

			if (T){
				library(eha) ## make "dllogis" and "pllogis" available to the working environment
		
				custom.llogis <- list(name="llogis",
					pars=c("shape","scale"),
					location="scale",
					transforms=c(base::log, base::log),
					inv.transforms=c(base::exp, base::exp),
					inits=function(t){ c(1, median(t)) })
				#browser()
				model_loglogistic <- flexsurvreg(survival_object~1,dist=custom.llogis);

				loglogistic_regression_model <- model_loglogistic;
				loglogistic_shape = model_loglogistic$res[1]
				loglogistic_scale = model_loglogistic$res[2]
				loglogistic_AIC = model_loglogistic[["AIC"]];

				loglogistic_shape_95_conf_l = model_loglogistic$res[1,2]
				loglogistic_shape_95_conf_h = model_loglogistic$res[1,3]
				loglogistic_scale_95_conf_l = model_loglogistic$res[2,2]
				loglogistic_scale_95_conf_h = model_loglogistic$res[2,3]
				successful_fits = c("loglogistic",successful_fits)


			}
			else{

				model_loglogistic = survreg(formula=survival_object~1,dist="loglogistic");
				loglogistic_regression_model <- model_loglogistic;
				loglogistic_scale = exp(coef(model_loglogistic));
				loglogistic_shape = 1/model_loglogistic$scale;

			}
			#},
			#error = function(e) {warning(paste("Problem fitting loglogistic:",e));}
		#);
	}

	if (("gompertz" %in% parametric_fits)){
		tryCatch(
			{
				model_gompertz <- flexsurvreg(survival_object~1,dist="gompertz")
				gompertz_regression_model <- model_gompertz;
				gompertz_shape <- model_gompertz$res[1]
				gompertz_rate <- model_gompertz$res[2]

				gompertz_shape_95_conf_l <- model_gompertz$res[1,2]
				gompertz_shape_95_conf_h <- model_gompertz$res[1,3]
				gompertz_rate_95_conf_l <- model_gompertz$res[2,2]	
				gompertz_rate_95_conf_h <- model_gompertz$res[2,3]

				gompertz_AIC = model_gompertz[["AIC"]]
				
				successful_fits = c("gompertz",successful_fits)
			},
			error = function(e){warning(paste("Problem fitting Gompertz fit: ",e$message));}
		);
	}

	fits = data.frame(
		lognormal_intercept=lognormal_intercept,
		lognormal_scale=lognormal_scale,
		lognormal_AIC=lognormal_AIC,

		weibull_scale=weibull_scale,
		weibull_shape=weibull_shape,
		weibull_AIC=weibull_AIC,
		weibull_shape_95_conf_l=weibull_shape_95_conf_l,
		weibull_shape_95_conf_h=weibull_shape_95_conf_h,
		weibull_scale_95_conf_l=weibull_scale_95_conf_l,
		weibull_scale_95_conf_h=weibull_scale_95_conf_h,

		inverse_gaussian_mu=inverse_gaussian_mu,
		inverse_gaussian_sigma2=inverse_gaussian_sigma2,
		inverse_gaussian_mu_95_conf_l=inverse_gaussian_mu_95_conf_l,
		inverse_gaussian_mu_95_conf_h=inverse_gaussian_mu_95_conf_h ,
		inverse_gaussian_sigma_95_conf_l=inverse_gaussian_sigma_95_conf_l,
		inverse_gaussian_sigma_95_conf_h=inverse_gaussian_sigma_95_conf_h,
		inverse_gaussian_AIC=inverse_gaussian_AIC,
		
		generalized_extreme_value_shape = generalized_extreme_value_shape,
		generalized_extreme_value_scale = generalized_extreme_value_scale,
		generalized_extreme_value_loc = generalized_extreme_value_loc ,
		generalized_extreme_value_AIC = generalized_extreme_value_AIC,

		weibull_gamma_frailty_scale=weibull_gamma_frailty_scale,
		weibull_gamma_frailty_shape=weibull_gamma_frailty_shape,
		weibull_gamma_frailty_theta=weibull_gamma_frailty_theta,
		weibull_gamma_frailty_AIC=weibull_gamma_frailty_AIC,
		weibull_gamma_frailty_shape_95_conf_l=weibull_gamma_frailty_shape_95_conf_l,
		weibull_gamma_frailty_shape_95_conf_h=weibull_gamma_frailty_shape_95_conf_h,
		weibull_gamma_frailty_scale_95_conf_l=weibull_gamma_frailty_scale_95_conf_l,
		weibull_gamma_frailty_scale_95_conf_h=weibull_gamma_frailty_scale_95_conf_h,
		weibull_gamma_frailty_theta_95_conf_l=weibull_gamma_frailty_theta_95_conf_l,
		weibull_gamma_frailty_theta_95_conf_h=weibull_gamma_frailty_theta_95_conf_h,

		gompertz_rate= gompertz_rate,
		gompertz_shape=gompertz_shape,
		gompertz_shape_95_conf_l=gompertz_shape_95_conf_l,
		gompertz_shape_95_conf_h=gompertz_shape_95_conf_h,
		gompertz_rate_95_conf_l=gompertz_rate_95_conf_l,
		gompertz_rate_95_conf_h=gompertz_rate_95_conf_h,
		gompertz_AIC=gompertz_AIC,

		gompertz_gamma_frailty_shape=gompertz_gamma_frailty_shape,
		gompertz_gamma_frailty_rate=gompertz_gamma_frailty_rate,
		gompertz_gamma_frailty_theta=gompertz_gamma_frailty_theta,
		gompertz_gamma_frailty_AIC=gompertz_gamma_frailty_AIC,
		gompertz_gamma_frailty_shape_95_conf_l=gompertz_gamma_frailty_shape_95_conf_l,
		gompertz_gamma_frailty_shape_95_conf_h=gompertz_gamma_frailty_shape_95_conf_h,
		gompertz_gamma_frailty_rate_95_conf_l=gompertz_gamma_frailty_rate_95_conf_l,
		gompertz_gamma_frailty_rate_95_conf_h=gompertz_gamma_frailty_rate_95_conf_h,
		gompertz_gamma_frailty_theta_95_conf_l=gompertz_gamma_frailty_theta_95_conf_l,
		gompertz_gamma_frailty_theta_95_conf_h=gompertz_gamma_frailty_theta_95_conf_h,

		gompertz_alt_gamma_frailty_initial=gompertz_alt_gamma_frailty_initial,
		gompertz_alt_gamma_frailty_scale=gompertz_alt_gamma_frailty_scale,
		gompertz_alt_gamma_frailty_theta=gompertz_alt_gamma_frailty_theta,
		gompertz_alt_gamma_frailty_AIC=gompertz_gamma_frailty_AIC,
		gompertz_alt_gamma_frailty_initial_95_conf_l=gompertz_alt_gamma_frailty_initial_95_conf_l,
		gompertz_alt_gamma_frailty_initial_95_conf_h=gompertz_alt_gamma_frailty_initial_95_conf_h,
		gompertz_alt_gamma_frailty_scale_95_conf_l=gompertz_alt_gamma_frailty_scale_95_conf_l,
		gompertz_alt_gamma_frailty_scale_95_conf_h=gompertz_alt_gamma_frailty_scale_95_conf_h,
		gompertz_alt_gamma_frailty_theta_95_conf_l=gompertz_alt_gamma_frailty_theta_95_conf_l,
		gompertz_alt_gamma_frailty_theta_95_conf_h=gompertz_alt_gamma_frailty_theta_95_conf_h,

		gompertz_makeham_shape=gompertz_makeham_shape,
		gompertz_makeham_rate=gompertz_makeham_rate,
		gompertz_makeham_constant=gompertz_makeham_constant,
		gompertz_makeham_AIC=gompertz_makeham_AIC,
		gompertz_makeham_shape_95_conf_l=gompertz_makeham_shape_95_conf_l,
		gompertz_makeham_shape_95_conf_h=gompertz_makeham_shape_95_conf_h,
		gompertz_makeham_rate_95_conf_l=gompertz_makeham_rate_95_conf_l,
		gompertz_makeham_rate_95_conf_h=gompertz_makeham_rate_95_conf_h,
		gompertz_makeham_constant_95_conf_l=gompertz_makeham_constant_95_conf_l,
		gompertz_makeham_constant_95_conf_h=gompertz_makeham_constant_95_conf_h,

		loglogistic_scale = loglogistic_scale, loglogistic_shape=loglogistic_shape,
		loglogistic_AIC=loglogistic_AIC,
		loglogistic_shape_95_conf_l=loglogistic_shape_95_conf_l,
		loglogistic_shape_95_conf_h=loglogistic_shape_95_conf_h,
		loglogistic_scale_95_conf_l=loglogistic_scale_95_conf_l,
		loglogistic_scale_95_conf_h=loglogistic_scale_95_conf_h,

		generalized_f_mu = generalized_f_mu, 
		generalized_f_sigma = generalized_f_sigma, 
		generalized_f_P = generalized_f_P, 
		generalized_f_Q = generalized_f_Q, 
		generalized_f_AIC = generalized_f_AIC,


		parametric_fit_lower_bound = parametric_fit_lower_bound,
		parametric_fit_upper_bound = parametric_fit_upper_bound
		);
		#browser();
	if (dim(fits)[1] > 1)
		browser()
	return(list(fits=fits,
		    successful_fits = successful_fits,
			models=list(
	    lognormal=lognormal_regression_model,
	    weibul=weibull_regression_model,
	    weibull_gamma_frailty=weibull_gamma_frailty_regression_model,
	    gompertz_gamma_frailty=gompertz_gamma_frailty_regression_model,
	    gompertz_alt_gamma_frailty=gompertz_alt_gamma_frailty_regression_model,
	    gompertz_makeham=gompertz_makeham_regression_model,
	    gompertz=gompertz_regression_model,
	    generalized_extreme_value=generalized_extreme_value_regression_model,
	    loglogistic=loglogistic_regression_model))); 
}
#Does not handle censoring!
ns_one_sided_ks_test=function(deaths,parameters,parameterization){
	if (parameterization == "lognormal"){
 		p = ks.test(deaths,"plnorm",parameters$lognormal_intercept,parameters$lognormal_scale);
 	}
 	else if (parameterization == "loglogistic"){
 		p = ks.test(deaths,"pllogis",scale=parameters$loglogistic_scale,shape=parameters$loglogistic_shape);
 	}
 	else if (parameterization == "cropped_loglogistic"){
 		p = ks.test(deaths,"pllogis",scale=parameters$cropped_loglogistic_scale,shape=parameters$cropped_loglogistic_shape);
 	}
 	else if (parameterization == "weibull"){
 		p = ks.test(deaths,"pweibull",scale=parameters$weibull_scale,shape=parameters$weibull_shape);				
 	}
 	else if (parameterization == "cropped_weibull"){
 		#print(t)
 		p = ks.test(deaths,"pweibull",scale=parameters$cropped_weibull_scale,shape=parameters$cropped_weibull_shape);				
 	}
 	else if (parameterization == "inverse_gaussian"){
 		p = ks.test(deaths,"ns_pinvgauss",mu=parameters$inverse_gaussian_mu,sigma2=parameters$inverse_gaussian_sigma2);
 	}
 	else if (parameterization == "gompertz"){
 		#browser();
 		fs_pgompertz = flexsurv::pgompertz
 		p = ks.test(deaths,"fs_pgompertz",shape=parameters$gompertz_shape,rate=parameters$gompertz_rate)
 	}else if (parameterization == "gompertz_gamma_frailty"){
 		p = ks.test(deaths,"pmgompertz",shape=parameters$gompertz_gamma_frailty_shape,rate=parameters$gompertz_gamma_frailty_rate,theta=parameters$gompertz_gamma_frailty_theta)
 	}
 	else if (parameterization == "cropped_gompertz"){
 		s = parameters$cropped_gompertz_shape;
 		r = parameters$cropped_gompertz_rate
 		fs_pgompertz = flexsurv::pgompertz
 		p = ks.test(deaths,"fs_pgompertz",shape=parameters$cropped_gompertz_shape,rate=parameters$cropped_gompertz_rate)
 
 	}
 	else if (parameterization == "generalized_f"){
 
 		p = ks.test(deaths,"pgenf",
 			mu=parameters$generalized_f_mu,
 			sigma=parameters$generalized_f_sigma,
 			Q=parameters$generalized_f_Q,
 			P=parameters$generalized_f_P)
 	}
 	else if (parameterization == "weibull_gamma_frailty"){
 		p = ks.test(deaths,"pmllogis",scale=parameters$weibull_gamma_frailty_scale,shape=parameters$weibull_gamma_frailty_shape,theta=parameters$weibull_gamma_frailty_theta);
 	}
 	else if (parameterization == "gompertz_gamma_frailty"){
 		p = ks.test(deaths,"pmgompertz",shape=parameters$gompertz_gamma_frailty_shape,rate=parameters$gompertz_gamma_frailty_rate,theta=parameters$gompertz_gamma_frailty_theta);
 	}
 	else if (parameterization == "gompertz_makeham"){
 		#browser()
 		eha_pmakeham = eha::pmakeham
 		p = ks.test(deaths,"eha_pmakeham",shape=c(parameters$gompertz_makeham_shape,parameters$gompertz_makeham_constant),1/parameters$gompertz_makeham_rate);
 	}
	else stop(paste("Unknown parameterization: ",parameterization));
	return (p);
}

ns_get_parametric_survival = function(tt,parameters,parameterization,survival_by="survival"){
	if (survival_by=="quantile"){
		suffix="q";
	}else if (survival_by=="survival"){
		suffix="p";
	}else if (survival_by=="random_sample"){
		suffix="r";
	}else if (survival_by=="pdf"){
		suffix="d";
	}else if (survival_by=="hazard"){
		return( ns_get_parametric_survival(tt,parameters,parameterization,"pdf")/ns_get_parametric_survival(tt,parameters,parameterization,"survival"))
	}
	else stop(paste("Unknown survival by specification = ",survival_by));
	
	pfs_gompertz = flexsurv::pgompertz
	qfs_gompertz = flexsurv::qgompertz
	dfs_gompertz = flexsurv::dgompertz
		
	#	browser()
	if (parameterization == "lognormal"){
		p = do.call(paste0(suffix,"lnorm"),args=list(tt,parameters$lognormal_intercept,parameters$lognormal_scale));
	}
	else if (parameterization == "loglogistic"){
		p = do.call(paste0(suffix,"llogis"),args=list(tt,scale=parameters$loglogistic_scale,shape=parameters$loglogistic_shape));
	}
	else if (parameterization == "cropped_loglogistic"){
		p = do.call(paste0(suffix,"llogis"),args=list(tt,scale=parameters$cropped_loglogistic_scale,shape=parameters$cropped_loglogistic_shape));
	}
	else if (parameterization == "weibull"){
		p = do.call(paste0(suffix,"weibull"),args=list(tt,scale=parameters$weibull_scale,shape=parameters$weibull_shape));				
	}
	else if (parameterization == "cropped_weibull"){
		#print(t)
		p = do.call(paste0(suffix,"weibull"),args=list(tt,scale=parameters$cropped_weibull_scale,shape=parameters$cropped_weibull_shape));				
	}
	else if (parameterization == "inverse_gaussian"){
		#rowser()
		p = do.call(paste0(suffix,"ns_invgauss"),args=list(tt,mu=parameters$inverse_gaussian_mu,sigma2=parameters$inverse_gaussian_sigma2));
	}
	else if (parameterization == "gompertz"){
		
		p = do.call(paste0(suffix,"fs_gompertz"),args=list(tt,shape=parameters$gompertz_shape,rate=parameters$gompertz_rate))
	}else if (parameterization == "gompertz_gamma_frailty"){
		p = do.call(paste0(suffix,"mgompertz"),args=list(tt,shape=parameters$gompertz_gamma_frailty_shape,rate=parameters$gompertz_gamma_frailty_rate,theta=parameters$gompertz_gamma_frailty_theta))
	}
	else if (parameterization == "cropped_gompertz"){
		s = parameters$cropped_gompertz_shape;
		r = parameters$cropped_gompertz_rate
		p = do.call(paste0(suffix,"fs_gompertz"),args=list(tt,shape=parameters$cropped_gompertz_shape,rate=parameters$cropped_gompertz_rate))

	}
	else if (parameterization == "generalized_f"){

		p = do.call(paste0(suffix,"genf"),args=list(tt,
			mu=parameters$generalized_f_mu,
			sigma=parameters$generalized_f_sigma,
			Q=parameters$generalized_f_Q,
			P=parameters$generalized_f_P))
	}
	else if (parameterization == "weibull_gamma_frailty"){
		p = do.call(paste0(suffix,"mllogis"),args=list(tt,scale=parameters$weibull_gamma_frailty_scale,shape=parameters$weibull_gamma_frailty_shape,theta=parameters$weibull_gamma_frailty_theta));
	}
	else if (parameterization == "gompertz_gamma_frailty"){
		p = do.call(paste0(suffix,"mgompertz"),args=list(tt,shape=parameters$gompertz_gamma_frailty_shape,rate=parameters$gompertz_gamma_frailty_rate,theta=parameters$gompertz_gamma_frailty_theta));
	}
	else if (parameterization == "gompertz_makeham"){
		#browser()
		p = do.call(paste0("eha::",suffix,"makeham"),args=list(tt,shape=c(parameters$gompertz_makeham_shape,parameters$gompertz_makeham_constant),1/parameters$gompertz_makeham_rate));
	}
	else stop(paste("Unknown parameterization: ",parameterization));
	if (survival_by=="survival")
		return(1-p);
	return (p);
}
