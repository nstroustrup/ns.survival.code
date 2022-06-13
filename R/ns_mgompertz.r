
#' @export
pmgompertz <- function (q, shape = 1, rate = 1,theta=1,log.p=F) 
{
    if (any(c(theta,rate) <= 0)) {
        warning(paste("Non-positive shape, rate, or theta: ",shape, rate, theta))
        return(NaN)
    }
    if (shape == 0){
    	if (!log.p){
    		return (ifelse(q <= 0, 0, 1-(1+theta*rate*q)^(-1/theta)));
    	}
    	else return (ifelse(q <= 0, -Inf, log(1-(1+theta*rate*q)^(-1/theta))));
    }
    if (!log.p){
    	ret <- ifelse(q <= 0, 0, 1-(1-theta*rate/shape*(1-exp(shape*q)))^(-1/theta))
    }
    else{
    	ret <- ifelse(q <= 0, -Inf, log(1-(1-theta*rate/shape*(1-exp(shape*q)))^(-1/theta)))
    }
    return(ret)
}

#' @export
qmgompertz<- function (q, shape = 1, rate = 1,theta=1,log.p=F) 
{	
    if (any(c(rate,theta) <= 0)) {
        warning(paste("Non-positive shape, rate, or theta: ",shape, rate, theta))
        return(NaN)
    }
    if (shape == 0){
    	if (!log.p){
    		ret <- ifelse(q <= 0, 0, ((1-q)^(-theta)-1)/(rate*theta) );
    	}
    	else{
    	stop("log.p");
    	}
    	return(ret)
    }
   if (log.p == F){
    ret <- ifelse(q <= 0, 0, (1/shape)*log( shape/(theta*rate)*((1-q)^(-theta)-1)+1) );
   }
   else{
    stop("log.p")
   }
	
    return(ret)
}


#' @export
dmgompertz <- function (x, shape = 1, rate = 1,theta=1,log=F) 
{
    if (any(c( rate,theta) <= 0)) {
        warning(paste("Non-positive shape, rate, or theta: ",shape, rate, theta))
        return(NaN)
    }
    if (shape==0){
    #	print("WHA")
    	ret <- ifelse(x < 0, -Inf, log(rate)+(-(theta+1)/theta)*log(1+rate*theta*x))
    }
    else{
   	 ret <- ifelse(x < 0, -Inf, log(rate)+shape*x+(-(theta+1)/theta)*log(rate*theta/shape*(exp(shape*x)-1)+1))
   }
    if(log==F)
    	return(exp(ret))
    return(ret)
}

if (0){
a=.1;
b = 1;
theta = .000000000001;
#theta = .1;
	x = seq(from=0,to=10,by=.1)
	tt = dmgompertz(x,b,a,theta,log=T)
	tt1 = flexsurv::dgompertz(x,b,a,log=T)
	xlim=c(0,4);
	indx = x >= xlim[1] & x <= xlim[2];
	vals = c(tt[indx],tt1[indx])
	ylim = c(min(vals),max(vals))
	plot(x,tt,type="l",col="red",ylim=ylim,xlim=xlim);
	lines(x,tt1,lty=2);
}else if (0){
	q = seq(from=.99,to=.01,by=-.05);
	tt = qmgompertz(q,b,a,theta);
	tt1 = qgompertz(q,b,a);
	ylim = c(0,1);
	xlim = c(min(c(tt,tt1)),max(c(tt,tt1)));
	plot(tt,q,type="l",col="red",ylim=ylim,xlim=xlim);
	lines(tt1,q,lty=2)
	
}else if (0){
	x = seq(from=0,to=10,by=.1)
	tt = pmgompertz(x,b,a,theta,log=F)
	tt1 = flexsurv::pgompertz(x,b,a,log=F)
	xlim=c(0,4);
	indx = x >= xlim[1] & x <= xlim[2];
	vals = c(tt[indx],tt1[indx])
	ylim = c(min(vals),max(vals))
	plot(x,tt,type="l",col="red",ylim=ylim,xlim=xlim);
	lines(x,tt1,lty=2);
}

#' @export
pns_makeham <- function(q, shape,rate,makeham,lower.tail = TRUE, log.p = FALSE){
			if (any(c(shape,1/rate,makeham) <= 0)) {
				warning(paste("Non-positive shape, rate: ",shape, rate))
				return(NaN)
			}		
			return(eha::pmakeham(q,shape=c(shape,makeham),scale=1/rate,lower.tail,log.p));
		}
#' @export
qns_makeham <- function(p, shape,rate,makeham,lower.tail = TRUE, log.p = FALSE){
			if (any(c(shape,1/rate,makeham) <= 0)) {
					warning(paste("Non-positive shape, rate: ",shape, rate))
					return(NaN)
			}
			return(eha::qmakeham(p,shape=c(shape,makeham),scale=1/rate,lower.tail,log.p));
		}
#' @export
dns_makeham <- 	function(x, shape,rate,makeham,log = FALSE,debug=T){
			if(debug)
			params <<- list(x=x,shape=shape,rate=rate,makeham=makeham)
			print(c(params[["shape"]],params[["rate"]],params[["makeham"]]))
			if (any(c(shape,1/rate,makeham) <= 0)) {
				warning(paste("Non-positive shape, rate: ",shape, rate))
				return(NaN)
			}
			return(eha::dmakeham(x,shape=c(shape,makeham),scale=1/rate,log));
		}

#' @export
custom.mgompertz <- list(name="mgompertz",
	pars=c("shape","rate","theta"),
	location="rate",
	transforms=c(base::identity, base::log,base::log),
	inv.transforms=c(base::identity, base::exp,base::exp),
	inits=function(t){ c(0,1 / mean(t),1) })

#' @export
custom.makeham <- list(name="ns_makeham",
	pars=c("shape","rate","makeham"),
	location="rate",
	transforms=c(base::log, base::log,base::log),
	inv.transforms=c(base::exp, base::exp,base::exp),
	inits=function(t){ c(0.002,1/mean(t),.0001) }
	)