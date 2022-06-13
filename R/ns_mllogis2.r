#Implements the closed form distribution of a weibull distribution with gamma frailty

 min_mllogis_theta <<- .00001;


#' @export
pmllogis <- function (q, shape = 1, scale = 1,theta=1,log.p=F) 
{

    if (any(c(shape, scale) <= 0) || theta < 0) {
    #	stop("INVALID VALUE");
        warning(paste("Non-positive shape, scale: ",shape, scale, theta))
        return(rep(NaN,length(q)))
    }
  
   
    if (!log.p){
	ret <- ifelse(q <= 0, 0, 

		ifelse(rep(theta,length(q)) > min_mllogis_theta,
			1-(1+theta*(q/scale)^(shape))^(-1/theta),
			pweibull(q,shape=shape,scale=scale,log.p=log.p)
			)

		)
   
    }
    else{
	ret <- ifelse(q <= 0, -Inf, 
			ifelse(theta > min_mllogis_theta,
				log(1-(1+theta*(q/scale)^(shape))^(-1/theta)),
				pweibull(q,shape=shape,scale=scale,log.p=log.p)
			)
		)
    }
  #  if (unique(ret) == 0)
  #  	browser();
    
 #   print(unique(ret))
    return(ret)
}

#' @export
qmllogis <- function (q, shape = 1, scale = 1,theta=1,log.p=F) {
#stop("WJA")
    if (any(c(shape, scale) <= 0) || theta < 0) {
        warning(paste("Non-positive shape, scale, or theta: ",shape, scale, theta))
       return(rep(NaN,length(q)))
    }
   
if (log.p == F){
    ret <- ifelse(q <= 0, 0, 
		ifelse(rep(theta,length(q)) >= min_mllogis_theta,
			scale*(((1-q)^(-theta) - 1)/theta)^(1/(shape)),
			qweibull(q,shape=shape,scale=scale,log.p=log.p)
		)
	   )
   }
   else{
    stop("log.p")
}
	
    return(ret)
}


#' @export
dmllogis <- function (x, shape = 1, scale = 1,theta=1,log=F) 
{
	
    if (any(c(shape, scale) <= 0) || theta < 0) {
    #	browser();#browser()
    #	stop("BAD PARAM")
        warning(paste("Non-positive shape, scale, or theta: ",shape, scale, theta))
        return(rep(NaN,length(x)))
    }

if (log == F){
    s = x/scale;
    ret <- ifelse(x <= 0, 0, 
		ifelse(rep(theta,length(x)) >= min_mllogis_theta,
			((shape/scale)*s^(shape-1))*(1+theta*s^shape)^(-(theta+1)/theta) ,
			dweibull(x,shape=shape,scale=scale,log=log)
		)
	)
   }
   else{
    s = x/scale;
    ret <- ifelse(x <= 0, 0, 
	ifelse(rep(theta,length(x)) >= min_mllogis_theta,
		log(shape/scale)+(shape-1)*log(s)+(-(theta+1)/theta)*log(1+theta*s^shape),
		dweibull(x,shape=shape,scale=scale,log=log)
		)
	)
   }
    #if (any(is.infinite(ret)) || any(is.na(ret)) || any(is.na(ret)) < 0)
    #	browser();
    #print(unique(ret))
    return(ret)
}

custom.mllogis <- list(name="mllogis",
pars=c("shape","scale","theta"),
location="scale",
transforms=c(base::log, base::log,base::log),
inv.transforms=c(base::exp, base::exp,base::exp),
inits=function(t){ c(12, stats::median(t),1) })


if (0){
t = seq(0,.01,.00001)
tt = t;
for (i in 1:length(t))
	tt[i] = dmllogis(.5,9,1,t[i])
plot(t,tt,type='l');
}

rmllogis  = function(n, shape = 1, scale = 1,theta=1){
	qmllogis(runif(n),shape,scale,theta);
}
