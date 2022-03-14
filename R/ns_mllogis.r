#Implements the closed form distribution of a weibull distribution with gamma frailty

 min_mllogis_theta <<- .00001;


pmllogis <- function (q, shape = 1, scale = 1,theta=1,log.p=F) 
{
    if (any(c(shape, scale) <= 0) || theta < 0) {
        warning(paste("Non-positive shape, scale: ",shape, scale, theta))
        return(NaN)
    }
  
   
    if (!log.p){
	ret <- ifelse(q <= 0, 0, 1-(1+theta*(q/scale)^(shape))^(-1/theta))
    }
    else{
	ret <- ifelse(q <= 0, -Inf, log(1-(1+theta*(q/scale)^(shape))^(-1/theta)))
    }
   indx = theta < min_mllogis_theta
   ret[indx] <- pweibull(q[indx],shape=shape[indx],scale=scale[indx],log.p=log.p);
    
    
    return(ret)
}

qmllogis <- function (q, shape = 1, scale = 1,theta=1,log.p=F) 
{	

    if (any(c(shape, scale) <= 0) || theta < 0) {
        warning(paste("Non-positive shape, scale, or theta: ",shape, scale, theta))
        return(NaN)
    }
     
	if (log.p == F){
	    ret <- ifelse(q <= 0, 0, scale*(((1-q)^(-theta) - 1)/theta)^(1/(shape)))
	   }
	   else{
	    stop("log.p")
	 #   ret <- ifelse(q <= 0, -Inf, scale*(((1-exp(q))^(-theta) - 1)/theta)^(1/(shape)))
	   }
   indx = theta < min_mllogis_theta
   ret[indx] = qweibull(q[indx],shape=shape[indx],scale=scale[indx],log.p=log.p);
   
	
    return(ret)
}


dmllogis <- function (x, shape = 1, scale = 1,theta=1,log=F) 
{
    if (any(c(shape, scale) <= 0) || theta < 0) {
        warning(paste("Non-positive shape, scale, or theta: ",shape, scale, theta))
        return(NaN)
    }

if (log == F){
    s = x/scale;
    ret <- ifelse(x <= 0, 0, ((shape/scale)*s^(shape-1))*(1+theta*s^shape)^(-(theta+1)/theta) )
   }
   else{
    s = x/scale;
    ret <- ifelse(x <= 0, 0, log(shape/scale)+(shape-1)*log(s)+(-(theta+1)/theta)*log(1+theta*s^shape))
   }
 indx = theta < min_mllogis_theta
 ret[indx] = dweibull(x[indx],shape=shape[indx],scale=scale[indx],log=log);
    
	
    return(ret)
}

custom.mllogis <- list(name="mllogis",
	pars=c("shape","scale","theta"),
	location="scale",
	transforms=c(base::log, base::log,base::log),
	inv.transforms=c(base::exp, base::exp,base::exp),
	inits=function(t){ c(1, stats::median(t),1) }
)


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

