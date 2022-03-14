pmgompertz2 <- function (q, initial = 1, scale = 1,theta=1,log.p=F) {
	theta[theta<.0000001] = .0000001
	return(pmgompertz(q,shape=initial/scale,rate=scale,theta=theta,log.p=log.p))
}

qmgompertz2<- function (q, initial = 1, scale = 1,theta=1,log.p=F){
	
	theta[theta<.0000001] = .0000001
	return(qmgompertz(q,shape=initial/scale,rate=scale,theta=theta,log.p=log.p))
}


dmgompertz2 <- function (x, initial = 1, scale = 1,theta=1,log=F){

	theta[theta<.0000001] = .0000001
	return(dmgompertz(x,shape=initial/scale,rate=scale,theta=theta,log=log))
}


scale = 1;
initial = .000000000001;
theta = 0.001;
if (0){
	x = seq(from=0,to=10,by=.1)
	tt = dmgompertz2(x,initial,scale,theta,log=T)/pmgompertz2(x,initial,scale,theta,log=T)
	tt1 = flexsurv::dgompertz(x,initial/scale,scale,log=T)/flexsurv::pgompertz(x,initial/scale,scale,log=T)
	xlim=c(.01,8);
	indx = x >= xlim[1] & x <= xlim[2];
	vals = c(tt[indx],tt1[indx])
	ylim = c(min(vals),max(vals))
	plot(x,tt,type="l",col="red",ylim=ylim,xlim=xlim,log="y");
	lines(x,tt1,lty=2);
	
	x = seq(from=0,to=10,by=.1)
	tt = dmgompertz2(x,initial,scale*4,theta,log=T)/pmgompertz2(x,initial,scale,theta,log=T)
	
	lines(x,tt,type="l",col="red",ylim=ylim,xlim=xlim,log="y");
	
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

custom.mgompertz2 <- list(name="mgompertz2",
	pars=c("initial","scale","theta"),
	location="scale",
	transforms=c(base::identity, base::log,base::log),
	inv.transforms=c(base::identity, base::exp,base::exp),
	inits=function(t){ c(.00001, mean(t),.01) })
