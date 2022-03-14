ns_area_between_km = cmpfun(function(t_a,t_b,s_a,s_b,distance_norm=1){
	#calculates the area between two kaplan meier survival curves
	#if the last event is censored for one curve,
	#the difference between the two curves is only calculated up to that point.
	
	s_a_1 = s_a==1;
	s_b_1 = s_b==1;
	
	s_a = s_a[!s_a_1];
	t_a = t_a[!s_a_1];
	s_b = s_b[!s_b_1];
	t_b = t_b[!s_b_1];
	
	last_a_censored = length(which(s_a==0)) == 0
	last_b_censored = length(which(s_b==0)) == 0
	
	#browser()
	d = rbind(data.frame(t=t_a,s_a=s_a,s_b=NA),
		  data.frame(t=t_b,s_a=NA, s_b=s_b) );
	d = d[order(d$t),]
	#print(d)
	
	#if curve a has the first death
	#then curve b needs to have its survival set to 1 at that time.
	#if curve be has the first death
	#ten curve a nees to have its survival set to 1 at that time.
	last_s_a = NA
	last_s_b = NA
	min_t_a_i = which.min(t_a)
	min_t_b_i = which.min(t_b)
	min_t_a = t_a[min_t_a_i];
	min_t_b = t_b[min_t_b_i];
	
	if (min_t_a == min_t_b){
		last_s_a = last_s_b = 1;
	}else if (min_t_a > min_t_b){
		last_s_a = 1
		last_s_b = s_a[min_t_b_i]
	}else{
		last_s_a = s_a[min_t_a_i]
		last_s_b = 1;
	}
	
	sum = 0;
	
	last_t = dim(d)[1];
	#print(paste("Last t:",last_t))
	if (last_a_censored)
		last_t = max(which(!is.na(d$s_a)))
	if (last_b_censored){
		last_t = min( max(which(!is.na(d$s_b))),
			      last_t);
	}
	#print(paste("Last t:",last_t))
	
	#duplicate t's will be handled correctly.
	#the first of a duplicate t encountered will produce the correct result
	#and subsequent duplicate ts will have a dt = 0, thus producing no effect the final sum
	#but updating last_s_a and last_s_b correctly.
	
	if (distance_norm == 1){
		for (i in 2:last_t){

			s_a_i = last_s_a;
			s_b_i = last_s_b;

			if (!is.na(d$s_a[i]))
				last_s_a = d$s_a[i];

			if (!is.na(d$s_b[i]))
				last_s_b = d$s_b[i];

			sum = sum + abs(s_a_i-s_b_i)*(d$t[i]-d$t[i-1]);
		}
		
	}
	else if (distance_norm == 2){
		for (i in 2:last_t){

			s_a_i = last_s_a;
			s_b_i = last_s_b;

			if (!is.na(d$s_a[i]))
				last_s_a = d$s_a[i];

			if (!is.na(d$s_b[i]))
				last_s_b = d$s_b[i];

			sum = sum + (s_a_i-s_b_i)^2*(d$t[i]-d$t[i-1]);
		}
		sum = sqrt(sum);
	}
	else if (distance_norm == 3){
		for (i in 2:last_t){
			s_a_i = last_s_a;
			s_b_i = last_s_b;

			if (!is.na(d$s_a[i]))
				last_s_a = d$s_a[i];

			if (!is.na(d$s_b[i]))
				last_s_b = d$s_b[i];

			sum = max(sum,abs(s_a_i-s_b_i)*(d$t[i]-d$t[i-1]));
		}
	}
	return(sum);
})


ns_calculate_a_b_n =cmpfun(
function (s){

	#total number of deaths or censored events
	N0 = length(s[,1])
	ss = s;
	ordered_events = s[order(s[,1]),]
	t = ordered_events[,1]
	c = ordered_events[,2]
	
	#number of animals that left observation prior to time index
	NN = 0:(N0-1);
	#number of animals in observation prior to time index
	N = N0:1;
	#parameters to be calculated
	res = matrix(nrow=N0,ncol=5);
	colnames(res) = c("t","a","b","d","N");
	#res = data.frame(t=dd,
	#		 a=dd,
	#		 b=dd,
	#		 d=dd,
	#		 N=dd,
	#		 );
	#iterate to calculate alpha and beta
	prev_b = 0;
	prev_a = 0
	j = 1;
	total_skipped=0;
	#browser();
	while(1){
		res[j-total_skipped,"t"] = t[j];
		cur_z = t == t[j]
		cur_e = which(cur_z);
		
		res[j-total_skipped,"N"] = N[min(cur_e)]
		
		total_events_skipped_this_round = length(cur_e)-1;
		
		a_N = N0-N[min(cur_e)]-1
		
		if (a_N<0){
			res[j-total_skipped,"a"] = 0;
		}else{
			res[j-total_skipped,"a"] = sum(1/(N0-0:a_N))-prev_b
		}
		#all events which occur at the current time
		cur_d = which(cur_z & c==1);
		if (length(cur_d) == 0){
			res[j-total_skipped,"b"] = prev_b;
			res[j-total_skipped,"d"] = 0; #no death events this round
			
			j = max(cur_e)+1;
			
			total_skipped=total_events_skipped_this_round+total_skipped;
			if (j > N0)
				break;
			next();
		}
			
		#do the sum in vetorized form by allocating c(1,2,3...)
		#and subtracting it from the number under observation at
		#time j (which is N0-NN[cur_d])
		
		res[j-total_skipped,"b"] = prev_b + sum(1/(N0-NN[min(cur_d)]-0:(length(cur_d)-1)));
		res[j-total_skipped,"d"] = length(cur_d);
		#go to the next new time
		prev_b = res[j-total_skipped,"b"];
		
		total_skipped=total_events_skipped_this_round+total_skipped;
		
		j = max(cur_e)+1;
		
		if (j > N0)
			break;
	}
		#if (s[1,1] == s[2,1])
		#browser()
	if (total_skipped == 0){
		return (res);
	}
	#indexing silliness to avoid implicit conversion of one-row matricies to numeric objects
	rr = c(rep(T,N0-total_skipped),rep(F,total_skipped))
	#browser()
	return(subset(res,rr));
},options=list(optimize=3) )


ns_find_after_i = function(needle,haystack,haystack_length,i){
	#browser()
	if (i+1 > haystack_length)	
		return(integer(0))
	for (j in (i+1):haystack_length){
		if (haystack[j]==needle)
			return(j)
		if (haystack[j]>needle)
			return(integer(0))
	}
	return(integer(0))
}


ns_min_larger_after_i = function(bound,haystack,haystack_length,i){
	if (i+1 > haystack_length)	
		return(integer(0))
	for (j in (i+1):haystack_length){
		if (haystack[j]>bound)
			return(j)
	}
	return(integer(0))
}

ns_min_eq_or_larger_after_i = function(bound,haystack,haystack_length,i){
	if (i+1 > haystack_length)	
		return(integer(0))
	for (j in (i+1):haystack_length){
		if (haystack[j]>=bound)
			return(j)
	}
	return(integer(0))
}

ns_merge_abn_lists = function(abn1,abn2){
	abn1_length = dim(abn1)[1];
	abn2_length = dim(abn2)[1];
	t = unique(c(abn1[,"t"],abn2[,"t"]));
	t = t[order(t)];
			 
	res = matrix(nrow=length(t),ncol=12);
	colnames(res) = c('t','N1','d1','b1','a1','N2','d2','b2','a2','n','u','Y');
	
	
	#The maximum i where abn1$t < t[i]
	i1=0
	i2=0
	p_abN1 = list(a=0,b=0,N=abn1[1,"N"],d=0)
	p_abN2 = list(a=0,b=0,N=abn2[1,"N"],d=0)
	N1_alive = abn1[1,"N"];
	N2_alive = abn2[1,"N"];
	out_j = 1;
	#browser();
	for (j in 1:length(t)){
		#print(c(i1,i2))
		current_uncensored_event=F
		valid_t1s = abn1[,"t"]>=i1;
		valid_t2s = abn2[,"t"]>=i2;
		
		#if s1 has an event at the current time
		#update b,d,and N
		#di_1 = which(abn1$t == t[j])
		di_1 = ns_find_after_i(t[j],abn1[,"t"],abn1_length,i1)
		if (length(di_1) > 0){
			p_abN1[['d']] = abn1[di_1,'d'];
			p_abN1[['b']] = abn1[di_1,'b'];
			p_abN1[['N']] = abn1[di_1,'N'];
			next_i1 = di_1;
		}
		else{
			#if the other curve has an event at the current time
			#get the N immediately after the previous event in s1
			#di_1 = which(abn1$t > t[j] )
			di_1 = ns_min_larger_after_i(t[j],abn1[,"t"],abn1_length,i1);
			if (length(di_1) > 0){
				p_abN1[['N']] = abn1[di_1,"N"]
			}else p_abN1[['N']] = 0;
			p_abN1[['d']] = 0;
			next_i1=i1;
		}
		#repeat as before for s2
		#di_2 = which(abn2[,"t"] == t[j])
		di_2 = ns_find_after_i(t[j],abn2[,"t"],abn2_length,i2)
		if (length(di_2) > 0){
			p_abN2[['d']] = abn2[di_2,"d"];
			p_abN2[['b']] = abn2[di_2,"b"];
			p_abN2[['N']] = abn2[di_2,"N"];
			next_i2 = di_2;
		}
		else{
			#di_2 = which(abn2$t > t[j])
			di_2 = ns_min_larger_after_i(t[j],abn2[,"t"],abn2_length,i2);
			if (length(di_2) > 0){
				p_abN2[['N']] = abn2[di_2,"N"];
			}else{
				p_abN2[['N']] = 0;
			}
			#browser();
			p_abN2[['d']] = 0;
			next_i2 = i2
		}
		#get the a that holds up until the next event in s1
		#ai_1 = which(abn1[,"t"] >= t[j] & valid_t1s)
		ai_1 = ns_min_eq_or_larger_after_i(t[j],abn1[,"t"],abn1_length,i1);
		if (length(ai_1) > 0)
			p_abN1[['a']] = abn1[ai_1,"a"];
			
		#repeat as before for s2
		#ai_2 = which(abn2[,"t"] >= t[j] & valid_t2s)
		ai_2 = ns_min_eq_or_larger_after_i(t[j],abn2[,"t"],abn2_length,i2);
		if (length(ai_2) > 0)
			p_abN2[['a']] = abn2[ai_2,"a"];
		
	
		if (p_abN1[['d']] == 0 && p_abN2[['d']] == 0)
			next;
		res[out_j,'t'] = t[j];
		res[out_j,'a1'] = p_abN1[['a']]
		res[out_j,'b1'] = p_abN1[['b']]
		
		res[out_j,'N1'] = p_abN1[['N']]
		res[out_j,'d1'] = p_abN1[['d']]
		
		res[out_j,'a2'] = p_abN2[['a']]
		res[out_j,'b2'] = p_abN2[['b']]
		res[out_j,'N2'] = p_abN2[['N']]
		res[out_j,'d2'] = p_abN2[['d']]
		i1 = next_i1
		i2 = next_i2
		out_j = out_j+1;	
	}
	res = res[1:(out_j-1),]
	return(res);
}

ns_calc_u = function(abn_merged){
	nT = dim(abn_merged)[1];
	#browser()
	prev_u = 0;
	j = 1;
	u = rep(NA,nT)
	for (j in 1:nT){
		
		#do the sum in vectorized form by allocating c(1,2,3...)
		#and subtracting it from the number under observation at
		#time j 
		l1=l2=0;
		if (abn_merged[j,'d1']>0)
			l1 = sum(1/(abn_merged[j,'N1']- 0:(abn_merged[j,'d1']-1)));
		if (abn_merged[j,'d2']>0)
			l2 = sum(1/(abn_merged[j,'N2']- 0:(abn_merged[j,'d2']-1)));
		u[j] = prev_u + abn_merged[j,'n']*( l1-l2 ) ;
		
		#go to the next new time
		prev_u = u[j]
	}
	
	#browser()
	return (u);
}
ns_P = function(V,R){

	#if all the animals are dead at the final time point, then R is 1.
	#the calculation becomes
	#1-pnorm(V/0) + pnorm(V/0)*exp(-2*V*V)
	#which is simply
	#exp(-2*V*V)
	
	return(
		1-pnorm(V/sqrt(R-R*R))+
		  pnorm(V*(2*R-1)/sqrt(R-R*R))*
			exp(-2*V*V)
	 )

}

ns_ks_test_ = function(s1,s2,one_sided=F,return_timeseries=F){
	N1 = length(s1[,1])
	N2 = length(s2[,1])
	##if (s2[1,1]==s2[2,1])
	#	browser()
	#calculate alphas and betas for each curve
	abn1 = ns_calculate_a_b_n(s1);
	abn2 = ns_calculate_a_b_n(s2);
	#print(dim(s1))
		#browser()
	abn = ns_merge_abn_lists(abn1,abn2)
	abn[,'n'] = 1/sqrt(1/(N1*exp(-1*abn[,'a1'])) + 1/(N2*exp(-1*abn[,'a2'])))
	
	abn[,'u'] = ns_calc_u(abn);
	abn[,'Y'] = .5*abn[,'u']*(exp(-1*abn[,'b1'])+exp(-1*abn[,'b2']))
	V = max(c(0,abn[,'Y']));
	
	#this is the standard KS statistic
	A = max(abs(abn[,'Y']));
	
	N = dim(abn)[1]
	R = 1-.5*(exp(-1*abn[N,'b1'])+exp(-1*abn[N,'b2']))
	if (one_sided){
		p = ns_P(V,R);
		test_statistic = V;
	}else{
		p = 2*ns_P(A,R);
		test_statistic = A;
		if (p > .2){
			warning("The approximation for the two-sided KS test is good only for p < .2.  Returning .2");
			p = .2;
		}
	}
	if (!return_timeseries)
		return (list(p_value=p,Y=test_statistic,R=R,N1=N1,N2=N2));
	return (list(p_value=p,Y=test_statistic,R=R,N1=N1,N2=N2,
			ts=data.frame(t=abn[,"t"],
					a=abn[,"a1"],
					b=abn[,"b1"],
					a=abn[,"a2"],
					b=abn[,"b2"],
					u=abn[,"u"],
					Y=abn[,'Y'])
					));
}

ns_ks_test = cmpfun(ns_ks_test_,options=list(optimize=3) )

ns_ks_test_group = function (formula,reference_group,one_sided=F,return_timeseries=F){
	f = terms(formula);
	r = eval.parent(attr(f,"variables"))
	p = attr(f,"factors")
	#browser()
	group = r[[which(p==1)]]
	deaths = r[[which(p==0)]]
	if (class(deaths) != "Surv")
		stop("Must provide a Surv object");
	ref = subset(deaths,group==reference_group);
	if (dim(ref)[1]==0)
		stop("Could not find any members of the reference group");
	to_test = unique(group);
	to_test = to_test[to_test != reference_group];
	stats = NULL;
	ts = NULL;
	for (i in 1:length(to_test)){
		t = subset(deaths,group==to_test[[i]]);
		res = ns_ks_test(ref,t,one_sided=one_sided,return_timeseries=return_timeseries);
		stats = rbind(data.frame(group=to_test[[i]],
					p_value = res$p_value,
					Y = res$Y,
					R = res$R,
					N1 = res$N1,
					N2 = res$N2),
					stats);
		if (return_timeseries){			
			res$ts$group = to_test[[i]];
			ts = rbind(res$ts,ts);
		}
	}
	if (return_timeseries){
		return (list(stats = stats,ts=ts));
	}else return (stats);
}



#call this function to test the ks code.
ns_test_ks_code = function(){
	s1 = Surv(c(28,89,175,195,309,377,393,421,447,462,709,744,770,1106, 1206),
		  c(1 , 1,  1,  1,  1,  0,  0,  0,  0,  1,  0,  0,  0,   0,    0));
	s2 = Surv(c(34, 88, 137, 199, 280,291,299,300,309,351,358,369,369,370,375,382,392,429,451,1119),
		  c( 1,  1,   1,   1,   1,  1,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,   0));

	res = ns_ks_test(s2,s1);
	if (abs(res$p_value - 0.002226671) > .000001){
		die("Test Failed");
	}else print("Test Passed.")
}
#ns_test_ks_code()

ns_test_area_code = function(){

	s1 =  c(1,1,1,1, 1,.5,.5,0)
	s2 =  c(1,1,1,1,.5,.5,.5,0)
	s1t = c(1,2,3,4,5,6,7,8,9)
	s2t = s1t
	res = ns_area_between_km(s1t,s2t,s1,s2);
	if (res != .5){
		print (paste("Error: ",res,"!= .5"));
	}else print("Test 1: OK");
			
	s1 =  c(1,1,1,1, 1,.5,.5,0)
	s2 =  c(1,1,1,.5,.5,.5,.5,0)
	s1t = c(1,2,3,4,5,6,7,8,9)
	s2t = s1t;
	res = ns_area_between_km(s1t,s2t,s1,s2);
	if (res != 1){
		print (paste("Error: ",res,"!= 1"));
	}else print("Test 2: OK");
	
	s1 =  c(1,1,1,1, 1,.5,.5)
	s2 =  c(1,1,1,1,.5,.5,.5)
	s1t = c(1,2,3,4,5  ,6,7)
	s2t = c(1,2,3,4,5.5,6,7)	
	res = ns_area_between_km(s1t,s2t,s1,s2);
	if (res != .25){
		print (paste("Error: ",res,"!= .25"));
	}else print("Test 3: OK");

	
	s1 =  c(.5,.5,.5,.5,0,0,0,0,0)
	s2 =  c(1 ,.5,.5,.5,0,0,0,0,0)
	s1t = c(1,2,3,4,5,6,7,8,9)
	s2t = s1t;
	res = ns_area_between_km(s1t,s2t,s1,s2);
	if (res != .5){
		print (paste("Error: ",res,"!= .5"));
	}else print("Test 4: OK");
	
	s1 =  c(.5,.5,.5,.5,0)
	s2 =  c(1 ,.5,.5,.5,0)
	s1t = c(1,2,3,4,5)
	s2t = s1t;
	res = ns_area_between_km(s1t,s2t,s1,s2);
	if (res != .5){
		print (paste("Error: ",res,"!= .5"));
	}else print("Test 5: OK");
	
	s1 =  c(1,.5,.5,.5)
	s2 =  c(1 ,.5,.5,.5,0)
	s1t = c(1,2,3)
	s2t = c(1,2,3,4)
	res = ns_area_between_km(s1t,s2t,s1,s2);
	if (res != 0){
		print (paste("Error: ",res,"!= 0"));
	}else print("Test 6: OK");
}