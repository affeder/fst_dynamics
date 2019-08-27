#ABC procedure
require(stringr)

abcProc <- function(true.stats.ref, sim.stats.ref, sim.par.ref, bigRet = FALSE, tols = .01){

#Which quantiles do we want returned
quant.rets <- c(.025, .25,  .5, .75, .975)
#How big should the first round posteriors be?
r1.posteriors <- c(.025, .975)
#First round tolerance
r1.tolerance <- tols
#Second round tolerance    
r2.tolerance <- tols
#true stats ref 1 by numstats * numtimepoints
true.stats.ref
#sim stats ref x by numstats * numtimepoints
sim.stats.ref
#sim par ref # x by 3
sim.par.ref
#sim ref
sim.ref <- cbind(sim.par.ref, sim.stats.ref)

round1params <- "thetaEst|bestS"
round2params <- "FST|diffInHets|GSTp"

#first, we'll subset it to round 1 parameters
sub.true.r1 <- true.stats.ref %>% select(matches(round1params))
sub.sim.r1 <-  sim.stats.ref %>% select(matches(round1params))

#Fit the first ABC
abcfits.r1 <- abc(target = sub.true.r1, 
                  param = par.sim, 
                  sumstat = sub.sim.r1, 
                  tol = r1.tolerance, method = "rejection")
all.r1 <- tbl_df(abcfits.r1$unadj.values) %>%
    gather(var, val)

    #Determine the posteriors for s and m
    lims.r1 <- all.r1 %>% group_by(var) %>%
        summarize(q.l = quantile(val, r1.posteriors[1]), 
                  q.h = quantile(val, r1.posteriors[2])) %>% mutate(round = 1)
    theta.q.l <- (lims.r1 %>% filter(var == "theta") %>% select(q.l))
    theta.q.h <- (lims.r1 %>% filter(var == "theta") %>% select(q.h))
    s1.q.l <- (lims.r1 %>% filter(var == "s1") %>% select(q.l))
    s1.q.h <- (lims.r1 %>% filter(var == "s1") %>% select(q.h))
    #restrict our prior to these limits
    fullsims <- sim.ref %>% 
        filter(theta >= (theta.q.l$q.l) & theta <= (theta.q.h$q.h)) %>%
        filter(s1 >= (s1.q.l$q.l) & s1 <= (s1.q.h$q.h))
    par.sim.r2 <- fullsims %>% select(theta, s1, m)
    #And only match round2params for the prior and the trial itself
    sub.sim.r2 <- fullsims %>% select(matches(round2params))
    sub.true.r2 <- true.stats.ref %>% select(matches(round2params))
    abcfits.r2 <- abc(target = sub.true.r2, 
                      param = par.sim.r2, 
                      sumstat = sub.sim.r2, 
                      tol = r2.tolerance, method = "rejection")
    all.r2 <- tbl_df(abcfits.r2$unadj.values) %>%
                 gather(var, val)
#Determine the posteriors for s and m
    lims.toRet.r1 <- all.r1 %>% group_by(var) %>%
        summarize(q.025 = quantile(val, c(.025)),
                  q.25 = quantile(val, c(.25)), 
                  q.5 = quantile(val, c(.5)),
                  q.75 = quantile(val, c(.75)),
                  q.975 = quantile(val, c(.975))) %>% mutate(round = 1)
    lims.toRet.r2 <- all.r2 %>% group_by(var) %>%
        summarize(q.025 = quantile(val, c(.025)),
                  q.25 = quantile(val, c(.25)), 
                  q.5 = quantile(val, c(.5)),
                  q.75 = quantile(val, c(.75)),
                  q.975 = quantile(val, c(.975))) %>% mutate(round = 2)
    toRet <- rbind(lims.toRet.r1, lims.toRet.r2)
    names(toRet) <- c("var", "q.025", "q.25", 
                      "q.50", "q.75", "q.975", "round")
if(bigRet == TRUE){
    return(list(all.r1, all.r2, par.sim.r2))
}else{
    return(toRet)
}
}




abcProc.q <- function(true.stats.ref, sim.stats.ref, sim.par.ref, par.true, 
	  bigRet = FALSE, tols = .01, quants = c(.025, .25, .5, .75, .975)){

#Which quantiles do we want returned

#How big should the first round posteriors be?
r1.posteriors <- c(.025, .975)
#First round tolerance
r1.tolerance <- tols
#Second round tolerance
r2.tolerance <- tols
#true stats ref 1 by numstats * numtimepoints
true.stats.ref
#sim stats ref x by numstats * numtimepoints
sim.stats.ref
#sim par ref # x by 3
sim.par.ref
#sim ref
sim.ref <- cbind(sim.par.ref, sim.stats.ref)

round1params <- "thetaEst|bestS"
round2params <- "FST|diffInHets|GSTp"

#first, we'll subset it to round 1 parameters
sub.true.r1 <- true.stats.ref %>% select(matches(round1params))
sub.sim.r1 <-  sim.stats.ref %>% select(matches(round1params))

#Fit the first ABC
abcfits.r1 <- abc(target = sub.true.r1,
                  param = par.sim,
                  sumstat = sub.sim.r1,
                  tol = r1.tolerance, method = "rejection")
		  all.r1 <- tbl_df(abcfits.r1$unadj.values) %>%
		  gather(var, val)

    #Determine the posteriors for s and m
        lims.r1 <- all.r1 %>% group_by(var) %>%
	        summarize(q.l = quantile(val, r1.posteriors[1]),
	                  q.h = quantile(val, r1.posteriors[2])) %>% mutate(round = 1)
      theta.q.l <- (lims.r1 %>% filter(var == "theta") %>% select(q.l))
      theta.q.h <- (lims.r1 %>% filter(var == "theta") %>% select(q.h))
      s1.q.l <- (lims.r1 %>% filter(var == "s1") %>% select(q.l))
      s1.q.h <- (lims.r1 %>% filter(var == "s1") %>% select(q.h))
  #restrict our prior to these limits
      fullsims <- sim.ref %>%
            filter(theta >= (theta.q.l$q.l) & theta <= (theta.q.h$q.h)) %>%
            filter(s1 >= (s1.q.l$q.l) & s1 <= (s1.q.h$q.h))
     par.sim.r2 <- fullsims %>% select(theta, s1, m)
   #And only match round2params for the prior and the trial itself
      sub.sim.r2 <- fullsims %>% select(matches(round2params))
      sub.true.r2 <- true.stats.ref %>% select(matches(round2params))
      abcfits.r2 <- abc(target = sub.true.r2,
      param = par.sim.r2,
      sumstat = sub.sim.r2,
      tol = r2.tolerance, method = "rejection")
    all.r2 <- tbl_df(abcfits.r2$unadj.values) %>%
             gather(var, val)
    #Determine the posteriors for s and m

#par.true <- par.true %>% gather(var, true) %>% filter(var != "i")


avdists.r1 <- all.r1 %>% filter(str_detect(var, 's1|theta')) %>%
    left_join(par.true, by = 'var') %>% 
    mutate(dist = abs(val - true), distsq = (val-true)^2) %>%
    group_by(var) %>% 
    summarize(med.dist = median(dist), mean.dist = mean(dist), 
   	      med.distsq = median(distsq), mean.distsq = mean(distsq))
avdists.r2 <- all.r2 %>% filter(var == "m") %>%
    left_join(par.true, by = 'var') %>% 
    mutate(dist = abs(val - true), distsq = (val-true)^2) %>%
    group_by(var) %>% 
    summarize(med.dist = median(dist), mean.dist = mean(dist), 
   	      med.distsq = median(distsq), mean.distsq = mean(distsq))

retQuants <- function(mat, quants){
	  qs <- as.matrix(quantile((mat %>% select(val))$val, quants))	
	  colnames(qs) <- "tmp"
          tbl_df(qs) %>%
              mutate(q = quants, var = unique(mat$var))
}

    toRet.r1 <- all.r1 %>% group_by(var) %>%
        do(retQuants(., quants)) %>%
        filter(str_detect(var, 's1|theta'))
    toRet.r2 <- all.r2 %>% group_by(var) %>%
        do(retQuants(., quants)) %>%
        filter(str_detect(var, 'm'))
    toRet <- rbind(toRet.r1, toRet.r2)
    names(toRet) <- c("val", "q", "var")

toRet <- toRet %>%
    mutate(med.dist = ifelse(var == "m", avdists.r2$med.dist, avdists.r1$med.dist),
           mean.dist = ifelse(var == "m", avdists.r2$mean.dist, avdists.r1$mean.dist),
           med.distsq = ifelse(var == "m", avdists.r2$med.distsq, avdists.r1$med.distsq),
           mean.distsq = ifelse(var == "m", avdists.r2$mean.distsq, avdists.r1$mean.distsq))

    if(bigRet == TRUE){
    	      return(list(all.r1, all.r2, par.sim.r2))
    }else{
        return(toRet)
    }
}


abcProc.l <- function(true.stats.ref, sim.stats.ref, sim.par.ref, par.true, 
	  bigRet = FALSE, tols = .01, quants = c(.025, .25, .5, .75, .975)){

#Which quantiles do we want returned

#How big should the first round posteriors be?
r1.posteriors <- c(.025, .975)
#First round tolerance
r1.tolerance <- tols
#Second round tolerance
r2.tolerance <- tols
#true stats ref 1 by numstats * numtimepoints
true.stats.ref
#sim stats ref x by numstats * numtimepoints
sim.stats.ref
#sim par ref # x by 3
sim.par.ref
#sim ref
sim.ref <- cbind(sim.par.ref, sim.stats.ref)

#round1params <- "thetaEst|bestS"
#round2params <- "FST|diffInHets|GSTp"
round1params <- "thetaEst|percDR"
round2params <- "FST|diffInHets|GSTp|numShared"

#first, we'll subset it to round 1 parameters
sub.true.r1 <- true.stats.ref %>% select(matches(round1params))
sub.sim.r1 <-  sim.stats.ref %>% select(matches(round1params))

#Fit the first ABC
abcfits.r1 <- abc(target = sub.true.r1,
                  param = par.sim,
                  sumstat = sub.sim.r1,
                  tol = r1.tolerance, method = "rejection")
		  all.r1 <- tbl_df(abcfits.r1$unadj.values) %>%
		  gather(var, val)

    #Determine the posteriors for s and m
        lims.r1 <- all.r1 %>% group_by(var) %>%
	        summarize(q.l = quantile(val, r1.posteriors[1]),
	                  q.h = quantile(val, r1.posteriors[2])) %>% mutate(round = 1)
      theta.q.l <- (lims.r1 %>% filter(var == "theta") %>% select(q.l))
      theta.q.h <- (lims.r1 %>% filter(var == "theta") %>% select(q.h))
      s1.q.l <- (lims.r1 %>% filter(var == "s1") %>% select(q.l))
      s1.q.h <- (lims.r1 %>% filter(var == "s1") %>% select(q.h))
  #restrict our prior to these limits
      fullsims <- sim.ref %>%
            filter(theta >= (theta.q.l$q.l) & theta <= (theta.q.h$q.h)) %>%
            filter(s1 >= (s1.q.l$q.l) & s1 <= (s1.q.h$q.h))
     par.sim.r2 <- fullsims %>% select(theta, s1, m)
   #And only match round2params for the prior and the trial itself
      sub.sim.r2 <- fullsims %>% select(matches(round2params))
      sub.true.r2 <- true.stats.ref %>% select(matches(round2params))
      abcfits.r2 <- abc(target = sub.true.r2,
      param = par.sim.r2,
      sumstat = sub.sim.r2,
      tol = r2.tolerance, method = "rejection")
    all.r2 <- tbl_df(abcfits.r2$unadj.values) %>%
             gather(var, val)
    #Determine the posteriors for s and m

#par.true <- par.true %>% gather(var, true) %>% filter(var != "i")


avdists.r1 <- all.r1 %>% filter(str_detect(var, 's1|theta')) %>%
    left_join(par.true, by = 'var') %>% 
    mutate(dist = abs(log(val,10) - log(true,10)), 
    	   distsq = (log(val,10)-log(true,10))^2) %>%
    group_by(var) %>% 
    summarize(med.dist = median(dist), mean.dist = mean(dist), 
   	      med.distsq = median(distsq), mean.distsq = mean(distsq))
avdists.r2 <- all.r2 %>% filter(var == "m") %>%
    left_join(par.true, by = 'var') %>% 
    mutate(dist = abs(log(val,10) - log(true,10)), 
    	   distsq = (log(val,10)-log(true,10))^2) %>%
    group_by(var) %>% 
    summarize(med.dist = median(dist), mean.dist = mean(dist), 
   	      med.distsq = median(distsq), mean.distsq = mean(distsq))

retQuants <- function(mat, quants){
    tbl_df((as.matrix(quantile((mat %>% select(val))$val, quants)))) %>%
            mutate(q = quants, var = unique(mat$var))
	    }

toRet.r1 <- all.r1 %>% group_by(var) %>%
    do(retQuants(., quants)) %>%
        filter(str_detect(var, 's1|theta'))
	toRet.r2 <- all.r2 %>% group_by(var) %>%
	    do(retQuants(., quants)) %>%
	        filter(str_detect(var, 'm'))
		toRet <- rbind(toRet.r1, toRet.r2)
		names(toRet) <- c("val", "q", "var")

toRet <- toRet %>%
    mutate(med.dist = ifelse(var == "m", avdists.r2$med.dist, avdists.r1$med.dist),
           mean.dist = ifelse(var == "m", avdists.r2$mean.dist, avdists.r1$mean.dist),
           med.distsq = ifelse(var == "m", avdists.r2$med.distsq, avdists.r1$med.distsq),
           mean.distsq = ifelse(var == "m", avdists.r2$mean.distsq, avdists.r1$mean.distsq))

    if(bigRet == TRUE){
    	      return(list(all.r1, all.r2, par.sim.r2))
    }else{
        return(toRet)
    }
}







