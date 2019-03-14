

n <- 30
relGens <- c(1:4, seq(5, 100, by = 5))
sampGens <- c(5, 20, 50, 100)
trueMu <- 10^(-5)


stat.sim.full <- read.table(paste("../raw/stats/n",n,"_prior.txt", sep = ""), sep = " ", header = FALSE, nrows = 1500000)
stat.sim.full <- tbl_df(stat.sim.full)
    colnamesToAdd <-  c(apply(expand.grid(c("diffInHets", "FST", "GSTp",
                                        "thetaEst", "percDR", "numShared"),
                                          paste(relGens)),
                 1, paste, collapse = "_"), "N", "s1", "m", "n", "i")
    names(stat.sim.full) <- colnamesToAdd
    par.sim <- stat.sim.full %>% 
        select(N, s1, m) %>% mutate(theta = N*trueMu) %>% select(theta, s1, m)

#read in s data (easier to just store this separately)   
    bestS <- read.table(paste("../raw/svals/svals_n",n,"_gen5.txt",
                                      sep = ""), sep = " ", header = FALSE, nrows = 1500000)


    seqRow <- seq(10, 70, 5)
    relGensIt <- foreach(i = 1:length(seqRow))%do%{
        c(5, seqRow[i],seqRow[i] + 30 )
    }

    colnames(bestS) <- unlist(lapply(relGensIt, function(x){paste0(x, collapse = "_")}))

bestS <- tbl_df(bestS)

stat.sim.sub <- stat.sim.full %>% subSetDat(sampGens) %>%
    mutate(bestS = as.vector(bestS)$`5_20_50`)

#In the actual run for the paper, I here import a second prior so the sample size remains large for the second fit.
# I do not generate this in the sample code, but I hope to update this shortly
## r2 <- tbl_df(read.table(paste("../out/prior_n30_r2.txt", sep = ""),header = FALSE))
## relGens <- c(5, 20, 50, 100)
## colnamesToAdd <-  c(apply(expand.grid(c("diffInHets", "FST", "GSTp",
##                                         "thetaEst", "percDR", "numShared"),
##                                           paste(relGens)),
##                  1, paste, collapse = "_"), "N", "s1", "m", "n", "i")
## names(r2) <- colnamesToAdd
## par.r2 <- r2 %>% 
##     select(N, s1, m) %>% mutate(theta = N*trueMu) %>% select(theta, s1, m)
## #New posterior
## r2_ep <- tbl_df(read.table(paste("../d_out/n30_prior_ep.txt", sep = ""),header = FALSE))
## sampgens <- c(1:4, seq(5, 100, by = 5))

## colnamesToAdd <-  c(apply(expand.grid(c("diffInHets", "FST", "GSTp",
##                                         "thetaEst", "percDR", "numShared"),
##                                           paste(sampgens)),
##                  1, paste, collapse = "_"), "N", "s1", "m", "n", "i")
## names(r2_ep) <- colnamesToAdd
## par.r2_ep<- r2_ep%>% 
##     select(N, s1, m) %>% mutate(theta = N*trueMu) %>% select(theta, s1, m)

## bestS <- tbl_df(read.table(paste("../d_out/svals_n30_ep.txt",
##                                       sep = ""), sep = " ", header = FALSE))

## colnames(bestS) <- unlist(lapply(relGensIt, function(x){paste0(x, collapse = "_")}))

## r2_ep <- r2_ep %>% subSetDat(sampGens) %>%
##     mutate(bestS = as.vector(bestS)$`5_20_40`) %>%
##     mutate(theta = par.r2_ep$theta, 
##            s1 = par.r2_ep$s1,
##            m = par.r2_ep$m)




realDat <- function(inDat, tol1 = .01, tol2 = .01){

formatted <- inDat %>% mutate(pop = ifelse(loc == loc[1], "A", "B")) %>%
    group_by(pop) %>%
    mutate(mut = ifelse(mut != "WT", 
               paste0(which(mut == allmuttypes), "-C"), "WT" )) %>%
    select(pop, mut, sampNames) %>% ungroup()

    sampGens <- c(5, 20, 50, 100)
    sampNames <- paste("samp", sampGens, sep = "")
    num.ref <- formatted %>% computeStats(sampNames) 
    relPerc <- num.ref %>% select(matches('percDR'))
    stat.true <- num.ref %>% mutate(bestS = bespokeMinimization(relPerc, sampGens))
    return(abcProc.r2(stat.true, stat.sim.sub, par.sim, bigRet = TRUE, tol = tol1, tol2 = tol2))
}


abcProc.r2 <- function(true.stats.ref, sim.stats.ref, sim.par.ref, bigRet = FALSE, tol1 = .01, tol2 = .01){

#Which quantiles do we want returned
quant.rets <- c(.025, .25,  .5, .75, .975)
#How big should the first round posteriors be?
r1.posteriors <- c(.025, .975)
#First round tolerance
r1.tolerance <- tol1
#Second round tolerance    
r2.tolerance <- tol2
#true stats ref 1 by numstats * numtimepoints
true.stats.ref
#sim stats ref x by numstats * numtimepoints
sim.stats.ref
#sim par ref # x by 3
sim.par.ref
#sim ref
sim.ref <- cbind(sim.par.ref, sim.stats.ref)

round1params <- "thetaEst|bestS"
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
    #And only match round2params for the prior and the trial itself

    par.sim.r2 <- fullsims %>% select(theta, s1, m)
    sub.sim.r2 <- fullsims %>% select(matches(round2params))

    par.sim.r2.ep <- r2_ep %>% select(theta, s1, m)
    sub.sim.r2.ep <- r2_ep %>% select(matches(round2params))

    ## par.sim.r2.t1 <- r2 %>% select(theta, s1, m)
    ## sub.sim.r2.t1 <- r2 %>% select(matches(round2params))
   
    par.sim.r2 <- rbind(par.sim.r2, par.sim.r2.ep)#, par.sim.r2.t1)
    sub.sim.r2 <- rbind(sub.sim.r2, sub.sim.r2.ep)#, sub.sim.r2.t1)
    print(dim(par.sim.r2))


    sub.true.r2 <- true.stats.ref %>% select(matches(round2params))
#sub.sim.r2 <- r2 %>% select(matches(round2params))
#par.sim.r2 <- par.r2 %>% select(theta, s1, m)
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
    return(list(all.r1, all.r2, par.sim.r2, fullsims))
}else{
    return(toRet)
}
}


realDat.no_r2 <- function(inDat, tol1 = .01, tol2 = .01){

formatted <- inDat %>% mutate(pop = ifelse(loc == loc[1], "A", "B")) %>%
    group_by(pop) %>%
    mutate(mut = ifelse(mut != "WT", 
               paste0(which(mut == allmuttypes), "-C"), "WT" )) %>%
    select(pop, mut, sampNames) %>% ungroup()

    sampGens <- c(5, 20, 50, 100)
    sampNames <- paste("samp", sampGens, sep = "")
    num.ref <- formatted %>% computeStats(sampNames) 
    relPerc <- num.ref %>% select(matches('percDR'))
    stat.true <- num.ref %>% mutate(bestS = bespokeMinimization(relPerc, sampGens))
    return(abcProc.no_r2(stat.true, stat.sim.sub, par.sim, bigRet = TRUE, tol = tol1, tol2 = tol2))
}


abcProc.no_r2 <- function(true.stats.ref, sim.stats.ref, sim.par.ref, bigRet = FALSE, tol1 = .01, tol2 = .01){

#Which quantiles do we want returned
quant.rets <- c(.025, .25,  .5, .75, .975)
#How big should the first round posteriors be?
r1.posteriors <- c(.025, .975)
#First round tolerance
r1.tolerance <- tol1
#Second round tolerance    
r2.tolerance <- tol2
#true stats ref 1 by numstats * numtimepoints
true.stats.ref
#sim stats ref x by numstats * numtimepoints
sim.stats.ref
#sim par ref # x by 3
sim.par.ref
#sim ref
sim.ref <- cbind(sim.par.ref, sim.stats.ref)

round1params <- "thetaEst|bestS"
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
    #And only match round2params for the prior and the trial itself

    par.sim.r2 <- fullsims %>% select(theta, s1, m)
    sub.sim.r2 <- fullsims %>% select(matches(round2params))

#    par.sim.r2.ep <- r2_ep %>% select(theta, s1, m)
#    sub.sim.r2.ep <- r2_ep %>% select(matches(round2params))

    ## par.sim.r2.t1 <- r2 %>% select(theta, s1, m)
    ## sub.sim.r2.t1 <- r2 %>% select(matches(round2params))
   
#    par.sim.r2 <- rbind(par.sim.r2, par.sim.r2.ep)#, par.sim.r2.t1)
#    sub.sim.r2 <- rbind(sub.sim.r2, sub.sim.r2.ep)#, sub.sim.r2.t1)
    print(dim(par.sim.r2))


    sub.true.r2 <- true.stats.ref %>% select(matches(round2params))
#sub.sim.r2 <- r2 %>% select(matches(round2params))
#par.sim.r2 <- par.r2 %>% select(theta, s1, m)
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
    return(list(all.r1, all.r2, par.sim.r2, fullsims))
}else{
    return(toRet)
}
}



## realDat.inv <- function(runType, tol = .01){

##     sampGens <- c(5, 20, 50, 100)
##     trueDat <- readInReal(runType)
##     sampNames <- paste("samp", sampGens, sep = "")
##     num.ref <- trueDat %>% computeStats(sampNames) 
##     relPerc <- num.ref %>% select(matches('percDR'))
##     stat.true <- num.ref %>% mutate(bestS = bespokeMinimization(relPerc, sampGens))
##     return(abcProc.r2.inv(stat.true, stat.sim.sub, par.sim, bigRet = TRUE, tol = tol))
## }

## abcProc.r2.inv <- function(true.stats.ref, sim.stats.ref, sim.par.ref, bigRet = FALSE, tol = .01){

## #Which quantiles do we want returned
## quant.rets <- c(.025, .25,  .5, .75, .975)
## #How big should the first round posteriors be?
## r1.posteriors <- c(.025, .975)
## #First round tolerance
## r1.tolerance <- tol
## #Second round tolerance    
## r2.tolerance <- tol
## #true stats ref 1 by numstats * numtimepoints
## true.stats.ref
## #sim stats ref x by numstats * numtimepoints
## sim.stats.ref
## #sim par ref # x by 3
## sim.par.ref
## #sim ref
## sim.ref <- cbind(sim.par.ref, sim.stats.ref)

## round1params <- "thetaEst|bestS"
## round2params <- "FST|diffInHets|GSTp|numShared"

## #first, we'll subset it to round 1 parameters
## sub.true.r1 <- true.stats.ref %>% select(matches(round1params))
## sub.sim.r1 <-  sim.stats.ref %>% select(matches(round1params))

## #Fit the first ABC
## abcfits.r1 <- abc(target = sub.true.r1, 
##                   param = par.sim, 
##                   sumstat = sub.sim.r1, 
##                   tol = r1.tolerance, method = "rejection")
## all.r1 <- tbl_df(abcfits.r1$unadj.values) %>%
##     gather(var, val)

##     #Determine the posteriors for s and m
##     lims.r1 <- all.r1 %>% group_by(var) %>%
##         summarize(q.l = quantile(val, r1.posteriors[1]), 
##                   q.h = quantile(val, r1.posteriors[2])) %>% mutate(round = 1)
##     theta.q.l <- (lims.r1 %>% filter(var == "theta") %>% select(q.l))
##     theta.q.h <- (lims.r1 %>% filter(var == "theta") %>% select(q.h))
##     s1.q.l <- (lims.r1 %>% filter(var == "s1") %>% select(q.l))
##     s1.q.h <- (lims.r1 %>% filter(var == "s1") %>% select(q.h))
##     #restrict our prior to these limits
##     fullsims <- sim.ref %>% 
##         filter(theta >= (theta.q.l$q.l) & theta <= (theta.q.h$q.h)) %>%
##         filter(s1 >= (s1.q.l$q.l) & s1 <= (s1.q.h$q.h))
##     #And only match round2params for the prior and the trial itself
## #    par.sim.r2 <- fullsims %>% select(theta, s1, m)
## #    sub.sim.r2 <- fullsims %>% select(matches(round2params))
## sub.true.r2 <- true.stats.ref %>% select(matches(round2params))
## sub.sim.r2 <- r2 %>% select(matches(round2params))
## par.sim.r2 <- par.r2 %>% select(theta, s1, m)
##     abcfits.r2 <- abc(target = sub.true.r2, 
##                       param = par.sim.r2, 
##                       sumstat = sub.sim.r2, 
##                       tol = r2.tolerance, method = "rejection")
##     all.r2 <- tbl_df(abcfits.r2$unadj.values)
##     return(abcfits.r2)
## }
