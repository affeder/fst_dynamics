
nval <- 100
trueMu <- .00001
Nval <- 1/trueMu

#vary n
nvar <- tbl_df(expand.grid(
    s = c(1),
    m = lseq(1e-04, 1e-01, length.out = 10),
    n = c(nval),
    asym = c(0.1, 0.5, 1),
   #on cluster, i = 1:100
    i = 1:2,
    offs = c(30),
    addhund = c(1),
    N = 1/trueMu,
    #on cluster, tol = 0.001
    tol = .1,
    runtype = c("m.asym", "N.asym", "s.asym")))

par.true.full <- nvar %>% arrange(n)

len <- 100
sampGens <- seq(5, 100, by = 5)

truths.sim <- foreach(s = par.true.full$s,
                      m = par.true.full$m,
                      i = par.true.full$i,
		      n = par.true.full$n,
		      N = par.true.full$N,
		      asym = par.true.full$asym,
		      runtype = par.true.full$runtype,
                      .combine='rbind',
                      .packages = c('foreach', 'dplyr', 'abc')) %do% {

    if(runtype == "N.asym"){
    	toRet <- subSim.asym(N, asym*N, trueMu, s, s, m, m, len, relgens = sampGens, recenterCond = NA) %>%
             do(sampleMe(., n))  %>%
             do(computeStats(., paste("samp", sampGens, sep = "")))  %>%
	     mutate(N = N, n = n, s1 = s, m = m, i = i, runtype = runtype, asym = asym)
	     }

    if(runtype == "m.asym"){
    toRet <- subSim.asym(N, N, trueMu, s, s, m,asym*m, len, relgens = sampGens, recenterCond = NA) %>%
             do(sampleMe(., n))  %>%
             do(computeStats(., paste("samp", sampGens, sep = "")))  %>%
	     mutate(N = N, n = n, s1 = s, m = m, i = i, runtype = runtype, asym = asym)
	     }

    if(runtype == "s.asym"){
    	toRet <- subSim.asym(N, N, trueMu, s, asym*s, m, m, len, relgens = sampGens, recenterCond = NA) %>%
             do(sampleMe(., n))  %>%
             do(computeStats(., paste("samp", sampGens, sep = "")))  %>%
	     mutate(N = N, n = n, s1 = s, m = m, i = i, runtype = runtype, asym = asym)
	     }

    toRet
 }                                                              



stat.true.full <- tbl_df(truths.sim)
par.true <- stat.true.full %>% select(N, s1, m, i, asym, runtype) %>%
    mutate(theta = N*trueMu) %>% select(theta, s1, m, i, asym, runtype)

par.true.full <- par.true.full %>% mutate(s1 = s)

stat.true.full <- left_join(tbl_df(truths.sim), par.true.full,
          by = c("N", "n", "s1", "m", "i", "asym", "runtype"))

## This has been read in previously
## tmp <- readInPrior(seq(5, 100, by = 5), n = nval)
## stat.sim.full <- tmp[[1]]
## par.sim <- tbl_df(tmp[[2]])
## svals <- tbl_df(tmp[[3]])

foreach(i = 1:nrow(stat.true.full), .combine = "rbind")%do%{

    print(i)
    tmp.true <- (stat.true.full %>% select(-N, -n, -s1, -m, -i, -runtype, -offs, -addhund, -tol, -s, -asym))[i,]
    par.true <- (stat.true.full %>% mutate(theta = N*trueMu) %>% select(theta, s1, m))[i,] %>% gather("var", "true")

    offs <- stat.true.full[i,]$offs
    addhund <- stat.true.full[i,]$addhund
    tolval <- stat.true.full[i,]$tol

    #compute stats

    firsttp <- sampGens[min(which((tmp.true %>% select(matches('percDR'))) >= .95))]


    if(is.na(firsttp) | firsttp + offs > 100){
        badadd <- (stat.true.full %>% select(N, n, s1, m, i, runtype, offs, addhund, tol, s, asym))[i,] %>%
    		 mutate(val = NA, q = NA, var = NA, med.dist = NA, mean.dist = NA, med.distsq = NA,
		 mean.distsq = NA)

     write.table(badadd, paste0("../../out/sumStats-asym-n",nval,"-long-ah.txt"),
            row.names = FALSE,
            col.names = FALSE, quote = FALSE, append=TRUE )

     }else{

    selGens <- c(5, firsttp, firsttp + 30)#20
    relGens <- c(5, firsttp, firsttp + offs)

    if(addhund == 1 & firsttp + offs != 100){
        relGens <- c(relGens, 100)
    }

    relGenReg <- paste0(paste0("_",relGens, "$"), collapse = "|")
    stats.true <- tmp.true %>% select(matches(relGenReg)) 

    #compute s
    relperc <- as.matrix(stats.true %>% select(matches('percDR')) %>% rowwise())
    svalsToAdd <- apply(relperc, 1, bespokeMinimization, relGens)
    stats.true <- tbl_df(stats.true %>% mutate(bestS = svalsToAdd))

    #Here, subset the true matrix
    sToAdd <- c(as.matrix((svals %>% select(paste0(selGens, collapse = "_")))))
    stat.sim.sub <- stat.sim.full %>% select(matches(relGenReg)) %>% 
         		mutate(bestS = sToAdd)

    abcout <- abcProc.q(stats.true, stat.sim.sub, par.sim, par.true, tols = tolval, 
      	  quants = c(seq(0, 1, by= .05), .025, .975)) 

    stats.to.add <- (stat.true.full %>% select(N, n, s1, m, i, runtype, offs, addhund, tol, s, asym))[i,]
    abc.to.add <- abcout %>% ungroup()

    toPrint <- tbl_df(cbind(stats.to.add, abc.to.add))

     write.table(toPrint, paste0("../../out/sumStats-asym-n",nval,"-long.txt"),
            row.names = FALSE,
            col.names = FALSE, quote = FALSE, append=TRUE )

	    }
}




