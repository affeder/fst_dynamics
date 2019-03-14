
foreach(nval = c(30, 100, 500))%do%{
    
trueMu <- .00001
Nval <- 1/trueMu

#vary n
nvar <- tbl_df(expand.grid(
    s = c(.3, 1, 3),
    #on cluster, length.out = 10
    m = lseq(1e-04, 1e-01, length.out = 10),
    n = c(nval),
    #on cluster, i = 1:200
    i = 1:2,
    offs = c(30),
    addhund = c(1),
    N = 1/trueMu,
    #on cluster, tol = 0.001
    tol = .1,
    runtype = "varyn"))

par.true.full <- nvar %>% arrange(n)

len <- 100
sampGens <- c(1:4, seq(5, 100, by = 5))

truths.sim <- foreach(s = par.true.full$s,
                      m = par.true.full$m,
                      i = par.true.full$i,
		      n = par.true.full$n,
		      N = par.true.full$N,
		      runtype = par.true.full$runtype,
                      .combine='rbind',
                      .packages = c('foreach', 'dplyr', 'abc')) %do% {
    toRet <- subSim(N, N, trueMu, s, s, m, len, relgens = sampGens, recenterCond = NA) %>%
             do(sampleMe(., n))  %>%
             do(computeStats(., paste("samp", sampGens, sep = "")))  %>%
	     mutate(N = N, n = n, s1 = s, m = m, i = i, runtype = runtype)
    toRet
 }                                                              

stat.true.full <- tbl_df(truths.sim)
par.true <- stat.true.full %>% select(N, s1, m, i) %>%
    mutate(theta = N*trueMu) %>% select(theta, s1, m, i)

par.true.full <- par.true.full %>% mutate(s1 = s)

stat.true.full <- left_join(tbl_df(truths.sim), par.true.full,
          by = c("N", "n", "s1", "m", "i", "runtype"))

tmp <- readInPrior(c(1:4, seq(5, 100, by = 5)), n = nval, nos = FALSE)

stat.sim.full <- tmp[[1]]
par.sim <- tbl_df(tmp[[2]])
svals <- tbl_df(tmp[[3]])

foreach(i = 1:nrow(stat.true.full), .combine = "rbind")%do%{

    print(i)
    tmp.true <- (stat.true.full %>% select(-N, -n, -s1, -m, -i, -runtype, -offs, -addhund, -tol, -s))[i,]
    par.true <- (stat.true.full %>% mutate(theta = N*trueMu) %>% select(theta, s1, m))[i,] %>% gather("var", "true")

    offs <- stat.true.full[i,]$offs
    addhund <- stat.true.full[i,]$addhund
    tolval <- stat.true.full[i,]$tol

    #compute stats
    firsttp <- sampGens[min(which((tmp.true %>% select(matches('percDR'))) >= .95))]

    selGens <- c(5, firsttp, firsttp + 30)
    relGens <- c(5, firsttp, firsttp + offs)
 
    if(addhund == 1){
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

    abcout <- abcProc.l(stats.true, stat.sim.sub, par.sim, par.true, tols = tolval, 
      	  quants = c(seq(0, 1, by= .05), .025, .975)) 

    stats.to.add <- (stat.true.full %>% select(N, n, s1, m, i, runtype, offs, addhund, tol, s))[i,]
    abc.to.add <- abcout %>% ungroup()

    toPrint <- tbl_df(cbind(stats.to.add, abc.to.add))

    write.table(toPrint, paste0("../../out/sumStats-varyn-n",nval,".txt"),
            row.names = FALSE,
            col.names = FALSE, quote = FALSE, append=TRUE )

}

}


