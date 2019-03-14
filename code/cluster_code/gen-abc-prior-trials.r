#Explore timing and stat subsetting

#Let's precompute a bunch of different trials we're interested in with full info
#i.e., all time points and all statistics
trueMu <- .00001

#Number of trials:
#Note: on cluster, numtrials = 2*10^6 
# 1000 unparallelized trials took 28 minutes on a late 2013 iMac 
numtrials <- 10000

ms <- 10^runif(numtrials, -5, -.3)
thetas <- 10^runif(numtrials, -1, 1)
ss <-  10^runif(numtrials, -1, 2)

#thetas <- 10^runif(numtrials, -1, 0.2787536)
#ss <-  10^runif(numtrials, -0.1972263, 1.5)
trialnums <- 1:numtrials

len <- 100

sampGens = c(1:4, seq(5, 100, by = 5))

ptime <- system.time(foreach( theta = thetas,
                      s = ss,
                      m = ms,
                      trialnum = trialnums, 
                      .combine='rbind',
                      .packages = c('foreach', 'dplyr'),
                      .errorhandling = 'remove') %do% {
       toRet <- subSim(floor(theta/trueMu), floor(theta/trueMu), trueMu,
         s, s, m, len, relgens = sampGens, recenterCond = NA) %>%
           mutate(N = floor(theta/trueMu), s1 = s, m = m, trial = trialnum)
       conn <- file(sprintf("../../raw/prior/nf_output_%d.txt" , 
       	    Sys.getpid()), open = "a" )  
       write.table(toRet, conn, col.names=FALSE, row.names = FALSE)
       close(conn)
   })
print(ptime)




