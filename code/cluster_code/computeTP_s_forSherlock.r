
nvals <- c(30, 100, 500)

trueMu <- .00001
relGens <- c(1:4, seq(5, 100, by = 5))
sampGens <- c(1:4, seq(5, 100, by = 5))


foreach(n = nvals, .packages = c('foreach', 'dplyr', 'abc'),
                    .errorhandling = 'remove')%do%{

            pathToFiles <- paste("../../raw/stats/n",n,"/", sep = "")

            #Here's our list of files
            fileCheck <- list.files(pathToFiles)

         foreach(fopen = fileCheck,
                 .packages = c('foreach', 'dplyr', 'abc'),
                 .errorhandling = 'remove') %do% {

	#read in data
                     stat.sim.full <- read.table(paste(pathToFiles, fopen, sep = ""),
                                                 sep = " ", header = FALSE)
                     stat.sim.full <- tbl_df(stat.sim.full)
                     colnamesToAdd <-  c(apply(expand.grid(c("diffInHets", "FST", "GSTp",
                                        "thetaEst", "percDR", "numShared"),
                                        paste(relGens)),
                                        1, paste, collapse = "_"), "N", "s1", "m", "n", "i")
                     names(stat.sim.full) <- colnamesToAdd
                     par.sim <- stat.sim.full %>% 
                         select(N, s1, m) %>% mutate(theta = N*trueMu) %>% 
                         select(theta, s1, m)

                     foreach(gen = c(2, 5), .packages = c('foreach', 'dplyr', 'abc'),
                             .errorhandling = 'remove')%do%{

##################################################
#effect of sampling time
##################################################

                        seqRow <- seq(10, 70, 5)
                        relGensIt <- foreach(i = 1:length(seqRow))%do%{
                            c(gen, seqRow[i],seqRow[i] + 30 )
                        }

#################################################
# Here, we'll precompute all the bestS values and
# store them in a place where we can retrieve them
#################################################

#I wrote a twostep minimization function that first coarsely searches
#a grid and then uses R's optimize once it has narrowed the search space
 ## bespokeMinimization <- function(freqs, sampGens){

 ##     if(sum(as.numeric(freqs)) == 0){ return(0) }
 ##         coarse <- c(seq(.01, 1, by = .01), seq(2, 20, by = 1))
 ## 	coarsefits <- sAtGensVect(coarse, freqs, relGens = sampGens)
 ## 	fitval <- coarse[min(which(coarsefits == min(coarsefits)))]
 ## 	return(optimize(sAtGensVect, c(fitval * .5, fitval * 2), 
 ## 				     freqs, relGens = sampGens)$minimum)  

 ## }

                        allSVals <- foreach(i = 1:length(relGensIt), .combine = "cbind") %do% {
                            sampGens = relGensIt[[i]]
                            bestSvect <-  stat.sim.full %>% subSetDat(sampGens) %>%
                                select(matches('percDR')) %>%
                                rowwise()  %>%
                                do(tbl_df(bespokeMinimization(freqs = ., sampGens = sampGens))) %>%
                                ungroup()
                            names(bestSvect) <- paste("bestS_", paste(relGensIt[[i]], collapse = "-"), sep = "")
                            return(bestSvect)
                        }


                        con.out <- file(
                            paste("../../raw/svals/n",n,"/sout_",gsub("stats_output_|.txt", "", fopen),"-gen", gen,".txt", sep = ""), open = "a" )  
                                 write.table(allSVals, con.out, col.names=FALSE, row.names = FALSE)
                                 close(con.out)
                        
                             }
                 }
}                    
