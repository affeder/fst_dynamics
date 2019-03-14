#Let's precompute a bunch of different trials we're interested in with full info
#i.e., all time points and all statistics
trueMu <- .00001

pathToFiles <- "../../raw/prior/"
#Here's our list of files
fileCheck <- list.files(pathToFiles)

sampGens = c(1:4, seq(5, 100, by = 5))

nvals <- c(30, 100, 500)

foreach(fopen = fileCheck,
        .packages = c('foreach', 'dplyr'),
        .errorhandling = 'remove') %do% {
            #Create path to file
            file <- paste(pathToFiles, fopen, sep = "")
            ## Create connection
            con.in <- file(description=file, open="r")
            com <- paste("cut -d' ' -f31 ", file) 
            trialIDs <- system(command = com, intern = TRUE)
            trialnums <- table(as.numeric(trialIDs))
            #Ok, so, using trialIDs as a guide, we're going to read in one
            #trial at a time. 
            for(i in 1:(length(trialnums)-1)){
                #Read in just enough to compute one trials worth of stats
                #we do it this way because at scale, we quickly exceed mem limits
                tmp <- scan(file=con.in, nlines=trialnums[i], quiet=TRUE, 
                            sep = " ", what = "string")
                #Format those things into a format that we can compute with
                toAdd <- tbl_df(matrix(tmp, nrow = trialnums[i], byrow = TRUE))
                names(toAdd) <- c("pop", "mut", paste("gen", sampGens, sep = ""),
                                  "rec", "N", "s1", "m", "trialnum")
                toAdd <- toAdd %>% mutate_each(funs(as.numeric),
                                               matches("gen|^N$|s1|^m$|trialnum")) %>%
                    mutate(theta = N*trueMu)
                #Compute
		for(n in nvals){
		    toPrint <- toAdd %>% do(sampleMe(., n)) %>%
                    	    do(computeStats(., paste("samp", sampGens, sep = ""))) %>%
			    mutate(N = unique(toAdd$N), s1 = unique(toAdd$s1),
                            m = unique(toAdd$m), n = n, trial = unique(toAdd$trialnum))
                #Print
			con.out <- file(
	sprintf("../../raw/stats/n%d/stats_output_%d.txt" ,
	n, Sys.getpid()), open = "a" )  
                	write.table(toPrint, con.out, col.names=FALSE, row.names = FALSE)
			close(con.out)
			}
   		}
}



