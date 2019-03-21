#justsims.r is going to have functions JUST to run the simulation
#abc-functions.r is going to have all of the functions that deal with
#  computing statistics on the simulated runs, and helper plotting functions

source('justsims.r')



select <- dplyr::select

#These functions are needed to compute GST'
computeMax <- function(x, K){
    cb <- function(xval){
        return(xval - floor(xval))
    }

    jf <- function(xval){
        return(ceiling(1/xval))
    }
    x <- x*K
    toRet.gt <- (K*(K-1)-2*(K-1)*cb(x)*(1-cb(x))-floor(x)*(floor(x)-1) - 2*cb(x)*floor(x))/(K*(K-1)+2*cb(x)*(1-cb(x)) - floor(x)*(floor(x) -1) -2*cb(x)*floor(x))
    toRet.lt <- ((K-1)*(1 - x * (jf(x) -1) * (2 - jf(x) * x))/(K-1 + x * (jf(x) -1) *(2 - jf(x)* x)))
    toRet <- rep(NA, length(x))
    toRet[x >= 1] <- toRet.gt[x >= 1]
    toRet[x < 1] <- toRet.lt[x < 1]
    return(toRet)
}

readInPrior <- function(sampGens, n, nos = FALSE, gen2 = FALSE){

    trueMu <- .00001
    numrows <- 3000000

    stat.sim.full <- read.table(paste("../../raw/stats/n",n,"_prior.txt", sep = ""), 
    		  sep = " ", header = FALSE, nrows = numrows) 

    stat.sim.full <- tbl_df(stat.sim.full)
    colnamesToAdd <-  c(apply(expand.grid(c("diffInHets", "FST", "GSTp",
                              "thetaEst", "percDR", "numShared"),
                               paste(sampGens)),
                       1, paste, collapse = "_"), "N", "s1", "m", "n", "i")
    names(stat.sim.full) <- colnamesToAdd
    par.sim <- stat.sim.full %>%
              select(N, s1, m) %>% mutate(theta = N*trueMu) %>% select(theta, s1, m)

    if(nos == FALSE){

    	   if(gen2 == TRUE){
	           svals <- read.table(paste0("../../raw/svals/svals_n",n,"_gen2.txt"),
           	      sep = " ", header = FALSE, nrows = numrows)
        	seqRow <- seq(10, 70, 5)
	        relGensIt <- foreach(i = 1:length(seqRow))%do%{
		        c(2, seqRow[i],seqRow[i] + 30 )
			}
		}else{
		   svals <- read.table(paste0("../../raw/svals/svals_n",n,"_gen5.txt"),
           	      sep = " ", header = FALSE, nrows = numrows)
		seqRow <- seq(10, 70, 5)
                relGensIt <- foreach(i = 1:length(seqRow))%do%{
	  	      c(5, seqRow[i],seqRow[i] + 30 )
		    }
      	    }

	    names(svals) <- unlist(lapply(relGensIt, function(x){paste0(x, collapse = "_")}))
	 }else{
	 	 svals = NA
	 }
    return(list(stat.sim.full, par.sim, svals))
    }



#Sample me takes a tbl with generations (i.e., 'gen30') and returns samples
# of size n under the name (i.e., 'samp30')
sampleMe <- function(single, n){

    Acols <- single %>% filter(pop == "A") %>% select(starts_with('gen')) 
    Bcols <- single %>% filter(pop == "B") %>% select(starts_with('gen')) 

    colsToAddA <- apply(Acols, 2, function(x){ rmultinom(1, n, x) })
    colsToAddB <- apply(Bcols, 2, function(x){ rmultinom(1, n, x) })

    colsToAdd <- rbind(colsToAddA, colsToAddB)
    colnames(colsToAdd) <- gsub("gen", "samp", colnames(colsToAdd))

    tbl_df(cbind( single, colsToAdd))
}


computeStats <- function(single, anamesToLoop){

    #For now, we can throw in whatever stats we think will be helpful and 
    # then pare down after

    allNames <- single %>%
        filter(pop == "A") %>%
            select(mut)

    WTinds <- grep("WT", allNames$mut)

    ARows <- single %>%
        filter(pop == "A") 

    BRows <- single %>%
        filter(pop == "B") 

    allHets <- foreach(colname = anamesToLoop, .combine = "c") %do% {

        As <- ARows %>% select(matches(paste("^", colname, "$", sep = "")))
        Bs <- BRows %>% select(matches(paste("^", colname, "$", sep = "")))

        An <- sum(As)
        Bn <- sum(Bs)

        As <- As/An
        Bs <- Bs/Bn

#        DRinds <- setdiff(1:nrow(As), c(aWTInd, bWTInd))
        DRinds <- setdiff(1:nrow(As), WTinds)

        freqs <- (As + Bs)/2
#        DRfreqs <- freqs[DRinds,] * (An +Bn)

        countFreqs <- As[DRinds,] *An + Bs[DRinds,]* Bn
        DRfreqs <- sum(countFreqs)/(An + Bn)


        Het.A <- 1 - (sum(As[-c(WTinds),]^2) + sum(As[c(WTinds),])^2)
        Het.B <- 1 - (sum(Bs[-c(WTinds),]^2) + sum(Bs[c(WTinds),])^2)
        Het.AB <- 1 - (
            sum(As[-c(WTinds),] * Bs[-c(WTinds),]) +
            sum(As[c(WTinds),]) * sum(Bs[c(WTinds),]) )

        Het.AB.2 <- 1 - sum(As[DRinds, ] * Bs[DRinds,]) -
            sum(As[c(WTinds),] * Bs[c(WTinds),])
       
        fst1 <- (Het.AB - Het.A)/Het.AB
        fst2 <- (Het.AB - Het.B)/Het.AB
        
        #migration stats:
        diffInHets <- abs(Het.A - Het.B)
        FST <- (2*Het.AB - Het.A - Het.B)/(Het.A + Het.B + 2*Het.AB)
        if(is.nan(FST)){ FST <- 0 }
        GSTp <- FST/computeMax(max(As + Bs)/2, 2)
        if(is.nan(GSTp)){ GSTp <- 0 }

        if(length(countFreqs) == 0){ thetaEst <- 0
                                 }else{
                                     thetaEst <- optimal.theta(countFreqs)
                                 }

        #What is the percent WT?
        percDR <- DRfreqs 

        #Number of shared haplotypes
        numShared <- sum(Bs[DRinds,] > 0 & As[DRinds,] > 0)

        c(diffInHets, FST, GSTp, thetaEst, percDR, numShared)
        
    }

    sampNums <- as.numeric(gsub("samp|gen", "", anamesToLoop))
    Hets <- as.data.frame(t(allHets))
    colnames(Hets) <-  apply(expand.grid(c("diffInHets", "FST", "GSTp", "thetaEst",
                                           "percDR", "numShared"), paste(sampNums)),
                             1, paste, collapse = "_")

    return(Hets)

}


# subSetDat pares down the data to the correct time points 
subSetDat <- function(tblToSel, timepoints, varsToExc = "none"){

    Tregexp <- paste("_(", paste(timepoints, collapse = "|"), ")$", sep = "")
    rightTimes <- tblToSel %>% select( matches(Tregexp))

    excregexp <- paste("(", paste(varsToExc, collapse = "|"), ")_", sep = "")
    rightCols <-  rightTimes %>% select(-matches(excregexp))

    return(rightCols)
}


#I wrote a twostep minimization function that first coarsely searches
#a grid and then uses R's optimize once it has narrowed the search space
bespokeMinimization <- function(freqs, sampGens){

    if(sum(as.numeric(freqs)) == 0){ return(0) }
    coarse <- c(seq(.01, 1, by = .01), seq(2, 20, by = 1))
    coarsefits <- sAtGensVect(coarse, freqs, relGens = sampGens)
    fitval <- coarse[min(which(coarsefits == min(coarsefits)))]
    return(optimize(sAtGensVect, c(fitval * .5, fitval * 2),
                                freqs, relGens = sampGens)$minimum)

}


sAtGensVect <- function(sval, obs, relGens = relGens){

    assumedN <- 10^5
    vectorizedFreqs <- exp(as.matrix(sval) %*% relGens) / (exp(as.matrix(sval) %*% relGens) + matrix(rep(2*assumedN * sval, length(relGens)), nrow = length(sval)))

    vectorizedFreqs[which(is.nan(vectorizedFreqs))] <- 1
    obsFreqs <- matrix(rep(as.numeric(obs), length(sval)), nrow = length(sval), byrow = TRUE)

return(apply((vectorizedFreqs - obsFreqs)^2, 1, sum))

}


