#justsims.r is going to have functions JUST to run the simulation
#abc-functions.r is going to have all of the functions that deal with
#  computing statistics on the simulated runs, and helper plotting functions

library(foreach)
library(tidyverse)
#library(grid)
#library(gridExtra)

mutation <- function(variants, mu){

    #of the remaining WT strains, a proportion mu should become new types

    #How many WT individuals are there now?
    WTind <- grep("WT", rownames(variants))

    if(length(WTind) > 0){
        newMuts <- variants[WTind,]*mu
        names(newMuts) <- paste("NEW", names(newMuts), sep = "-")

                                        #Add the new mutations
        variants <- rbind(variants,  as.matrix(newMuts))
                                        #and subtract them from the WT frequency
        variants[WTind, ] <- variants[WTind,] - newMuts

    }

    #return the results
    return(variants)
    
}


selection <- function(variants, s.home, s.migrant, compname){

    #this function should return the deterministic proportions of offspring
    #the weighting function will be the multiplier against which the frequencies
    #are multiplied 

    #Any migrants get the migrant selective effect
    weighting <- rep(1+s.migrant, length(variants))

    #WT is neutral
    weighting[grep("WT", rownames(variants))] <- 1

    #mutants originating in this compartment get a home selective effect
    weighting[grep(paste("[0-9]+-", compname, sep = ""), rownames(variants))] <- 1 + s.home
    
    weighted <- variants*weighting
    weighted <- weighted/sum(weighted) #normalize to frequencies
    names(weighted) <- names(variants)
    
    return(weighted)
    #We'll deal with forgetting variants in the main function
}

migration <- function(recipient.pop, donor.pop, m){

    #recipient.pop and donor.pop need not have the same muts, so 
    # we will need to be a little more careful

    #merge the two (weighted) populations based on rownames
    preadd <- merge(recipient.pop*(1-m), donor.pop*m,
                    by = "row.names", all = TRUE)

    #We'll now reformat the row.names (first col) back into names
    postadd <- as.matrix(apply(preadd[,2:3], 1, sum, na.rm = TRUE))
    rownames(postadd) <- preadd[,1]

    return(postadd)
}
    

multinomsamp <- function(freqs, N){
    #This function is going to do resampling, and then
    # for each new type it draws, it will create a new haplotype id

    aorb <- NULL

    v.freqs <- as.vector(freqs)
    names(v.freqs) <- rownames(freqs)

    #First, just sample
    newVars <- rmultinom(1, N, v.freqs)

    #If there are new mutations, label them with currmut
    newMutsToAdd <- 0
    newInd <- grep("NEW", rownames(newVars))

    for(addInd in newInd){

        newMutsToAdd <- newVars[addInd]
        newMutType <- rownames(newVars)[addInd]

        if(grepl("-a", newMutType)){ aorb <- "a" }
        if(grepl("-b", newMutType)){ aorb <- "b" }

        if(newMutsToAdd > 0){ 
            appendMat <- c()
            while(newMutsToAdd > 0){
                appendMat <- append(appendMat,
                    setNames(1, paste(c(nextMut(), aorb), collapse = "-")))
                newMutsToAdd <- newMutsToAdd - 1
            }
            newVars <- rbind(newVars, as.matrix(appendMat))
        }
    }
    if(length(newInd) > 0){
    #now, delete the "NEW" placeholders
        newVars <- as.matrix(newVars[-newInd, ])
    }
    return(newVars)
}


nextMut <- function(){
    #<<- is a global operator
    currMut <<- currMut + 1
    return(currMut)
}



IBD.breakdown <- function(allpairs, aInds, bInds, aWTInd, bWTInd){

#same compartment, same state
    p.CM1.SM1 <- sum(diag(allpairs)[-c(aWTInd, bWTInd)])
    if(length(c(aWTInd, bWTInd)) == 0 ){
        p.CM1.SM1 <- sum(diag(allpairs))
    }

#same compartment, different states
    p.CM1.SM0 <-
        sum(allpairs[aInds,aInds]) +  sum(allpairs[bInds,bInds]) +
            sum(diag(allpairs)[c(aWTInd, bWTInd)])

#we also want to subtract out the diagonals here
    toSubtract <- 0
    if(length(bInds) == 1){
        toSubtract <- toSubtract + allpairs[bInds, bInds]
    }else{
        toSubtract <- toSubtract + sum(diag(allpairs[bInds, bInds]))
    }
    if(length(aInds) == 1){
        toSubtract <- toSubtract + allpairs[aInds, aInds]
    }else{
        toSubtract <- toSubtract + sum(diag(allpairs[aInds, aInds]))
    }

    p.CM1.SM0 <- p.CM1.SM0 - toSubtract

#different compartments
    p.CM0.SM0 <-
        sum(allpairs[aInds,bInds]) + sum(allpairs[bInds,aInds])

    props <- c(p.CM1.SM1,  p.CM1.SM0,  p.CM0.SM0)

    return(props)

}


oneGeneration <- function(P.a.variants, P.b.variants, N.a, N.b, mu, s1, s2, m){

         #Ok, deterministic infinite number of offspring
        # Note: these are now frequencies (not counts)
    postsel.freqs.a <- selection(P.a.variants, s1, s2, "a")
    postsel.freqs.b <- selection(P.b.variants, s1, s2, "b")

        #Next, migration. We can simply take m variants from popa and move them
        #to popb and vice versa. (symmetric migration)
    postmig.freqs.a <- migration(postsel.freqs.a, postsel.freqs.b, m)
    postmig.freqs.b <- migration(postsel.freqs.b, postsel.freqs.a, m)
        
        #Third: mutation (a fixed fraction of WT turn into new haplos)
    postmut.freqs.a <- mutation(postmig.freqs.a, mu)
    postmut.freqs.b <- mutation(postmig.freqs.b, mu)

        #Now, we do our sampling
    P.a.variants <- multinomsamp(postmut.freqs.a, N.a)
    P.b.variants <- multinomsamp(postmut.freqs.b, N.b)

        #If there are categories that are 0, delete them
    if(sum(P.a.variants == 0) > 0){
        P.a.variants <- as.matrix(P.a.variants[-which(P.a.variants == 0),])
    }
    if(sum(P.b.variants == 0) > 0){
        P.b.variants <- as.matrix(P.b.variants[-which(P.b.variants == 0),])
    }

    return(list(P.a.variants = P.a.variants, P.b.variants = P.b.variants))

}

recenterFunc <- function(newFreqs){
    helper <- function(x){
        matches <- grep("WT", rownames(x), invert = TRUE)
        if(length(matches) > 0){
            return(sum(x[matches,] > (.005) * N.a) > 0)
        }
        return(FALSE)
    }
    return( sum(unlist(lapply(newFreqs, helper))) > 0)
}

subSim <- function(N.a, N.b, mu, s1, s2, m, len, relgens = NA, recenterCond = NA){

     ## N.a <- 10^5
     ## N.b <- 10^5
     ## mu <- .00001
     ## s1 <- 1
     ## s2 <- 1
     ## m <- .0001
     ## len  <- 500

    #Set up data structure
    P.a.variants <- as.matrix(setNames(c(N.a), "WT-a"))
    P.b.variants <- as.matrix(setNames(c(N.b), "WT-b"))
    newFreqs <- list(P.a.variants = P.a.variants, P.b.variants = P.b.variants)
    currMut <<- 0
    
    #Ok, if we have a recenterCondition, run into the condition is met
    runsUntilRecenter <- 0
    if(!is.na(recenterCond)){
        while(!recenterFunc(newFreqs)){
            newFreqs <- oneGeneration(newFreqs$P.a.variants, newFreqs$P.b.variants,
                                      N.a, N.b, mu, s1, s2, m)
            runsUntilRecenter <- runsUntilRecenter + 1
        }
    }

    stotmp <- list()
    
    #Runs the simulation
    for(i in 1:(len+1)){
        
        newFreqs <- oneGeneration(newFreqs$P.a.variants, newFreqs$P.b.variants,
                                  N.a, N.b, mu, s1, s2, m)

        if(sum(grepl(paste("^", i, "$", sep = ""), relgens)) > 0 |
           sum(is.na(relgens)) > 0){

            preAdd <- merge(newFreqs$P.a.variants, newFreqs$P.b.variants,
                            by = "row.names", all = TRUE)
            preAdd[is.na(preAdd[,2]),2] <- 0
            preAdd[is.na(preAdd[,3]),3] <- 0

            preAdd[,2] <- preAdd[,2]/sum(preAdd[,2])
            preAdd[,3] <- preAdd[,3]/sum(preAdd[,3])

            stotmp[[paste(i)]] <- preAdd
        }

    }

    #Get all unique names
    allnames <- unique(unlist(lapply(stotmp, function(x){return(x[,1])})))

    allnames <- sort(allnames)
    #Create matrices to keep track of their time of occurrence 
    popa <- matrix(data = 0, nrow = length(allnames), ncol = length(stotmp))
    popb <- matrix(data = 0, nrow = length(allnames), ncol = length(stotmp))

    for(i in names(stotmp)){
        tmpnames <- stotmp[[i]][,1]
        for(j in 1:length(tmpnames)){
            colInd <- which(i == names(stotmp))
            popa[which(allnames == tmpnames[j]), colInd] <- stotmp[[i]][j,2]
            popb[which(allnames == tmpnames[j]), colInd] <- stotmp[[i]][j,3]
        }
    }

    if(sum(is.na(relgens)) > 0){
        relgens <- 1:ncol(popa)
    }
    
    colnames(popa) <- paste("gen", relgens, sep = "")
    colnames(popb) <- paste("gen", relgens, sep = "")

    #Put the different variants together to return in a long format
    toReturn <- rbind(data.frame(pop = "A", mut = allnames, popa, runsUntilRecenter),
                 data.frame(pop = "B", mut = allnames, popb, runsUntilRecenter))

    return(toReturn)
    
}




