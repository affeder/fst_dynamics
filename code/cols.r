#Set graphical values


baseCols <- c(rgb(51, 102, 170, maxColorValue = 255), 
              rgb(238, 51, 51, maxColorValue = 255), 
              rgb(153, 34, 136, maxColorValue = 255))


baseCols <- rev(c(rgb(51, 102, 170, maxColorValue = 255),
              rgb(17, 170, 153, maxColorValue = 255), 
              rgb(255, 238, 51, maxColorValue = 255), 
              rgb(238, 51, 51, maxColorValue = 255), 
              rgb(153, 34, 136, maxColorValue = 255)))

## baseCols <- c(rgb(0, 114, 178, maxColorValue = 255),
##               rgb(86, 180, 233, maxColorValue = 255), 
##               rgb(204, 121, 167, maxColorValue = 255), 
##               rgb(0, 158, 115, maxColorValue = 255), 
##               rgb(213, 94, 0, maxColorValue = 255),
##               rgb(230, 159, 0, maxColorValue = 255),
##               rgb(240, 228, 66, maxColorValue = 255))

#plot(1:5, col = baseCols, pch = 16, cex = 10)

redCol <- brewer.pal(3, "Set1")[1]
blueCol <-  brewer.pal(3, "Set1")[2]
yelCol <-  rgb(253, 228, 161, 1, max = 255)
purpleCol <-  brewer.pal(3, "Set1")[5]

makeColVals <- function(mutNames, buff = 2){

    #Things in A are going to be red -> black
    mutNames <- unique(mutNames)$mut

    redCol <- brewer.pal(3, "Set1")[1]
    blueCol <-  brewer.pal(3, "Set1")[2]

    aPal <- colorRampPalette(c("white", redCol))
    bPal <- colorRampPalette(c("white", blueCol))
    cPal <- colorRampPalette(brewer.pal(8, "Set1")[1:6][c(4, 2, 3, 6)])

    aInds <- grep("[0-9]+-a", mutNames)
    bInds <- grep("[0-9]+-b", mutNames)
    cInds <- grep("[0-9]+-c", mutNames)
    
    toRet <- rep(NA, length(mutNames))

    offs <- buff
    
    aCols <- aPal(length(aInds)+offs)[-c(1:offs)]
    bCols <- bPal(length(bInds)+offs)[-c(1:offs)]

    toRet[aInds] <- sample(aCols, length(aCols))
    toRet[bInds] <- sample(bCols, length(bCols))
    toRet[cInds] <- cPal(length(cInds)+offs)[-c(1:offs)]

    toRet[grep("WT-c", mutNames)] <- "grey"
    toRet[grep("WT", mutNames)] <- "grey"
    toRet[grep("WT-a", mutNames)] <- "grey"
    toRet[grep("WT-b", mutNames)] <- "grey"

    names(toRet) <- mutNames
    return(toRet)
}
