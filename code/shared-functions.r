#Shared functions

heatMapHelper <- function(x, qs, testvals, varval = 'm'){

    x.spread <- reformatForHeatMap(x, qs, varval)

    intervals <- foreach(test = testvals, .combine = "rbind") %do% {
        x.spread %>% mutate(testm = test, inInt = low <= test & high >= test ) 
    }

    return(intervals)
}

reformatForHeatMap <- function(x, qs, varval = 'm'){

    qlow <- (1-qs)/2
    qhigh <- 1-(1-qs)/2

    xnames <- names(x)

    x.spread <- x %>% filter(var == varval) %>% 
        filter(near(q,qlow) | near(q,qhigh)) %>%
        spread(q, val) 

    xnames <- c(xnames[-which(xnames == "val" | xnames == "q")], "low", "high")
    names(x.spread) <- xnames

    return(x.spread)

}


fancy_scientific <- function(l) {
#Modified from: https://groups.google.com/forum/#!topic/ggplot2/a_xhMoQyxZ4
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e", "10^", l)
     # return this as an expression
     l <- gsub("0e\\+00","0",l)
     l <- gsub("\\+", "", l)
     parse(text=l)
 }
