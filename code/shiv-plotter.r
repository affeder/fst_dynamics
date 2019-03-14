source('abc-functions-2.r')
library(abc)
library(doParallel)
library(emdbook)
library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(cowplot)

trueMu <- .00001

n <- 30
relGens <- c(1:4, seq(5, 100, by = 5))
sampGens <- c(5, 20, 50, 100)

#set the colors
maf <- 1
allmuts <- unique(tbl_df(read.table(paste0("../output/shivdata-maf", maf,".txt", 
                                           sep = ""), 
                             header = TRUE, stringsAsFactors = FALSE))$mut)


colFunc <- colorRampPalette(brewer.pal(9, "Set1")[c(1:6, 8)])
cols.mullers <- c("grey", (colFunc(length(allmuts) - 1)))
names(cols.mullers) <- allmuts

colref <- names(cols.mullers)
coladj <- names(cols.mullers)

## plot(1:length(cols.mullers), col = cols.mullers, pch = 16)
## text(1:length(cols.mullers), labels = names(cols.mullers), col = cols.mullers, pch = 16)

nameswitch <- matrix(c("M184V-N255N", "I47N-L149L-M184V",
                       "L149L-M184V", "M184I-K223K", #"M184V-E204E-L205L-K249K",
                      "Y144Y-I178M-M184V", "S156*-M184I"), 
                     ncol = 2, byrow = TRUE)

for(i in 1:nrow(nameswitch)){
    fromval <-nameswitch[i,1]
    toval <-nameswitch[i,2]
    coladj[which(colref == fromval)] = 
        colref[which(colref == toval)]
    coladj[which(colref == toval)] = 
        colref[which(colref == fromval)]
}
names(cols.mullers) <- coladj


cols.mullers[which(names(cols.mullers) == "M184V")] <- baseCols[2]
cols.mullers[which(names(cols.mullers) == "M184V+N255N")] <- baseCols[1]
cols.mullers[which(names(cols.mullers) == "M184I")] <- baseCols[5]

outdists <- foreach(maf = c(1, 2, 5), .combine = 'rbind') %do%{

trueDat <- tbl_df(read.table(paste0("../output/shivdata-maf", maf,".txt", sep = ""), 
                             header = TRUE, stringsAsFactors = FALSE))

allT <- trueDat %>% mutate(pop = ifelse(loc == "PLASMA", "Plasma",
                         ifelse(loc == "LN", "Lymph Node",
                         ifelse(loc == "GUT", "Gut", NA)))) %>% 
                    mutate(mutname = mut) %>% select(-loc)

mutlevs <- c(unique(allT$mut)[-1], "WT")

####################################################################################################

d1 <- allT %>% mutate(mutname = factor(mutname, levels = mutlevs)) %>%
        gather(sampgen, freq, -mut, -pop, -mutname) %>% 
        mutate(sampgen = as.numeric(gsub("samp", "", sampgen))) %>%
        group_by(pop, sampgen) %>%
        arrange(pop, desc(mutname)) %>%
        mutate(freq = cumsum(freq)) %>%
        ungroup() %>%
        arrange(sampgen) 

d2 <- d1 %>% mutate(freq = 0) %>% arrange(desc(sampgen))
d <- rbind(d1, d2)

ypos <- d1 %>% arrange(pop, sampgen, desc(mutname)) %>% group_by(pop, sampgen) %>%
    mutate(lfreq = lead(freq, order_by = mutname)) %>% 
    mutate(lfreq = ifelse(is.na(lfreq), 0, lfreq)) %>% 
    mutate(diff =   freq - lfreq) %>% mutate(csum = cumsum(diff)) %>%
    mutate(ypos = freq - diff/2) %>%
        select(-freq, -lfreq, -diff, -csum) %>% ungroup()


xpos <- allT %>% gather(sampgen, freq, -mut, -pop, -mutname) %>% 
    mutate(sampgen = as.numeric(gsub("samp", "", sampgen))) %>%
    filter(freq > 0) %>% group_by(pop, mut) %>% mutate(ind = 1:n()) %>% 
    filter(ind == 1) %>% ungroup() %>% select(pop, mutname, sampgen) %>% 
    mutate(ftp = 1)


labpos <- inner_join(xpos, ypos, by = c("pop", "mutname", "sampgen")) %>% 
    mutate(hjust = ifelse(sampgen == 100, 1, 0)) %>%
    mutate(xpos = ifelse(sampgen == 100, sampgen -1.5, sampgen + 1.5))

    if(maf == 2){

        labpos <- labpos %>% filter(!(pop == "Lymph Node" & mutname == "L149L-M184V"))
        
    }


toPlot <- d 
mullers <- toPlot %>% 
#        mutate(mutname = factor(mutname, levels = o.mutnames)) %>%
        ggplot() +
        geom_polygon(mapping = aes(x = sampgen,  y = freq, fill = mutname) , alpha = 1) +
        theme_gen() + theme(legend.position = "none") +
#        theme(legend.position="top") + 
        guides(fill = guide_legend(ncol = 3, override.aes = list(size=1))) + 
        theme(legend.text = element_text(size=5)) + 
        labs( y= "Count", x = "Generation", fill = "") +
        scale_fill_manual(values = cols.mullers) + 
        facet_wrap(~as.factor(pop), ncol = 1) +
        scale_y_continuous(limits = c(0, 35)) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(), 
              axis.ticks.x=element_blank()) + 
        theme(strip.background = element_blank(), strip.text.x = element_blank()) + 
        geom_text(aes(x = 100, y= 34, label = pop), hjust = 1, 
                  data = toPlot %>% select(pop) %>% distinct()) + 
        theme(plot.margin = margin(0,.25,0,0, "cm")) + 
        geom_text(aes(x = xpos, y = ypos, label = mutname, 
                      hjust = hjust), 
                  data = labpos, size = 2, col = "white") + 
        geom_segment(aes(x = sampgen, y = ypos, xend = xpos, yend = ypos),
                  data = labpos, col = "white")


##########################################################################################

sampNames <- paste("samp", sampGens, sep = "")

plas.ln <- allT %>% filter(grepl("Plasma|Lymph", pop)) %>% 
    mutate(pop = str_replace_all(pop, setNames(c("A", "B"), c("Plasma", "Lymph Node"))))
plas.gut <- allT %>% filter(grepl("Plasma|Gut", pop)) %>% 
    mutate(pop = str_replace_all(pop, setNames(c("A", "B"), c("Plasma", "Gut"))))
gut.ln <- allT %>% filter(grepl("Gut|Lymph", pop)) %>% 
    mutate(pop = str_replace_all(pop, setNames(c("A", "B"), c("Gut", "Lymph Node"))))


t1 <- plas.ln %>% computeStats(sampNames) %>% select(starts_with('FST')) %>% 
    mutate(comps ="Plasma v LN")
t2 <- plas.gut %>% computeStats(sampNames) %>% select(starts_with('FST'))%>% 
    mutate(comps = "Plasma v Gut")
t3 <- gut.ln %>% computeStats(sampNames) %>% select(starts_with('FST'))%>% 
    mutate(comps = "Gut v LN")


cols.fst <- c(baseCols[2], baseCols[5], baseCols[1])
names(cols.fst) <- c("Plasma v LN", "Plasma v Gut", "Gut v LN")

fst_plot_dat <- bind_rows(t1, t2, t3) %>% gather(Generation, FST, -comps) %>% 
    mutate(Generation = as.numeric(gsub("FST_", "", Generation))) %>% 
    mutate(comps = factor(comps, levels = c("Plasma v LN", "Plasma v Gut", "Gut v LN")))


fstplot <- fst_plot_dat %>% 
    ggplot() + geom_line(aes(x = Generation, y = FST, group = comps, col = comps)) + 
    scale_color_manual(values = cols.fst) + 
    labs(lty = "", y = expression("F"['ST']), x = "Time (Generations)", col = "")  +
    theme_gen() + theme(legend.position=c(.25, .85)) + 
    theme(plot.margin = margin(0,.25,0,0, "cm"))


g1 <- plot_grid(mullers, fstplot, ncol = 1, align = "v", rel_heights = c(1, 2/5)) 

pdf(paste0("../graphs/mullers_and_fst_",maf,".pdf"), width = 3, height = 7.5)
print(g1)
dev.off()


    allmuttypes <- unique(trueDat$mut)
    ## tolval1 <- .00025
    ## tolval2 <- .00025

    #generally, I use a MUCH lower tolerance and a MUCH larger sample size
    # for illustrative purposes, here I use tol = 0.1 so it's tractable to run locally
    tolval = 0.1
    #The real function here is called realDat (it's in read-in-big-files-perc.r), and 
    # it reads in two much bigger priors, but otherwise uses the same summaries
pvg <- realDat.no_r2(trueDat %>% filter(loc != "LN"), tol1 = tolval, tol2 = tolval)
pvl <- realDat.no_r2(trueDat %>% filter(loc != "GUT"), tol1 = tolval, tol2 = tolval)
#pvg <- realDat(trueDat %>% filter(loc != "LN"), tol1 = tolval, tol2 = tolval)
#pvl <- realDat(trueDat %>% filter(loc != "GUT"), tol1 = tolval, tol2 = tolval)
    
parorder <- c("Population Mutation Rate", 
              "Selection strength",
              "Migration probability")

cols <- c(baseCols[2], baseCols[5])
    abcProcBigRet1 <- pvl
    abcProcBigRet2 <- pvg
    yval.s <- -5
    yval.s.2 <- -10
    yval <- yval.s
yval2 <- yval.s.2

histPlot1 <- bind_rows(abcProcBigRet1[[1]] %>% filter(grepl("theta|s1", var)), 
              abcProcBigRet1[[2]] %>% filter(grepl("m", var))) %>% 
              mutate(comp = "Plasma v LN")
    histPlot2 <- bind_rows(abcProcBigRet2[[1]] %>% filter(grepl("theta|s1", var)), 
              abcProcBigRet2[[2]] %>% filter(grepl("m", var))) %>%
              mutate(comp = "Plasma v Gut")
    histPlot <- rbind(histPlot1, histPlot2)
    histPlot <- histPlot %>% mutate(var = ifelse(var == "s1", "s", var))
    rangePlots <- histPlot %>% group_by(var, comp) %>%
        summarize(l95 = quantile(val, .025),
                  l50 = quantile(val, .25),
                  m50 = quantile(val, .5),
                  h50 = quantile(val, .75),
                  h95 = quantile(val, .975)) %>% 
        mutate(yval = ifelse(comp == "Plasma v Gut", yval, yval2)) %>%
        mutate(yval = ifelse(var == "m" & comp == "Plasma v Gut", yval.s, yval)) %>%
            mutate(yval = ifelse(var == "m" & comp == "Plasma v LN", yval.s.2, yval))
    rangePlots <- rangePlots %>% ungroup() %>%
        mutate(var = ifelse(var == "s1", "s", var))
rangevals <- tbl_df(apply(par.sim, 2, range) )  %>% 
   gather(var,val) %>%         mutate(var = ifelse(var == "s1", "s", var))

        
Nval <- (histPlot %>% filter(var == "theta") %>% summarize(Nval = median(val))/trueMu)$Nval

histPlot <- histPlot %>% 
    mutate(var = ifelse(var == "theta", "Population Mutation Rate",
                 ifelse(var == "s", "Selection strength", 
                 ifelse(var == "m", "Migration probability", var)))) %>%
    mutate(var = factor(var, levels = parorder))
rangevals <- rangevals %>% 
    mutate(var = ifelse(var == "theta", "Population Mutation Rate",
                 ifelse(var == "s", "Selection strength", 
                 ifelse(var == "m", "Migration probability", var)))) %>%
    mutate(var = factor(var, levels = parorder))
rangePlots <- rangePlots %>%
    mutate(var = ifelse(var == "theta", "Population Mutation Rate",
                 ifelse(var == "s", "Selection strength", 
                 ifelse(var == "m", "Migration probability", var)))) %>%
    mutate(var = factor(var, levels = parorder))

hists <- histPlot %>% 
    ggplot() +
        geom_histogram(mapping = aes(x = val, fill = comp, group = comp),
                       color = "black",
                       position = "identity", alpha = .5) + 
        facet_wrap(~var, scales = "free_x", ncol = 3) + 
        scale_x_log10(labels=fancy_scientific, 
                      breaks = c(.00001, .0001, .001, .01, .1, 1, 10)) +
        theme_gen() + 
        geom_segment(mapping = aes(x = l95, xend = h95, y = yval, yend = yval, 
                          col = comp), data = rangePlots, size = 1) +
        geom_segment(mapping = aes(x = l50, xend = h50, y = yval, yend = yval, 
                          col = comp), data = rangePlots, size = 2) + 
        geom_point(mapping = aes(x = val, y = 0), alpha = 0, data = rangevals)+
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) + 
        labs(y = "Count", fill = element_blank(), 
             x = expression(paste('      ',theta,"                                                   s                                                      m  ", sep = ""))) +
        guides(col=FALSE) + 
        theme(legend.position = c(0.9, 0.85))
#        theme(legend.position = "top")



pdf(paste0("../graphs/hists_maf",maf,"_tol1",tolval,"_tol2",tolval,".pdf"), width = 7.47, height = 3)
print(hists)
dev.off()

return(rangePlots %>% mutate(maf = maf))
}


Ns <- outdists %>% group_by(comp) %>% filter(grepl("Population",var)) %>% 
    mutate(N = m50/trueMu)  %>% select(comp, maf, N)

Migs <- left_join(outdists %>% filter(grepl("Migration", var)), Ns) %>% 
    mutate(l95 = l95*N, m50 = m50*N, h95 = h95*N) %>% select(-N) %>%
        mutate(var = "M")

tab1 <-   bind_rows(outdists, Migs) %>% mutate(var = paste(var)) %>%
    mutate(var = ifelse(grepl('Migration', var), 'm', ifelse(
                        grepl('Population', var), 'theta', ifelse(
                        grepl('Selection', var), 's', var)))) %>%
    mutate(ci = paste0("(", signif(l95, 3), ",", signif(h95, 3), ")")) %>%
    mutate(m50 = signif(m50, 3)) %>%
    select(comp,maf, var, m50, ci  )  %>%
    gather(x, par, -comp, -maf, -var) %>% 
    filter(maf == 2) %>% select(-maf) %>%
    spread(var, par) %>%
        arrange(comp,  desc(x))  %>% select(-x) %>%
    select(comp, m, M, s, theta)


tab2 <-   bind_rows(outdists, Migs) %>% mutate(var = paste(var)) %>%
    mutate(var = ifelse(grepl('Migration', var), 'm', ifelse(
                        grepl('Population', var), 'theta', ifelse(
                        grepl('Selection', var), 's', var)))) %>%
    mutate(ci = paste0("(", signif(l95, 3), ",", signif(h95, 3), ")")) %>%
    mutate(m50 = signif(m50, 3)) %>%
    select(comp,maf, var, m50, ci  )  %>%
    gather(x, par, -comp, -maf, -var) %>% 
    spread(var, par) %>%
    select(comp, maf, m, M, s, theta, x)  %>%
    arrange(maf, comp, desc(x)) %>% select(-x)



print(xtable(tab1), include.rownames=FALSE, file = "../output/tab1.txt")
print(xtable(tab1), include.rownames=FALSE)

print(xtable(tab2), include.rownames=FALSE, file = "../output/tab2.txt")
print(xtable(tab2), include.rownames=FALSE)

