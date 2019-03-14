
examp_plot_s <- function(trajs){

    genOrSamp = "gen"
    N <- 10^5
    mutNames <- trajs %>% select(mut)

    orderedFact <-
        paste(unique(mutNames %>% arrange(grepl("[0-9]+-[a|b]", mut)))$mut)

    datToPlot <- trajs %>%
        select(pop, mut, m, starts_with(genOrSamp)) %>%
            gather(gen, freq, -pop, -mut, -m)  %>%
        mutate(gen = as.numeric(gsub(genOrSamp, "", gen))) %>%
        mutate(Nm = paste0("Nm = ",N*m)) %>% select(-m) %>%
            mutate(pop = paste0(" Freq. in pop ", pop))

    expFST <- unique(trajs %>% select(m)) %>% mutate(Nm = paste0("Nm = ",N*m)) %>%
        mutate(FST = 1/(1 + 4*N*m)) %>% mutate(pop = "FST(A,B)")
  
    stats <- trajs %>% group_by(m) %>%
        do(computeStats(., paste("gen", sampGens, sep = "")))  %>% 
        select(starts_with("FST"))

    toAdd <- stats %>% gather(key = gen, value = FST, -m) %>% ungroup() %>%
        mutate(gen = as.numeric(gsub("FST_", "", gen))) %>%
        mutate(Nm = paste0("Nm = ",N*m)) %>% select(-m) %>%
        mutate(pop = "FST(A,B)") %>% mutate(mut = NA)

    toMerge1 <- datToPlot %>% mutate(mut = factor(mut, levels = orderedFact)) %>%
        mutate(FST = NA)
    toMerge2 <- toAdd %>% mutate(freq = NA)
    merged <- rbind(toMerge1, toMerge2)

    wid <- 2

    set.seed(1231298)
    colVals <- makeColVals(mutNames)

    p <- merged %>% 
        filter(freq != 0 | is.na(freq)) %>% ggplot() +
        geom_bar(mapping = aes(x = gen, y = freq, fill = as.factor(mut)), 
                 stat = "identity",  show.legend = FALSE, width = wid) +
        geom_line(mapping = aes(x = gen, y = FST), size = 1) +
        facet_grid(pop ~ Nm , scales = "free", switch = "y")+
        scale_fill_manual(values = colVals) + 
        theme(axis.title.y=element_blank()) +
        scale_y_continuous(position = "left") + 
        scale_color_manual(values = decols) +
        labs(x = "Generations", 
             fill = "Mutation", 
             y = expression(paste('      F'['ST'],"(A,B)           Freq. in Pop. B       Freq. in Pop. A ", sep = ""))) + 
        scale_x_continuous(breaks = seq(0, 100, by = 25), 
                           labels = c("", 25, 50, 75, 100)) + 
        theme_gen()  + 
        theme(strip.background =  element_blank(),
                   strip.text.y = element_blank(),
              strip.text.x = element_text(size = 11))

    return(p)
}




trueMu <- .00001
ms <- c(.000001, .00001, .0001, .001, .01, .1)
ss <- c(1)
numtrials <- 1
len <- 100
sampGens <- 1:100

#Theta = 0.1 (to demonstrate migrant-driven sweeps)
thetas <- c(.1)

runs <- expand.grid(thetas, ss, ms, 1:numtrials)
thetaToRun <- runs[,1]
sToRun <- runs[,2]
mToRun <- runs[,3]
trialsToRun <- 1:nrow(runs)

set.seed(65357)
trajs <- foreach( theta = thetaToRun,
                      s = sToRun,
                      m = mToRun,
                      i = trialsToRun,
                      .combine='rbind',
                      .packages = c('foreach', 'dplyr')) %do% {
       dat <- subSim(floor(theta/trueMu), floor(theta/trueMu), trueMu,
              s, s, m, len, relgens =sampGens, recenterCond = NA) %>% 
                  mutate(m = m, i = i)
    }
trajs <- tbl_df(trajs)


pdf("../d_graphs/examples_migrant.pdf", width = 7.87, height = 4.72)
print(examp_plot_s(trajs))
dev.off()




#Theta = 0.5 (to demonstrate local-driven sweeps)
thetas <- c(.5)

runs <- expand.grid(thetas, ss, ms, 1:numtrials)
thetaToRun <- runs[,1]
sToRun <- runs[,2]
mToRun <- runs[,3]
trialsToRun <- 1:nrow(runs)

set.seed(3633)
trajs <- foreach( theta = thetaToRun,
                      s = sToRun,
                      m = mToRun,
                      i = trialsToRun,
                      .combine='rbind',
                      .packages = c('foreach', 'dplyr')) %do% {
       dat <- subSim(floor(theta/trueMu), floor(theta/trueMu), trueMu,
              s, s, m, len, relgens =sampGens, recenterCond = NA) %>% 
                  mutate(m = m, i = i)
    }
trajs <- tbl_df(trajs)


pdf("../d_graphs/examples_local.pdf", width = 7.87, height = 4.72)
print(examp_plot_s(trajs))
dev.off()

