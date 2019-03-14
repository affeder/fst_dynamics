

set.seed(12412663)
trueMu <- .00001
ms <- c(.000001, .00001, .0001, .001, .01, .1)
thetas <- c(1)
ss <- c(0)
numtrials <- 1
runs <- expand.grid(thetas, ss, ms, 1:numtrials)
thetaToRun <- runs[,1]
sToRun <- runs[,2]
mToRun <- runs[,3]
trialsToRun <- 1:nrow(runs)#runs[,4]
len <- 1/trueMu*5
sampGens <- floor(seq(1, len, length.out = 500))

#generate the trials
cl <- makeCluster(4)
registerDoParallel(4)
trajs <- foreach( theta = thetaToRun,
                      s = sToRun,
                      m = mToRun,
                      i = trialsToRun,
                      .combine='rbind',
                      .packages = c('foreach', 'dplyr')) %dopar% {
       dat <- subSim(floor(theta/trueMu), floor(theta/trueMu), trueMu,
              s, s, m, len, relgens =sampGens, recenterCond = NA) %>% 
                  mutate(m = m, i = i)
   }
stopCluster(cl)
trajs <- tbl_df(trajs)


write.table(trajs, "../d_output/longtrajs-12412663.txt", quote = FALSE, 
            col.names = TRUE,   row.names = FALSE)



len <- 1/trueMu*5
sampGens <- floor(seq(1, len, length.out = 500))
    
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

wid <- 10000

set.seed(124909)
colVals <- makeColVals(mutNames, 500)

p <- merged %>% 
    filter(freq != 0 | is.na(freq)) %>% ggplot() +
        geom_bar(mapping = aes(x = gen, y = freq, fill = as.factor(mut)), 
                 stat = "identity",  show.legend = FALSE, width = wid) +
        geom_line(mapping = aes(x = gen, y = FST), size = 1) +
        facet_grid(pop ~ Nm , scales = "free", switch = "y")+
        scale_fill_manual(values = colVals) + 
        theme(axis.title.y=element_blank()) +
        scale_y_continuous(position = "left") + 
        geom_hline(aes(yintercept  = FST), data = expFST, col = "black",
                   lwd = .25, lty = "dashed", size = 2) + 
        scale_color_manual(values = decols) +
        labs(x =expression(paste("Generations (",N['e'],")", sep = "")),
                   fill = "Mutation", 
             y = expression(paste('      F'['ST'],"(A,B)           Freq. in Pop. B       Freq. in Pop. A ", sep = ""))) + 
        scale_x_continuous(breaks = N*c(0:5), labels = paste(c("",1, "",3,"",5))) + 
        theme_gen()  + 
        theme(strip.background =  element_blank(),
                   strip.text.y = element_blank(),
              strip.text.x = element_text(size = 11))


pdf("../d_graphs/examples_s0_3.pdf", width = 7.87, height = 4.72 , pointsize = 12)
print(p)
dev.off()




