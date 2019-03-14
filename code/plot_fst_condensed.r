trueMu <- .00001
ms <- c(.0001, .001, .01)
thetas <- c(1)
ss <- c(.3, 3)
numtrials <- 50
runs <- expand.grid(thetas, ss, ms, 1:numtrials)
thetaToRun <- runs[,1]
sToRun <- runs[,2]
mToRun <- runs[,3]
trialsToRun <- 1:nrow(runs)
len <- 100
sampGens <- 1:100
#generate the trials
cl <- makeCluster(4)
registerDoParallel(4)
fsts <- foreach( theta = thetaToRun,
                      s = sToRun,
                      m = mToRun,
                      i = trialsToRun,
                      .combine='rbind',
                      .packages = c('foreach', 'dplyr')) %dopar% {
       tbl_df(subSim(floor(theta/trueMu), floor(theta/trueMu), trueMu,
       s, s, m, len, relgens = 1:100, recenterCond = NA))  %>%
       mutate(muttype = ifelse(grepl("[0-9]+-a", mut), "1-a",
                          ifelse(grepl("[0-9]+-b", mut), "2-b", 
                          ifelse(grepl("WT-b", mut), "WT-b", 
                          ifelse(grepl("WT-a", mut), "WT-a", NA))))) %>%
        mutate(mut = muttype) %>% select(-muttype)  %>%
        group_by(pop, mut) %>% summarize_all(funs('sum')) %>% ungroup() %>%
        do(computeStats(., paste("gen", sampGens, sep = ""))) %>% 
        select(starts_with("FST")) %>%
        mutate(N = floor(theta/trueMu), s1 = s, m = m, i = i)
}          
stopCluster(cl)
fsts <- tbl_df(fsts)

write.table(fsts, "../output/fstvals.txt", quote = FALSE, row.names = FALSE,
            col.names = TRUE)

fsts <- tbl_df(read.table("../output/fstvals.txt", header = TRUE))

helper <- function(x){return(tbl_df(cbind(x, gen = -5:105))) }
eqFSTs <- 1 / (1 + 4*ms/trueMu)
eqFSTs <- tbl_df(cbind(ms/trueMu, eqFSTs))
names(eqFSTs) <- c("m", "freq")
eqFSTs <- eqFSTs %>% mutate(m = ifelse(m > 1, paste0("M = ", ceiling(m)), 
               paste0("M = ", round(m, 1)))) 
eqFSTs <- eqFSTs %>% group_by(m) %>% do(helper(.)) %>% 
    mutate(s1 = "s = 0\nequilibrium")  %>%
    mutate(sty = "pred")



Ft <- function(m, t){
    return(0.5 + (-0.5)*exp(-2*m*t))
}

ft <- function(s,t){
    return(exp(s*t)/(exp(s*t) + 10^5 - 1))
}

Fst <- function(m, s, ts){
    F <- Ft(m, ts)
    f <- ft(s,ts)
    return((((1 - 2*F)^2)*f)/(2 - 2*f*F*(1 - F) - f))
}


fstan <- foreach(s = c(.3, 3), .combine = "rbind")%do%{
    foreach(m = ms, .combine = "rbind")%do%{
        tbl_df(cbind(gen = 0:100, FST = Fst(m = m, s = s, t = c(0:100)))) %>%
            mutate(s1 = s, m = m)
    }
}


fstvals <- fsts %>% gather(gen, freq, -N, -s1, -m, -i) %>% 
    mutate(gen = as.numeric(gsub("FST_", "", gen))) %>%
    group_by(s1, m, gen) %>%
    mutate(av = mean(freq))



simvan <- left_join(fstvals, fstan) %>% filter(s1 > 0) %>% 
    gather(type, FST, -N, -s1, -m, -i, -gen) %>% 
    mutate(type = ifelse(type == "freq", "Simulation", 
                  ifelse(type == "av", "Average of simulations", 
                  ifelse(type == "FST", "Prediction", NA )))) %>%
    mutate(mode = factor(ifelse(type == "Prediction", "Predicted", "Simulated"),
                levels = c("Simulated", "Predicted")))  %>%
    mutate(ld = factor(ifelse(type == "Simulation", "light", "dark"),
                levels = c("light", "dark"))) %>%
    ungroup() %>%
    mutate(m = ifelse(N*m > 1, paste0("M = ", ceiling(N*m)), 
               paste0("M = ", round(N*m, 1)))) %>% 
    mutate(s1 = paste0("s = ", s1)) %>% 
    ggplot() + 
    geom_line(aes(x = gen, y = FST, group = paste(s1, i, type), 
                  lty = mode, alpha = ld)) +
    facet_grid(s1~m, scales = "free") +guides(alpha = FALSE, color = FALSE) +
    theme_gen() + 
    theme(legend.position = c(0.9, 0.9), 
          legend.background = element_rect(color = NA), 
          legend.title = element_blank()) + 
    labs(lty = "", y = expression("F"['ST']), x = "Time (Generations)")  



pdf("../graphs/sims.v.an.pdf", width = 7.87*.75, height = 4*.85)
print(simvan)
dev.off()










