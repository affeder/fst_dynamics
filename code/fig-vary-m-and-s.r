
#First, read in the data
svm <- tbl_df(read.table("../out/sumStats-svm.txt", header = FALSE))
names(svm) <- c("N", "n", "s1", "m", "i", "type", "offs", "addhund",  "tol",
                   "s", "val", "q", "var", 
                   "med_dist", "mean_dist", "med_dist_sq", "mean_dist_sq")

svm_joint <- svm


#Distance from truth
svmp <-  svm_joint%>%  filter(var == "m") %>% filter(q == .05) %>%
    group_by(N, n, s1, m, offs, addhund, tol) %>%
    summarize(dist = median(mean_dist_sq)) %>% 
    mutate(dist.norm = dist) %>%
#    mutate(dist.norm = ifelse(dist.norm > 2, 2, dist.norm)) %>% 
    ggplot() +
        geom_tile(aes(x = s1, y = N*m, fill = dist.norm)) +
        scale_y_log10() + scale_x_log10(breaks = c(.3, 1, 3)) + 
        labs(x = "Selection strength (s)", y = "Migration rate (M)", 
             fill = "Median MSE              ") +
        scale_fill_gradient2(low = baseCols[2], mid = "white", 
                             high = baseCols[5], midpoint= 0.5) + 
        theme_gen()



jh.approx <- function(s, m, N, f0){

#Ok, for a given s and m and N
#What is the waiting time for the establishing mutation? 
#m <- .01

    F <- .5 + (f0 - 0.5)*exp(-2*m*ts)
    f <- exp(s*ts)/(exp(s*ts) + 1/f0 - 1)

    frequencies <- tbl_df(cbind(ts, F, f, F*f, (1 - F)*f))
    names(frequencies) <- c("gen", "F", "f", "nl", "l")

#What does FST look like?
    fst.func <- function(tbl){

        f1 <- tbl$l
        f2 <- tbl$nl
        wt <- 1 - f1 - f2

        #Drawing from within a population
        pi.w <- 1 - f1^2 - f2^2 - wt^2

        #Drawing between two populations
        pi.b <- 1 - 2* f1 * f2 - wt^2
        fst <- (pi.b - pi.w)/pi.b
        if(is.na(fst)){fst = 0}
        return(tbl_df(tbl) %>% mutate(FST = fst))

    }


    fs <- frequencies %>% rowwise() %>% do(fst.func(.))  %>% filter(gen <= 100)
    toRet <- fs %>% gather(type, val, -gen)
    return(toRet)

}

ts <- 1:100
N <- 10^5
lout <- 15
ms <- lseq(.001, .1, length.out = lout)
ss <- lseq(.3, 3, length.out = lout)

fsts <-foreach(s = ss,.combine = "rbind") %do%{
    foreach(m = ms, .combine = "rbind") %do%{
        jh.approx(s, m, 10^5, 1/(s*N)) %>% mutate(s =s , m = m)
    }
}

nmlevs <- paste0("M=", rev(round(N*ms)))
tfcols <- c(redCol, blueCol)
names(tfcols) <- c("TRUE", "FALSE")
toPlot <- fsts %>% filter(type == "FST") %>%
    group_by(s, m) %>%
    mutate(maxfst = max(val)) %>% mutate(maxfstgen = which(val == maxfst)) %>% 
    mutate(fst30genlat = val[maxfstgen + 30]) %>%
    mutate(fstdec = maxfst - fst30genlat) %>% 
    mutate(dec = fstdec > .15) %>%
    ungroup()  %>% 
    mutate(s = factor(s)) %>%
    mutate(M = factor(paste0("M=",round(N*m)), levels = nmlevs)) 

trajs <- toPlot %>% 
        ggplot()  + 
        labs(x = "s", y = "M", 
             col = "0.15 FST decline\nin 30 gens", fill = "") + 
        facet_grid(M ~ s) +
    geom_line(aes(x = gen, y = val, col = dec)) +
   theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + 
   theme(strip.text.y = element_blank()) + 
   theme(strip.text.x = element_blank()) + 
       scale_color_manual(values = tfcols) +
    theme(panel.spacing = unit(0, "lines"))


svmp <-  svm_joint %>%  filter(var == "m") %>% filter(q == .05) %>%
    group_by(N, n, s1, m, offs, addhund, tol) %>%
    summarize(dist = median(mean_dist_sq)) %>% 
    mutate(dist.norm = dist)


s.coord <- (toPlot %>% filter(dec == TRUE))$s
m.coord <- (toPlot %>% filter(dec == TRUE))$m


smcoords <- toPlot %>% filter(dec == TRUE) %>% mutate(m = N*m) %>% 
    select(s, m) %>% distinct()
chullinds <- chull(as.matrix(smcoords))
smcoords[chullinds,] %>% mutate(s = as.numeric(s)) %>% 
    ggplot() + geom_polygon(aes(x = s, y = m), fill = NA, col = "black") +
        scale_x_log10() + scale_y_log10()


svmplot <- svmp %>%
    ggplot() +
        geom_tile(aes(x = s1, y = N*m, fill = dist.norm)) +
        scale_y_log10() + scale_x_log10(breaks = c(.3, 1, 3)) + 
        labs(x = "Selection strength (s)", y = "Migration rate (M)", 
             fill = "Median MSE              ") +
        scale_fill_gradient2(low = baseCols[2], mid = "white", 
                             high = baseCols[5], midpoint= 0.5) + 
        theme_gen() 

pdf(paste0("../graphs/fig-svm.pdf"), width = 4, height = 2.3)
print(svmplot)
dev.off()

pdf(paste0("../graphs/fig-trajs.pdf"), width = 4, height = 2.3)
print(trajs)
dev.off()

