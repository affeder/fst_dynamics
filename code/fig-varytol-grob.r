#Figure to show the effect of tolerance on estimation of m 

#First, read in the data
varytol <- tbl_df(read.table("../out/sumStats-varytol.txt", header = FALSE))
names(varytol) <- c("N", "n", "s1", "m", "i", "type", "offs", "addhund",  "tol",
                   "s", "val", "q", "var", 
                   "med_dist", "mean_dist", "med_dist_sq", "mean_dist_sq")

varytol <- varytol %>% filter(tol  >= .001)

Ms <- unique(varytol %>% mutate(M = N*m) %>% select(M))$M
factnames <- paste0("M = ", Ms)
lineinf <- tbl_df(cbind(Ms, factnames))
names(lineinf) <- c("M", "mname")
lineinf <- lineinf %>% mutate(mname = factor(mname, levels = factnames)) %>%
    mutate(M = as.numeric(M))

tm.50 <- varytol %>% filter(var == "m") %>%
    do(heatMapHelper(., qs = .95, testvals = lseq(.0001, .2, length.out = 100)))

xbr <- c(.001, .01, .1)
xlims <- c(.0005, .15)

s5.1<-   tm.50 %>% group_by(N, n, s1, m, type, offs, addhund, tol, s, testm) %>% 
    summarize(pInInt = mean(inInt)) %>% 
    mutate(mname = ifelse(N*m > 1, paste0("M = ", ceiling(N*m)), 
               paste0("M = ", round(N*m, 1)))) %>%
    mutate(mname = factor(mname, levels = factnames)) %>% 
    mutate(sname = paste0("s = ", s1)) %>%
    ggplot() + geom_tile(aes(x = tol, y = testm*N, fill = pInInt, color = pInInt)) + 
        scale_y_log10(breaks = c(100,1000, 10000)) + 
        scale_x_log10(breaks = xbr, limits = xlims) +
        geom_hline(aes(yintercept = M), data = lineinf) +
        facet_wrap(~ mname, nrow = 1) +
        scale_fill_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5) +
        scale_color_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5, guide = "none") +
        labs(y = "Test M", x = element_blank()) + 
          #,fill = "P(value in\n95% posterior)\n" ) + theme_linedraw() + 
                 theme_gen() + theme(legend.position="none") +
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank())



s5.2 <-  varytol %>% filter(var == "m") %>%
        do(heatMapHelper(., qs = .95, testvals = c(.001, .005, .01, .05))) %>%
        filter(testm == m) %>% 
    mutate(foundTrue = m >= low & m <= high) %>%
    group_by(N, n, s1, m, type, offs, addhund, tol, s) %>% 
    summarize(pFound = mean(foundTrue)) %>% 
    mutate(mname = ifelse(N*m > 1, paste0("M = ", ceiling(N*m)), 
               paste0("M = ", round(N*m, 1)))) %>%
    mutate(mname = factor(mname, levels = factnames)) %>% 
    ggplot(aes(x = tol, y = pFound)) + geom_point() + geom_line() +
    scale_x_log10(breaks = xbr, limits = xlims) + facet_grid(~ mname) +
    labs( y = "P(truth in 95% posterior)", x= element_blank()) +
    theme_gen()+
    theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) +
    theme(strip.background = element_blank(),
    strip.text.x = element_blank())



s5.3 <- varytol %>% filter(var == "m", q == .5) %>% #, addhund == 0) %>% 
    group_by(N, m, s1, n, offs, tol) %>% 
    mutate(med = mean(mean_dist_sq)) %>% ungroup() %>%
    mutate(M = N*m) %>%
    mutate(sname = paste0("s = ", s1)) %>%
    mutate(mname = ifelse(N*m > 1, paste0("M = ", ceiling(N*m)), 
               paste0("M = ", round(N*m, 1)))) %>%
    mutate(mname = factor(mname, levels = factnames)) %>%
    ggplot(aes(x = tol, y = med)) + 
    geom_jitter(aes(x = tol, y= mean_dist_sq), size = .1, col = "grey") +
    geom_point() + geom_line() +
    scale_x_log10(breaks = xbr, limits = xlims) +
    facet_grid(~mname) + 
    scale_color_manual(values = baseCols[c(2, 1, 5)]) + 
    labs(x = "Tolerance", y = "MSE", 
             col = element_blank()) +
                 theme_gen() + 
    theme(strip.background = element_blank(),
          strip.text.x = element_blank())



neededcols <- ncol(ggplotGrob(s5.1))

g5.2 <- ggplotGrob(s5.2)
while(ncol(g5.2) < neededcols){
    g5.2 <- gtable_add_cols(g5.2, unit(0,"mm"))
}

matchlay <- ggplotGrob(s5.1)$layout
tmplay <- g5.2$layout
panelinds <- grep("panel", tmplay[['name']] )
tmplay[panelinds,'l'] <- matchlay[grep("panel", matchlay[['name']] ), 'l']
tmplay[panelinds,'r'] <- matchlay[grep("panel", matchlay[['name']] ), 'r']

axisinds <- grep("axis-[tb]", tmplay[['name']] )
axreg <- paste(tmplay[['name']][axisinds], collapse = "|")
tmplay[axisinds,'l'] <- matchlay[grep(axreg, matchlay[['name']] ), 'l']
tmplay[axisinds,'r'] <- matchlay[grep(axreg, matchlay[['name']] ), 'r']

g5.2$layout <- tmplay


g5.3 <- ggplotGrob(s5.3)
while(ncol(g5.3) < neededcols){
    g5.3 <- gtable_add_cols(g5.3, unit(0,"mm"))
}

matchlay <- ggplotGrob(s5.1)$layout
tmplay <- g5.3$layout
panelinds <- grep("panel", tmplay[['name']] )
tmplay[panelinds,'l'] <- matchlay[grep("panel", matchlay[['name']] ), 'l']
tmplay[panelinds,'r'] <- matchlay[grep("panel", matchlay[['name']] ), 'r']

axisinds <- grep("axis-[tb]", tmplay[['name']] )
axreg <- paste(tmplay[['name']][axisinds], collapse = "|")
tmplay[axisinds,'l'] <- matchlay[grep(axreg, matchlay[['name']] ), 'l']
tmplay[axisinds,'r'] <- matchlay[grep(axreg, matchlay[['name']] ), 'r']

tmplay[which(tmplay$name == 'xlab-b'), 'l'] <- 9


g5.3$layout <- tmplay



pdf(paste0("../graphs/s5.pdf"), width = 7, height = 5)
grid.draw(rbind(ggplotGrob(s5.1),g5.2, g5.3, size = "first"))
dev.off()

