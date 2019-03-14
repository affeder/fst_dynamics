
#figure to show the effect of varying the time points sampled
varytp <- tbl_df(read.table("../out/sumStats-varytp.txt", header = FALSE))
names(varytp) <- c("N", "n", "s1", "m", "i", "type", "offs", "addhund",  "tol",
                   "s", "val", "q", "var", 
                   "med_dist", "mean_dist", "med_dist_sq", "mean_dist_sq")


factnames <- c(paste0("T = ", unique((varytp %>% select(offs))$offs)), " ")

tm.50 <- varytp %>% filter(var == "m", addhund == 0) %>%
    do(heatMapHelper(., qs = .95, testvals = lseq(.0001, .2, length.out = 100))) %>% 
    filter(n == 100) %>%
    group_by(N, n, s1, m, type, offs, addhund, tol, s, testm) %>% 
    summarize(pInInt = mean(inInt)) %>% 
        mutate(tpname = factor(paste0("T = ", offs), factnames)) %>% 
    mutate(M = N*m, estM = N*testm) %>%
    mutate(sname = paste0("s = ", s1))

lims <- c(5, 20000)/Nval

factnames <- paste0("T = ", unique((varytp %>% select(offs))$offs))
tm.50 <- varytp %>% filter(var == "m", addhund == 0) %>%
    do(heatMapHelper(., qs = .95, testvals = lseq(.0001, .2, length.out = 15)))

toPlot <- tm.50 %>% 
    group_by(N, n, s1, m, type, offs, addhund, tol, s, testm) %>% 
    summarize(pInInt = mean(inInt)) %>% 
    mutate(nname = factor(paste0("n = ", n), factnames)) %>%
    mutate(M = m, estM = testm) %>%
        mutate(sname = paste0("s = ", s1)) %>%
            mutate(tpname = factor(paste0("T = ", offs), factnames)) %>% 
        mutate(M = N*m, estM = N*testm) %>%
        mutate(sname = paste0("s = ", s1))



s3.1 <- toPlot %>% 
        ggplot() + geom_tile(aes(x = M/Nval, y = estM/Nval, fill = pInInt, color = pInInt)) + 
        scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, 
                      labels=fancy_scientific) + 
        scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, 
                      labels=fancy_scientific, 
                      sec.axis = sec_axis(~.*Nval, name = "True M (N*m)", 
                          breaks = c(10, 100, 1000, 10000)), lim = lims) +
        facet_grid(tpname ~ sname ) + 
        scale_fill_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5) +
        scale_color_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5, guide = "none") +
        geom_point(aes(x = M/Nval, y = M/Nval), pch = "-", size = 3) + 
        labs(x = "True m", y = "Estimated m", 
             fill = "P(value in\n95% posterior)\n" ) +
        theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
        theme_gen()  + theme( strip.background = element_blank(),
                                 strip.text.x = element_blank()) + 
        geom_text(aes(x = 10/Nval, y= 1, #y = 0.1,
                      label = sname), hjust = 0, 
           data = toPlot %>% ungroup() %>% select(sname) %>% distinct(), col = "white")



s3.2 <- varytp %>% filter(var == "m", q == .5, addhund == 0) %>% 
    group_by(N, m, s1, n, offs) %>% summarize(med = mean(mean_dist_sq)) %>% ungroup() %>%
    mutate(M = N*m) %>%
    mutate(sname = paste0("s = ", s1)) %>%
    mutate(tpname = factor(paste0("T = ", offs), factnames)) %>% 
    ggplot(aes(x = M/Nval, y = med, col = tpname)) + 
        geom_point() + geom_line() +
    scale_x_log10(breaks = c(10, 100, 1000, 10000)/Nval, limits = lims, labels=fancy_scientific) +
    facet_grid(~sname) + 
    scale_color_manual(values = baseCols[c(2, 1, 5)]) + 
    labs(x = "True m", y = "MSE", 
             col = element_blank()) +
                 theme_gen()



g2 <- ggplotGrob(s3.2)
g2 <- gtable_add_cols(g2, unit(0,"mm"))

tmplay <- g2$layout
tmplay[which(tmplay[['name']] == "guide-box"),'l'] <- 13
g2$layout <- tmplay


pdf(paste0("../graphs/s3.pdf"), width = 7, height = 6)
grid.draw(rbind(ggplotGrob(s3.1),g2, size = "first"))
dev.off()


#Ok, now let's see how each comparison looks when adding in 100 or not


s4 <- varytp %>% filter(var == "m", q == .5) %>% 
    group_by(N, m, s1, n, offs, addhund) %>% 
    summarize(med = mean(mean_dist_sq)) %>% ungroup() %>%
    mutate(M = N*m) %>%
    mutate(sname = paste0("s = ", s1)) %>%
    mutate(tpname = factor(paste0("T = ", offs), factnames))%>% 
    mutate(addhundname = factor(ifelse(addhund == 1, "Sample at t = 100", "No sample at t = 100")))%>% 
    ggplot(aes(x = M/Nval, y = med, col = addhundname)) + 
        geom_point() + geom_line() +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, labels=fancy_scientific) +
    facet_grid(tpname~sname) + 
    scale_color_manual(values = baseCols[c(2, 5)]) + 
    labs(x = "True m", y = "MSE", 
             col = element_blank()) +
                 theme_gen()


pdf(paste0("../graphs/s4.pdf"), width = 7, height = 4)
print(s4)
dev.off()
