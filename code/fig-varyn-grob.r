
#figure to show the effect of varying n
varyn <- tbl_df(read.table("../out/sumStats-varyn.txt", header = FALSE))

names(varyn) <- c("N", "n", "s1", "m", "i", "type", "offs", "addhund",  "tol",
                   "s", "val", "q", "var", 
                   "med_dist", "mean_dist", "med_dist_sq", "mean_dist_sq")


factnames <- c("n = 30", "n = 100", "n = 500")
tm.50 <- varyn %>% filter(var == "m") %>%
    do(heatMapHelper(., qs = .95, testvals = lseq(.0001, .2, length.out = 200)))

f2 <- tm.50 %>% filter(n == 30) %>%
    group_by(N, n, s1, m, type, offs, addhund, tol, s, testm) %>% 
    summarize(pInInt = mean(inInt)) %>% 
    mutate(nname = factor(paste0("n = ", n), factnames)) %>%
    mutate(M = N*m, estM = N*testm) %>%
    mutate(sname = paste0("s = ", s1)) %>%
    ggplot() + geom_tile(aes(x = M, y = estM, fill = pInInt, color = pInInt)) + 
        scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) + 
        scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)) +
        facet_wrap(~sname ) +
        scale_fill_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5) +
        scale_color_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5, guide = "none") +
        geom_point(aes(x = M, y = M), pch = "-", size = 3) + 
        labs(x = "True M", y = "Estimated M", 
             fill = "P(value in\n95% posterior)\n" ) +
        theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
            theme_gen()

Nval <- unique((tm.50 %>% select(N))$N)


toPlot <- tm.50 %>% filter(n == 30) %>%
    group_by(N, n, s1, m, type, offs, addhund, tol, s, testm) %>% 
    summarize(pInInt = mean(inInt)) %>% 
    mutate(nname = factor(paste0("n = ", n), factnames)) %>%
    mutate(M = m, estM = testm) %>%
        mutate(sname = paste0("s = ", s1))

f2 <- toPlot %>%  ggplot() + 
    geom_tile(aes(x = M, y = estM, fill = pInInt, color = pInInt)) + 
        scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, 
                      labels=fancy_scientific) + 
        scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, 
                      labels=fancy_scientific, 
                      sec.axis = sec_axis(~.*Nval, name = "True M (N*m)", breaks = c(10, 100, 1000, 10000))) +
#        geom_abline(slope = 1, intercept = 0) +
        facet_wrap(~sname ) +
        scale_fill_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5) +
        scale_color_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5, guide = "none") +
        geom_point(aes(x = M, y = M), pch = "-", size = 3) + 
        labs(x = "True m", y = "Estimated m", 
             fill = "P(value in\n95% posterior)\n" ) +
        theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
            theme_gen()  + theme( strip.background = element_blank(),
                                 strip.text.x = element_blank()) + 
        geom_text(aes(x = 10/Nval, y= .1, 
                      label = paste(c("A.", "B.", "C."), sname)), hjust = 0, 
           data = toPlot %>% ungroup() %>% select(sname) %>% distinct(), col = "white")



pdf(paste0("../graphs/f3.pdf"), width = 7, height = 2.75)
print(f2)
dev.off()


lims <- c(5, 20000)/Nval

tm.50 <- varyn %>% filter(var == "m") %>%
    do(heatMapHelper(., qs = .95, testvals = lseq(.0001, .2, length.out = 100)))

toPlot <- tm.50 %>% 
    group_by(N, n, s1, m, type, offs, addhund, tol, s, testm) %>% 
    summarize(pInInt = mean(inInt)) %>% 
    mutate(nname = factor(paste0("n = ", n), factnames)) %>%
    mutate(M = m, estM = testm) %>%
        mutate(sname = paste0("s = ", s1))

s2.1 <- toPlot %>%
    ggplot() + geom_tile(aes(x = M, y = estM, fill = pInInt, color = pInInt)) + 
        scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, 
                      labels=fancy_scientific) + 
        scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, 
                      labels=fancy_scientific, 
                      sec.axis = sec_axis(~.*Nval, name = "True M (N*m)", breaks = c(10, 100, 1000, 10000)), limits = lims) + 
        facet_grid(nname ~sname ) +
        scale_fill_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5) +
        scale_color_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5, guide = "none") + 
        geom_point(aes(x = M, y = M), pch = "-", size = 3) + 
        labs(x = "True m", y = "Estimated m", 
             fill = "P(value in\n95% posterior)\n" ) +
            theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
        theme_gen()  + theme( strip.background = element_blank(),
                                 strip.text.x = element_blank()) + 
        geom_text(aes(x = 10/Nval, y= .1, 
                      label = sname), hjust = 0, 
           data = toPlot %>% ungroup() %>% select(sname) %>% distinct(), col = "white")

toPlot.p2 <- varyn %>% filter(var == "m", q == .5) %>% 
    group_by(N, m, s1, n) %>% summarize(med = mean(mean_dist_sq)) %>% ungroup() %>%
    mutate(M = m) %>%
    mutate(sname = paste0("s = ", s1)) %>%
    mutate(nname = factor(paste0("n = ", n), factnames))

s2.2 <- toPlot.p2 %>%
    ggplot(aes(x = M, y = med, col = nname)) + 
        geom_point() + geom_line() +
    scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, limits = lims, 
                  labels=fancy_scientific, 
                  sec.axis = sec_axis(~.*Nval, name = "True M (N*m)", 
                      breaks = c(10, 100, 1000, 10000))
                  ) +
    facet_grid(~sname) + 
    scale_color_manual(values = baseCols[c(2, 1, 5)]) + 
    labs(x = "True M", y = "MSE", 
             col = element_blank()) +
    theme_gen() + theme( strip.background = element_blank(),
                                 strip.text.x = element_blank()) + 
    geom_text(aes(x = 10/Nval, y = 0.8, #y= 1.20, 
                      label = sname), col = "black", hjust = 0, vjust = 1, 
    data = toPlot.p2 %>% ungroup() %>% select(sname) %>% distinct())


g2 <- ggplotGrob(s2.2)
g2 <- gtable_add_cols(g2, unit(0,"mm"))

tmplay <- g2$layout
tmplay[which(tmplay[['name']] == "guide-box"),'l'] <- 13
g2$layout <- tmplay

pdf(paste0("../graphs/s2.pdf"), width = 7, height = 7)
grid.draw(rbind(ggplotGrob(s2.1),g2, size = "first"))
dev.off()
