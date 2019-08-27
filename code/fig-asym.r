
#figure to show the effect of varying n
#varyn <- tbl_df(read.table("../out/sumStats-varyn.txt", header = FALSE))
vary.na <- tbl_df(read.table("../out/sumStats-asym-n100-long-ah.txt", header = FALSE))
vary.nonna  <- tbl_df(read.table("../out/sumStats-asym-n100-long.txt", header = FALSE))
varyn <- bind_rows(vary.nonna, vary.na)


names(varyn) <- c("N", "n", "s1", "m", "i", "type", "offs", "addhund",  "tol",
                   "s", "asym", "val", "q", "var", 
                   "med_dist", "mean_dist", "med_dist_sq", "mean_dist_sq")


factnames <- c("n = 30", "n = 100", "n = 500")
tm.50 <- varyn %>% filter(var == "m") %>%
    do(heatMapHelper(., qs = .95, testvals = lseq(.0001, .2, length.out = 200)))
Nval <- unique((tm.50 %>% select(N))$N)
lims <- c(5, 20000)/Nval


toPlot <- tm.50 %>% 
    group_by(N, n, s1, m, type, offs, addhund, tol, s, testm, asym) %>% 
    summarize(pInInt = mean(inInt)) %>% 
    mutate(nname = factor(paste0("n = ", n), factnames)) %>%
    mutate(M = m, estM = testm) %>%
        mutate(sname = paste0("s = ", s1))


f2 <-toPlot %>%
    ggplot() + geom_tile(aes(x = M, y = estM, fill = pInInt, color = pInInt)) + 
        scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, 
                      labels=fancy_scientific) + 
        scale_x_log10(breaks = c(1, 10, 100, 1000, 10000)/Nval, 
                      labels=fancy_scientific, 
                      sec.axis = sec_axis(~.*Nval, name = "True M (N*m)", 
                                          breaks = c(10, 100, 1000, 10000)), limits = lims) + 
    facet_grid( type ~ asym ) +
        scale_fill_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5) +
        scale_color_gradient2(low = baseCols[5], mid = "white", high = baseCols[2],
                             midpoint= 0.5, guide = "none") + 
        geom_point(aes(x = M, y = asym*M), pch = "-", size = 3, col = baseCols[1],
                   data = toPlot %>% filter(type == "m.asym") %>% filter(asym*M >= 1e-4)) + 
        geom_point(aes(x = M, y = M), pch = "-", size = 3) + 
        labs(x = "True m", y = "Estimated m", 
             fill = "P(value in\n95% posterior)\n" ) +
        theme(axis.text.x=element_text(angle = -45, hjust = 0)) + 
        theme_gen()  + theme( strip.background = element_blank(),
                                 strip.text.x = element_blank()) + 
        geom_text(aes(x = 10/Nval, y= .1, 
                      label = paste(type)), hjust = 0, 
                  data = toPlot %>% ungroup() %>% select(asym, type) %>% distinct(), col = "white")


pdf(paste0("../graphs/fig_asym.pdf"), width = 7, height = 6)
print(f2)
dev.off()


