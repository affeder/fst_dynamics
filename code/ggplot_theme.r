theme_gen <- function(base_size = 11, base_family = "") {
  half_line <- base_size / 2

  theme(
    # Elements in this first block aren't used directly, but are inherited
    # by others
    line =  element_line(colour = "black", size = 0.5, linetype = 1,
        lineend = "butt"),
    rect =  element_rect(fill = NA, colour = "black",
        size = 0.5, linetype = 1),
    text =  element_text(
        family = base_family, face = "plain",
        colour = "black", size = base_size,
        lineheight = 0.9, hjust = 0.5, vjust = 0.5, angle = 0,
        margin = margin(), debug = FALSE
    ),

    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = rel(0.8), colour = "black"),
    axis.text.x = element_text(margin = margin(t = 0.8 * half_line / 2), vjust = 1),
    axis.text.y = element_text(margin = margin(r = 0.8 * half_line / 2), hjust = 1),
    axis.ticks = element_line(colour = "black"),
    axis.ticks.length =  unit(half_line / 2, "pt"),
    axis.title.x = element_text(
        margin = margin(t = 0.8 * half_line, b = 0.8 * half_line / 2)
                         ),
    axis.title.y = element_text(
                           angle = 90,
                           margin = margin(r = 0.8 * half_line, l = 0.8 * half_line / 2)
                         ),

    legend.background =  element_rect(fill = NA, colour = NA),
    legend.key =         element_rect(fill = NA, colour = NA),
    legend.key.size =    unit(1.2, "lines"),
    legend.key.height =  NULL,
    legend.key.width =   NULL,
    legend.text =        element_text(size = rel(0.8)),
    legend.text.align =  NULL,
    legend.title =       element_text(hjust = 0),
    legend.title.align = NULL,
    legend.position =    "right",
    legend.direction =   NULL,
    legend.justification = "center",
    legend.box =         NULL,

    panel.background =   element_rect(fill = NA, colour = NA),
    panel.border =       element_rect(fill = NA, colour = "black"),
    panel.grid.major =   element_line(colour =  NA),
    panel.grid.minor =   element_line(colour = NA, size = 0.25),
    panel.ontop    =     FALSE,
    panel.spacing = unit(0.2, "lines"),
      
    strip.background =   element_rect(fill = NA, colour = NA),
    strip.text =         element_text(colour = "black", size = rel(0.8)),
    strip.text.x =       element_text(margin = margin(t = half_line, b = half_line)),
    strip.text.y =       element_text(angle = -90, margin = margin(l = half_line, r = half_line)),
    strip.switch.pad.grid = unit(0.1, "cm"),
    strip.switch.pad.wrap = unit(0.1, "cm"),

    plot.background = element_rect(colour = "white"),
    plot.title = element_text(
        size = rel(1.2),
        margin = margin(b = half_line * 1.2)
    ),
    plot.margin =   margin(half_line, half_line, half_line, half_line),

    complete = TRUE
  )
}
