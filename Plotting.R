library(ggplot2)
library(ggpubr)
plot_aes = function(p, size = 20) {
  q = p + theme(
    plot.title = element_text(size = size),
    axis.title.x = element_text(size = size),
    axis.title.y = element_text(size = size),
    axis.text.x = element_text(size = size),
    axis.text.y = element_text(size = size),
    legend.title = element_text(size = size),
    legend.text = element_text(size = size),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = 'white', colour = 'black'),
    axis.line = element_line(colour = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    )
  ) + grids()
  theme_set(theme_gray(base_size = size))
  return(q) 
}
