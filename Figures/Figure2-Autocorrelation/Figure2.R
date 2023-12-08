# Set directory to source file loaction
setwd(paste0(figureDirectory, "Figure2-Autocorrelation" ) )

library(ggplot2)
library(cowplot)

load(file = "allAcfPermutatedPlot.RData")

plotAutocorrelation <- ggplot(data = allAcfPermutatedPlot, mapping=aes(x=lag)) +
  geom_hline(yintercept = 0, col = "black") +  
  geom_ribbon(aes(ymin =  acfMedian - acfSd, ymax =  acfMedian + acfSd, fill = "Sd"), alpha = .5) +
  geom_ribbon(aes(ymin =  acfBiolMedian - sdBiol, ymax =  acfBiolMedian + sdBiol, fill = "Sd"), alpha = .5) +
  geom_point(mapping = aes(x=lag, y = acfMedian, colour = "Random"), size = 3) +   geom_line(mapping = aes(x=lag, y = acfMedian, colour = "Random")) + 
  geom_point(mapping = aes(x=lag, y = acfBiolMedian, colour = "Biological"), size = 3) +   geom_line(mapping = aes(x=lag, y = acfBiolMedian, colour = "Biological")) + 
  scale_x_continuous(breaks = c(1,5,10,15,20),minor_breaks =  seq(0, 20, by = 1)) + 
  scale_color_manual(values = c("red", "black")) + scale_fill_manual(values = c("burlywood1")) + 
  ylab("Autocorrelation") + xlab("Positional lag of neighboring genes")+ theme(legend.position = "none")


png(filename = "Figure2.png", width = 2200, height = 1100, res = 300)
plot_grid(plotAutocorrelation, NULL,  align = "h", nrow = 1, label_size = 15, labels = c('A', 'B')) 

dev.off()
