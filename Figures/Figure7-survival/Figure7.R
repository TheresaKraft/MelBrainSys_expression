# Set directory to source file loaction
setwd(paste0(figureDirectory, "Figure7-survival" ) )

library(ggplot2)
library(cowplot)


load("plotListSurvival.RData")



png(filename = "Figrue7.png", height = 1000, width = 1000, res = 90)
plot_grid(plotList[[9]],plotList[[2]],plotList[[3]],plotList[[5]],plotList[[11]],plotList[[7]],
          plotList[[1]],plotList[[4]],plotList[[8]],plotList[[6]],plotList[[10]],plotList[[12]],
          align = "v", ncol = 3) 
dev.off()
