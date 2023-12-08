# Set directory to source file loaction
setwd(paste0(figureDirectory, "Figure6-candHeatmap" ) )

library(pheatmap)
library(stringr)
library(grid)
library(RColorBrewer)


load("geneCandidates.RData")



candGenes_biolFunction_summary <- read.csv("candGenes_biolFunction_summary_validation.csv")

rownames(candGenes_biolFunction_summary) <- candGenes_biolFunction_summary$Gene
candGenes_biolFunction_summary$Gene <- NULL


candGenes_biolFunction_summary$Function[candGenes_biolFunction_summary$Function == "apoptosis"] <- "other/unknown"
candGenes_biolFunction_summary$Function[candGenes_biolFunction_summary$Function == "cell adhesion"] <- "other/unknown"
candGenes_biolFunction_summary$Function[candGenes_biolFunction_summary$Function == "cell migration"] <- "other/unknown"
candGenes_biolFunction_summary$Function[candGenes_biolFunction_summary$Function == "differentiation"] <- "other/unknown"
candGenes_biolFunction_summary$Function[candGenes_biolFunction_summary$Function == "metabolism"] <- "other/unknown"
candGenes_biolFunction_summary$Function[candGenes_biolFunction_summary$Function == "muscle movement"] <- "other/unknown"
candGenes_biolFunction_summary$Function[candGenes_biolFunction_summary$Function == "transcription"] <- "other/unknown"


colorPalette <- brewer.pal(11,"Spectral")

ann_colors = list(  
  Function = c("cell growth" = colorPalette[3], "immune response" = colorPalette[6], "signal transduction" = colorPalette[9], "other/unknown" = "lightgrey",
               "transport" = colorPalette[11]), 
  Fischer = c("increased"= "brown3","decreased" = "cornflowerblue", "not listed" = "darkgrey"), 
  Chen = c("increased"= "brown3", "decreased" = "cornflowerblue","not listed" = "darkgrey"), 
  Biermann = c("increased"= "brown3", "decreased" = "cornflowerblue","not listed" = "darkgrey")
)


htmp <- pheatmap(geneCandidates, color = colorRampPalette(c('cornflowerblue','white','brown3'))(length(seq(-6,6, by = 0.01))), breaks = seq(-6,6, by = 0.01), 
                 clustering_distance_rows = "manhattan", clustering_distance_cols = "manhattan", legend_breaks = c(-6,-4.5,-3,0,3,4.5,6), 
                 main = "", legend_labels = c("-6","  Decreased \n  expression","-3","0","3", "  Increased \n  expression","6"), show_rownames=T, 
                 treeheight_row = 25, fontsize_row = 7,
                 annotation_row  = candGenes_biolFunction_summary, annotation_colors = ann_colors,legend = F,annotation_legend = F)




labelshtmp <- htmp$tree_col$labels[htmp$tree_col$order]
labelshtmp[str_detect(labelshtmp ,"BLun")] <- "blue3"
labelshtmp[str_detect(labelshtmp , "BLym")] <- "forestgreen"
labelshtmp[str_detect(labelshtmp , "BSki")] <- "gold3"
labelshtmp[str_detect(labelshtmp , "BSof")] <- "darkorchid"
labelshtmp[str_detect(labelshtmp , "BLiv")] <- "hotpink1"
labelshtmp[str_detect(labelshtmp , "BSmi")] <- "cyan3"

cols= labelshtmp
htmp$gtable$grobs[[5]]$gp=gpar(col=cols)





png(filename ="Figure6_noLegend.png", height = 3300, width = 1700, res = 300)
htmp
dev.off()






