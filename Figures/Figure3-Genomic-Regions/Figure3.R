# Set directory to source file loaction
setwd(paste0(figureDirectory, "Figure3-Genomic-Regions" ) )

library(ggplot2)
library(cowplot)


load(file = "barplotOverviewDF.RData")
load(file = "pathwayGenePlotDf.RData")
load(file = "immunePlotDf.RData")

colorsExp <- c(rep("blue3",8),rep("forestgreen",6),rep("gold3",3),rep("darkorchid",2),rep("hotpink1",1), "cyan3")

plotOverview <-ggplot(barplotOverviewDF, aes(color=Expression,y=value, x=reorder(PatientID,-ranking)), ylim = c(0,4)) + 
  geom_bar(position = "identity",stat="identity", fill="white",) +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill="white", size=.5, linetype="solid",color="black"),
        legend.title = element_blank(), legend.direction = "horizontal") +
  scale_y_continuous(breaks = c(-2000,-1500,-1000,-500,0,500,1000,1500,2000), labels = c("2000","1500","1000","500","0","500","1000","1500","2000"), limits = c(-2300,2300))+
  ggtitle("") + 
  ylab("Number of genes \n Decreased  \t \t Increased") +
  geom_vline(xintercept = c(8.5,14.5,17.5,19.5,20.5,21.5)) +  geom_hline(yintercept = 0)   + scale_color_manual(values= c("cornflowerblue","brown3") ) 


immunePlotDfTable <-   as.data.frame(table(immunePlotDf$association[immunePlotDf$significancePos == "X" | immunePlotDf$significanceNeg == "X"]))
freqImmuns <- immunePlotDfTable$Var1[immunePlotDfTable$Freq > 9]
immunePlotDfFreq <- immunePlotDf[immunePlotDf$association %in% freqImmuns,]

plotImmune <-    ggplot(immunePlotDfFreq, aes(fill=association,y=value, x=reorder(PatientID,-ranking)), ylim = c(0,4)) + 
  geom_bar(width = 0.7, stat="identity", position = position_dodge(width = 0.7), colour="white") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1,color = colorsExp),
        axis.title.x=element_blank(),legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill="white", size=.5, linetype="solid",color="black"),
        legend.title = element_blank(),legend.direction = "horizontal") +
  ggtitle("") +  ylab(" \n % affected genes \n Decreased \t \t  \t Increased") + 
  geom_text(aes(label=significancePos),stat='identity',position=position_dodge(0.68),vjust=0, size = 3, colour = "red",fontface = "bold")  +
  geom_text(aes(label=significanceNeg),stat='identity',position=position_dodge(0.68),vjust=1, size = 3, colour = "red",fontface = "bold")  +
  scale_y_continuous(breaks = c(-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5), labels = c("60","50","40","30","20","10","0","10","20","30","40","50"), limits = c(-0.60,0.5))+
  geom_vline(xintercept = c(8.5,14.5,17.5,19.5,20.5,21.5)) +  geom_hline(yintercept = 0) + scale_fill_manual(values= c("bisque3","blueviolet","darkolivegreen3","plum1") )


pathwayGenePlotDfTable <-   as.data.frame(table(pathwayGenePlotDf$association[pathwayGenePlotDf$significancePos == "X" | pathwayGenePlotDf$significanceNeg == "X"]))
freqPaths <- pathwayGenePlotDfTable$Var1[pathwayGenePlotDfTable$Freq >= 10]
pathwayGenePlotDfFreq <- pathwayGenePlotDf[pathwayGenePlotDf$association %in% freqPaths,]


pathwayGenePlotDfFreq$association[pathwayGenePlotDfFreq$association == "pi3k_akt"] <- "Pi3K/Akt"
pathwayGenePlotDfFreq$association[pathwayGenePlotDfFreq$association == "jak_stat"] <- "JAK-STAT"

plotPath <-    ggplot(pathwayGenePlotDfFreq, aes(fill=association,y=value, x=reorder(PatientID,-ranking)), ylim = c(0,4)) + 
  geom_bar(width = 0.7, stat="identity", position = position_dodge(width = 0.7), colour="white") + 
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),legend.justification=c(1,1), legend.position=c(1,1),
        legend.background = element_rect(fill="white", size=.5, linetype="solid",color="black"),
        legend.title = element_blank(), legend.direction = "horizontal") +
  ggtitle("") +  ylab("% affected genes \n Decreased  \t \t \t Increased ") + 
  geom_text(aes(label=significancePos),stat='identity',position=position_dodge(0.68),vjust=0, size = 3, colour = "red")  +
  geom_text(aes(label=significanceNeg),stat='identity',position=position_dodge(0.68),vjust=1, size = 3, colour = "red")  +
  scale_y_continuous(breaks = c(-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4), labels = c("50","40","30","20","10","0","10","20","30","40"), limits = c(-0.6,0.45))+
  geom_vline(xintercept = c(8.5,14.5,17.5,19.5,20.5,21.5)) +  geom_hline(yintercept = 0) + scale_fill_manual(values= c("darkorange","chartreuse","cornflowerblue","deeppink","aquamarine2","deepskyblue4"))





png(file ="Figure3.png", width = 2400, height = 3600, res = 300)
plot_grid(plotOverview, plotPath, plotImmune, align = "v", nrow = 3, rel_heights = c(31/100, 31/100, 38/100), label_size = 15, labels = c('A', 'B', 'C')) 

dev.off()


