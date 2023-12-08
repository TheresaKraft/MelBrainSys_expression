# Set directory to source file loaction
setwd(paste0(figureDirectory, "Figure4-histologicalRegions" ) )

library(stringr)
library(dendextend)
library(RColorBrewer)
library(ggpubr)



load("valuesManhD2PVPairs.RData")



patientColorManhD2 <- valuesManhD2PVPairs$hclust$labels[valuesManhD2PVPairs$hclust$order]
patientColorManhD2[str_detect(patientColorManhD2,"P42")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[1]
patientColorManhD2[str_detect(patientColorManhD2,"P106")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[2]
patientColorManhD2[str_detect(patientColorManhD2,"P16")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[3]
patientColorManhD2[str_detect(patientColorManhD2,"P18")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[4]
patientColorManhD2[str_detect(patientColorManhD2,"P4")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[5]
patientColorManhD2[str_detect(patientColorManhD2,"P101")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[6]
patientColorManhD2[str_detect(patientColorManhD2,"P78")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[7]
patientColorManhD2[str_detect(patientColorManhD2,"P74")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[8]
patientColorManhD2[str_detect(patientColorManhD2,"P111")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[9]
patientColorManhD2[str_detect(patientColorManhD2,"P77")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[10]
patientColorManhD2[str_detect(patientColorManhD2,"P8")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[11]
patientColorManhD2[str_detect(patientColorManhD2,"P3")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[12]
patientColorManhD2[str_detect(patientColorManhD2,"P13")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[13]
patientColorManhD2[str_detect(patientColorManhD2,"P108")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[14]
patientColorManhD2[str_detect(patientColorManhD2,"P107")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[15]
patientColorManhD2[str_detect(patientColorManhD2,"P39")] <- colorRampPalette(brewer.pal(8, "Oranges"))(17)[16]


tissueColorManhD2 <- valuesManhD2PVPairs$hclust$labels[valuesManhD2PVPairs$hclust$order]
tissueColorManhD2[str_detect(tissueColorManhD2,"_BLym")] <- "forestgreen"
tissueColorManhD2[str_detect(tissueColorManhD2,"_BLu")] <- "blue3"
tissueColorManhD2[str_detect(tissueColorManhD2,"_BLiv")] <- "deeppink2"
tissueColorManhD2[str_detect(tissueColorManhD2,"_BSki")] <- "gold3"
tissueColorManhD2[str_detect(tissueColorManhD2,"_BSof")] <- "darkorchid"
tissueColorManhD2[str_detect(tissueColorManhD2,"_BSmi")] <- "darkseagreen3"




png(file ="Figure4_A.png", width = 2000, height =1500, res = 300 )
par(mar = c(3,3,0,0))
plot(valuesManhD2PVPairs, hang = -1, cex = 0.9, print.pv = c("au"), main ="", axes = F,
     print.num = F) 
colored_bars(colors = cbind(patientColorManhD2, tissueColorManhD2),rowLabels = c("Patient","Tissue"))
dev.off()

load("diffsBetweenSamplesTable.RData")


plot2 <- ggplot(diffsBetweenSamplesTable, aes(x = PairComparison, y = OverlapPerc, color = PairComparison)) + geom_boxplot(color = "grey") +   
  labs(x = "Pair origin", y = "Overlap of HMM predictions in %") + geom_point(position = position_jitter(width = 0.2)) + 
  theme(legend.position = "none") + scale_color_manual(values = c("darkorange","forestgreen"))

png(file ="Figure4_B.png", width = 1000, height =1500, res = 300 )
plot2
dev.off()
