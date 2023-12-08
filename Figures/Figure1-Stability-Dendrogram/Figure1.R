# Set directory to source file loaction
setwd(paste0(figureDirectory, "Figure1-Stability-Dendrogram" ) )

library(RColorBrewer)
library(stringr)
library(pvclust)
library(dendextend)

load(file = "valuesManhD2PV.RData")


patientColorManhD2 <- valuesManhD2PV$hclust$labels[valuesManhD2PV$hclust$order]
patientColorManhD2[str_detect(patientColorManhD2,"P106")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[1]
patientColorManhD2[str_detect(patientColorManhD2,"P13")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[2]
patientColorManhD2[str_detect(patientColorManhD2,"P111")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[3]
patientColorManhD2[str_detect(patientColorManhD2,"P108")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[4]
patientColorManhD2[str_detect(patientColorManhD2,"P77")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[5]
patientColorManhD2[str_detect(patientColorManhD2,"P42")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[6]
patientColorManhD2[str_detect(patientColorManhD2,"P8")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[7]
patientColorManhD2[str_detect(patientColorManhD2,"P4")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[8]
patientColorManhD2[str_detect(patientColorManhD2,"P18")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[9]
patientColorManhD2[str_detect(patientColorManhD2,"P101")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[10]
patientColorManhD2[str_detect(patientColorManhD2,"P39")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[11]
patientColorManhD2[str_detect(patientColorManhD2,"P3")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[12]
patientColorManhD2[str_detect(patientColorManhD2,"P78")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[13]
patientColorManhD2[str_detect(patientColorManhD2,"P16")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[14]
patientColorManhD2[str_detect(patientColorManhD2,"P74")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[15]
patientColorManhD2[str_detect(patientColorManhD2,"P107")] <- colorRampPalette(brewer.pal(8, "Oranges"))(16)[16]

    
tissueColorManhD2 <- valuesManhD2PV$hclust$labels[valuesManhD2PV$hclust$order]
tissueColorManhD2[str_detect(tissueColorManhD2,"_B")] <- "darkgrey"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Lym")] <- "forestgreen"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Lu")] <- "blue3"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Liv")] <- "deeppink2"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Ski")] <- "gold3"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Sof")] <- "darkorchid"
tissueColorManhD2[str_detect(tissueColorManhD2,"_Smi")] <- "darkseagreen3"


labelsMandD2 <- valuesManhD2PV$hclust$labels

labelsMandD2[str_detect(labelsMandD2, "P13" ) |  str_detect(labelsMandD2, "P108") |
         str_detect(labelsMandD2,"P42" )] <- paste0("*",labelsMandD2[str_detect(labelsMandD2, "P13" ) |  str_detect(labelsMandD2, "P108") |
                                                           str_detect(labelsMandD2,"P42" )])

clusterColors <- c(rep("maroon2",8),rep("burlywood3",6),rep("yellow",23))


png(file ="Figure1.png", width = 3000, height = 2000, res = 300 )
par(mar = c(3.3,3,0,0))
plot(valuesManhD2PV, hang = -1, cex = 0.9, print.pv = c("au"), main ="", axes = F,    
     print.num = F, labels =labelsMandD2 )
colored_bars(colors = cbind(clusterColors,patientColorManhD2, tissueColorManhD2),rowLabels = c("Cluster","Patient","Tissue"))
dev.off()

