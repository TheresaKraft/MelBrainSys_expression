# Set directory to source file loaction
setwd(paste0(figureDirectory, "Figure5-GOenrichment" ) )

load("resultGO.RData")


load("numCandGens.RData")

resultPlus <- resultGO[resultGO$Expression.direction  == "Increased",]
resultPlus <- resultPlus[1:20,]

resultPlus$Description[resultPlus$Description == "calcium-dependent cell-cell adhesion via plasma membrane cell adhesion molecules"] <- "calcium-dependent cell-cell adhesion"

resultPlus$color <- ""
resultPlus$color[str_detect(resultPlus$Description, "synap") ] <- "synaptic processes"
resultPlus$color[ resultPlus$Description == "forebrain development" | resultPlus$Description == "gliogenesis" | resultPlus$Description == "astrocyte differentiation" | resultPlus$Description == "limbic system development" | 
                    resultPlus$Description == "glial cell differentiation"] <- "brain-specific cell development"
resultPlus$color[ resultPlus$Description == "behavior"| resultPlus$Description == "calcium ion-regulated exocytosis of neurotransmitter"| 
                    resultPlus$Description == "adult behavior"] <- "neuronal processes"

resultMinus <- resultGO[resultGO$Expression.direction  == "Decreased",]
resultMinus <- resultMinus[1:20,]


resultMinus$color <- ""
resultMinus$color[ resultMinus$Description == "humoral immune response" | resultMinus$Description == "leukocyte mediated immunity" | 
                     resultMinus$Description == "adaptive immune response" | resultMinus$Description == "immune response-regulating cell surface receptor signaling pathway" |
                     resultMinus$Description == "lymphocyte mediated immunity" | resultMinus$Description == "" ] <- "immune responses"
resultMinus$color[ resultMinus$Description == "regulation of lymphocyte activation" | resultMinus$Description == "regulation of leukocyte activation" | 
                     resultMinus$Description == "activation of immune response" |  resultMinus$Description == "positive regulation of immune response" |  
                     resultMinus$Description == "T cell activation" |  resultMinus$Description == "immune response-activating cell surface receptor signaling pathway" |  
                     resultMinus$Description == "immune response-activating signal transduction" |  resultMinus$Description == ""  ] <- "immune activation"
resultMinus$color[ resultMinus$Description == "leukocyte differentiation" | resultMinus$Description == "mononuclear cell differentiation" | 
                     resultMinus$Description == "" | resultMinus$Description == "" ] <- "immune-cell development"
resultMinus$color[ resultMinus$Description == "chemokine-mediated signaling pathway"| resultMinus$Description == "cellular response to chemokine" | 
                     resultMinus$Description == "response to chemokine" ] <- "chemokine-related processes"


plotPLus <- ggplot(resultPlus, aes(x = -log10(qvalue), y = reorder(Description,-qvalue), fill = color) ) + geom_bar(stat = "identity") + ylab("") + theme(legend.position = "none") + 
  scale_fill_manual(values = c("burlywood3","darkorange","chocolate4","brown2"))+ xlab(expression(-log[10]("q-value")))

plotMinus <- ggplot(resultMinus, aes(x = -log10(qvalue), y = reorder(Description,-qvalue), fill = color) ) + geom_bar(stat = "identity") + ylab("")  + theme(legend.position = "none")+ 
  scale_fill_manual(values = c("cyan3","darkgreen","darkolivegreen4","darkorchid4","deepskyblue4")) + xlab(expression(-log[10]("q-value")))



plotCand <- ggplot(numCandGens, aes(cutOff, genesCount, color = Expression)) + geom_point(cex = 3) +  geom_line(aes(cutOff,genesCount, color = Expression)) +
  scale_y_continuous(trans = "log10", name = "Number of gene candidates", breaks = c(1,10,100,1000,10000,10000,100000), labels = c(1,10,100,1000,10000,10000,"100000")) + 
  scale_x_continuous(name = "Number of patients cutoff", breaks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), minor_breaks = NULL)  +
  scale_color_manual(values = c("cornflowerblue","brown3")) + theme(legend.position = "none")+ 
  theme(legend.position = c(0.943,0.84), legend.background = element_rect(size=0.5, linetype="solid", 
                                                                          colour ="black"))  + labs(color = "Intracranial \nexpression") +
  geom_text_repel(aes(x = 16, y = 1), 
                  data = numCandGens[numCandGens$cutOff == 16 & numCandGens$genesCount == 1,], 
                  label ="CCL19", inherit.aes = F, color = "cornflowerblue",nudge_y = 2.1) + 
  geom_text_repel(aes(x = 15, y = 3), 
                  data = numCandGens[numCandGens$cutOff == 15 & numCandGens$genesCount == 3,], 
                  label ="MEOX\nCCL21", inherit.aes = F, color = "cornflowerblue",nudge_y = 2.1) + 
  geom_text_repel(aes(x = 14, y = 4), 
                  data = numCandGens[numCandGens$cutOff == 14 & numCandGens$genesCount == 4,], 
                  label ="CLEC10A", inherit.aes = F, color = "cornflowerblue",nudge_y = 1.5) + 
  geom_text_repel(aes(x = 14, y = 2), 
                  data = numCandGens[numCandGens$cutOff == 14 & numCandGens$genesCount == 2,], 
                  label ="ITHI2\nGAP43", inherit.aes = F, color = "brown3",nudge_y = 0.5,nudge_x = -3) 


lowerPlot <- plot_grid(plotPLus,plotMinus, labels = c("B","C"), rel_widths = c(4,5) )
finalPlot <- plot_grid(plotCand, lowerPlot, labels  = c("A",""),ncol = 1)






png(filename ="Figure5.png", height = 2300, width = 3000, res = 300)
finalPlot
dev.off()

