##### Ploting and analysis script for scRNAseq data for IgA paper
##### This analysis script expects a Seurat object that was previously created
##### Niklas Krupka, 2018-08-16

library(Seurat)
library(tidyverse)
library(stringr)   # To extract the WT vs KO identifiers in a simple way
library(cowplot)   # For the combined violin plot

##### Define parameters for tSNE plotting
NK_TSNEPlot <- function(seuratobject, ...){
  TSNEPlot(seuratobject, do.label = TRUE, pt.size = 0.1, do.return = TRUE, ...) +
    theme(aspect.ratio = 1)
}

##### Define a gene list that are interesting and should be shown in the
##### various plots
genelist <- c("Csf1r",  # Best Pan-macrophage marker in the cavities
              "Adgre1", # F4/80. high in LPM
              "Itgam",  # Cd11b, high in LPM
              "Icam2",  # CD102, high in LPM
              "Timd4",  # Some LPM express Tim4
              "H2-Aa",  # MHCII, high in SPM
              "Retnla", # RELMalpha, high in SPM
              "Ccr2",   # CCR2
              "Cd209a", # DC-SIGN  
              "Mki67", "Top2a", # Prolif
              "Map1lc3a", # Interesting gene: LC3+ Macs
              "Folr2",   # Interesting gene: Folate receptor+ Macs
              "Notch2", # Interesting gene
              "Ms4a1", "Cd79a", # BC markers
              "Fcer2a", "Fcmr", "Lifr", # Receptors and BC markers
              "Cd3e", "Cd4", "Nkg7")

##### Load Seurat object and markers (obtained by FindAllMarkers with
##### default arguments). Then calculate the Top10 markers per cluster
load(file = "xbp1.tSNE.RData")
xbp1.markers <- readRDS("./xbp1markers.rds")
xbp1.top10markers <- xbp1.markers %>% 
  group_by(cluster) %>% 
  top_n(10, avg_logFC)

##### Get the original info whether cell is from WT (-2) or KO (-1) dataset and
##### store in the metadata
cell_genotype <- str_sub(names(xbp1@ident), start = -2)
cell_genotype[cell_genotype == "-1"] <- "KO"
cell_genotype[cell_genotype == "-2"] <- "WT"
xbp1@meta.data$genotype <- cell_genotype


###########################################################################
##### ALL CALCULATIONS AND SETUP DONE AT THIS POINT. NOW IT'S ONLY PLOTTING
###########################################################################

##### Display numbers of cells of each condition in each cluster
table(xbp1@ident, xbp1@meta.data$genotype)
ggplot(data.frame(c = xbp1@ident, gt = xbp1@meta.data$genotype), 
        aes(x = c, fill = gt)) + 
  geom_bar()


##### Basic tSNE plots
xbp1 <- SetAllIdent(xbp1, id = "res.1")
NK_TSNEPlot(xbp1)
xbp1 <- SetAllIdent(xbp1, id = "genotype")
NK_TSNEPlot(xbp1) 
xbp1 <- SetAllIdent(xbp1, id = "res.1")

##### Violin plot of the indicated genes
##### This first creates a list of plots using Seurat, 
##### which is then modified by lapply and 
##### plotted in a nice grid using cowplot
plot.vln <- VlnPlot(object = xbp1, features.plot = genelist, 
                    point.size.use = 0, same.y.lims = FALSE,
                    do.return = TRUE, return.plotlist = TRUE) 
plot.vln <- lapply(plot.vln, function(x) x
                   + theme_void() 
                   + theme(plot.title = element_blank(), legend.position = "none") 
                   + coord_flip())
plot_grid(plotlist = plot.vln, nrow = 1, 
          labels = genelist, label_size = 7)


##### Feature Plot of the indicated genes
FeaturePlot(object = xbp1, features.plot = c("Csf1r", "Cd79a", "Cd3e", "Cd209a"),  
            min.cutoff = c("q6"),
            cols.use = c("lightgrey", "#b2182b"), 
            pt.size = 1, 
            no.axes = TRUE,
            nCol = 4)

##### Heatmap
DoHeatmap(object = xbp1, genes.use = xbp1.top10markers$gene, 
          slim.col.label = TRUE, remove.key = FALSE,
          group.cex = 7, cex.row = 6,
          col.low = "#2166ac", col.mid = "#f7f7f7", col.high = "#b2182b")
