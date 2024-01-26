# --- Fig.1

## ---  Fig.1g
library("Seurat")
library("dplyr")
library("ggplot2")
library("ggpubr")
data <- readRDS("./allStage.scRNAseq.expression.rds")

theme_set(theme_bw(base_size = 15))
col=c('#990F26','#CC7A88','#CC00FF','#E6B8BF','#99600F','#B3823E','#CCAA7A','#FFCC00','#E6D2B8','#FF6633','#FF6699',
       '#54990F','#78B33E','#A3CC7A','#CFE6B8','#0F8299','#B8DEE6','#7ABECC','#3E9FB3','#3D0F99','#653EB3','#967ACC','#C7B8E6','#666666','#999999')
options(repr.plot.width=8, repr.plot.height=6)
p <- DimPlot(data, reduction = "umap",label = T, repel = TRUE,group.by = "annotation",cols=col,label.size =5.5)+ggtitle("Cell type")
p