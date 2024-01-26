# --- Fig3

## --- Fig.3b
library(dplyr)
# read this talbe in the Supplementary Table 4. CellChat scRNA-seq D17 vs. D0 & PV upregulated interactions in selected sources and targets.
stsclr=read.csv("20230606_revision_ed3_onlyD17_PV_enriched_LRs_sc_overlap_spatial.csv",row.names = 1)

lrmerge2=stsclr
lrmerge2$receptor2=lrmerge2$receptor 
unique(lrmerge2$interaction_name_2)
unique(lrmerge2$interaction_name)
unique(lrmerge2$gs)
lrall=c(unique(lrmerge2$gs),unique(lrmerge2$gl1),unique(lrmerge2$gl2),unique(lrmerge2$gr1),unique(lrmerge2$gr2))
lrall %>% unique 
lrall %>% unique %>% length

#ligand
library(dplyr)
for (i in unique(lrmerge2$source)){
  print(i)
  celltype=filter(lrmerge2,source==i)
  unique(celltype$ligand) %>% unique %>% length %>% print
}

#receptor
for (i in unique(lrmerge2$target)){
  print(i)
  celltype=filter(lrmerge2,target==i)
  unique(celltype$receptor2) %>% unique %>% length %>% print
}




## --- Fig.3d
library("Seurat")
library("dplyr")
library("ggplot2")
library("ggpubr")
data <- readRDS("./allStage.scRNAseq.expression.rds")

Idents(data)="group"
D17=subset(data,idents= "D17")
Idents(D17)="anno_Hepmerge"

gene=unique(c("Fzd1","Fzd4","Fzd5","Lrp5","Lrp6","Met","Erbb2","Egfr","Igf1r","Igf2r","Fgfr1","Fgfr2","Fgfr3","Fgfr4","Tnfrsf12a"))
options(repr.plot.width=10, repr.plot.height=8)
p=DotPlot(D17,features=rev(gene),dot.scale=11)+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+ggtitle('Cholangiocyte receptors in D17')+theme(aspect.ratio=0.7)
plot(p2)




## --- Fig.3e, environment in Fig.3d.
Idents(D17)="anno_Hepmerge"
Immune=subset(D17,idents=c("KC","LAM","Monocyte","T cell","ILC1","NK","B cell","Neutrophil","pDC","cDC"))
gene=c("Tnfsf12")
options(repr.plot.width=10.5, repr.plot.height=5)
p=DotPlot(Immune,features=gene,dot.scale=11)+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p1=p+theme(panel.grid=element_blank())+theme_bw(base_size = 20)+RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+ggtitle('Tnfsf12 at D17')+theme(aspect.ratio=0.2)
p1




## --- Fig.3f, environment in Fig.3d.
gene=unique(c("Ccl2","Cxcl1","Cxcl2","Csf1","Cx3cl1","Cxcl2","Cxcl5","Cxcl16","Tgfb2","Pdgfb","Hbegf","Jag1"))
options(repr.plot.width=10, repr.plot.height=5)
p=DotPlot(D17,features=rev(gene),dot.scale=11)+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p1=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+ggtitle('Cholangiocyte periportal ligands in D17')+theme(aspect.ratio=0.5)
plot(p1)



## --- Fig.3i
library("Seurat")
library("dplyr")
library("ggplot2")
library(patchwork)
library(viridis)
data <- readRDS("Niche_22_slides_data_annotation_advanced_ed2.rds") # "Niche_22_slides_data_annotation_advanced_ed2.rds" this rds file contained domain annotation across 22 slides.

layerCounterByTime_niche <- function(obj, ligands) {
  require("Seurat")
  require("dplyr")
  
  # layerCounts 
  layerCounts <- function(obj.dat, pre) {
    if (missing(pre)) {
      pre <- ""
    }
    n=2
    for (i in unique(obj.dat$annotation_macro) %>%
         sort()) {
      tmp.dat <- subset(obj.dat, subset = annotation_macro == i)
      if (i == sort(unique(obj.dat$annotation_macro))[1]) {
        df <- apply(tmp.dat@assays$RNA@data, 1, mean) %>%
          as.data.frame(row.names = 1)
        colnames(df)[1] <- paste0(pre, "_domain_", i)
      } else {
        df[, n] <- apply(tmp.dat@assays$RNA@data, 1, mean)
        colnames(df)[n] <- paste0(pre, "_domain_", i)
        n=n+1
      }
    }
    return(df)
  }
  
  dat.ligands <- subset(obj, features = ligands)
  dat.ligands <- subset(dat.ligands, subset = annotation_macro %in% c("Hepa-domain","LPLC-domain","Chol-domain","Portal vein-area"))
  
  # time point
  time <- dat.ligands$time %>%
    unique() %>%
    sort()
  
  for (t in time) {
    # extract by time
    time.dat <- subset(dat.ligands, subset = time == t)
    if (t == time[1]) {
      df <- layerCounts(time.dat, t)
    } else {
      df_t <- layerCounts(time.dat, t)
      df <- cbind(df, df_t)
    }
  }
  
  anno_col <- data.frame(time = rep(time, each = 4))
  row.names(anno_col) <- colnames(df)
  return(list(df = df, anno_col = anno_col))
}

gene=c("Csf1")
info=layerCounterByTime_niche(data,gene)

# select the D0, D8 and D17 time points
info_gene=info$df[,c(1,3,5,6,7,9,10,11)]
info_gene

options(repr.plot.width = 8, repr.plot.height = 6)
ann_colors = list(time = c(D0="#440154FF",D8= "#414487FF",D17="#2A788EFF")) 
annotation_col=data.frame(time = factor(c(rep("D0",2),rep("D8",3),rep("D17",3)),
                                        levels=c("D0","D8","D17")), row.names = colnames(info_gene))
p=pheatmap::pheatmap(info_gene, show_colnames = T, scale = "row",
                      gaps_col = c(2,5),
                      labels_col = c("Hepa-domain","Chol-domain",rep(c("Hepa-domain","LPLC-domain","Chol-domain"),2)),  
                      cluster_rows =F,annotation_colors=ann_colors,
                      cluster_cols = F, annotation_col = annotation_col, color = colorRampPalette(c("navy",
                                                                                                    "white", "firebrick3"))(200),cellheight = 16,cellwidth = 20,fontsize = 10,main="Csf1 in different domain")
p
