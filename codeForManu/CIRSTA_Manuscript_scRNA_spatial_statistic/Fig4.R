# --- Fig4

## --- Fig.4a
library("Seurat")
library("dplyr")
library("ggplot2")
library("ggpubr")
data <- readRDS("./allStage.scRNAseq.expression.rds")

data$anno_LPLCsubtype2=as.character(data$anno_LPLCsubtype2)
data$anno_LPLCsubtype2[data$anno_LPLCsubtype2 %in% c("LAM","cDC","KC","Monocyte","LSEC","PC_LVEC","PP_LVEC","HSC","Effector T","B cell","T cell","Mesothelial cell","pDC","NK","ILC1","Fibroblast","Ery","VSMC","Neutrophil","LyEC")]="gray"
table(data$anno_LPLCsubtype2)
data$anno_LPLCsubtype2=factor(data$anno_LPLCsubtype2,levels=c("Hep","Proliferating hep","LPLC_1","LPLC_2","Cholangiocyte","gray"))

theme_set(theme_bw(base_size = 15))
col=c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF','#F39B7FFF','gray')
options(repr.plot.width=7, repr.plot.height=6)
p <- DimPlot(data, reduction = "umap",label = F, repel = T,group.by = "anno_LPLCsubtype2",cols=col)+ggtitle("Cell type")
p



## --- Fig.4b, environment in Fig.4a.
library(ggsci)
Idents(data)="anno_LPLCsubtype"
HepChol=subset(data,idents= c("Hep","Proliferating hep","LPLC_1","LPLC_2","Cholangiocyte"))
HepChol$anno_LPLCsubtype=factor(HepChol$anno_LPLCsubtype,levels=rev(c("Hep","Proliferating hep","LPLC_1","LPLC_2","Cholangiocyte")))

options(repr.plot.width=6, repr.plot.height=3)
col=rev(c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF','#F39B7FFF'))
p=VlnPlot(HepChol, features = c("Apoc2","Hnf4a","Spp1","Sox9","Epcam","Hnf1b","Krt19","Krt7"),stack = T,
          group.by="anno_LPLCsubtype",log = TRUE,pt.size=0,fill.by='ident',col=col)+xlab("") + ylab("") + ggtitle("") +scale_color_npg()
p



## --- Fig.4d
library("Seurat")
library("dplyr")
library("ggplot2")
data <- readRDS("./allStage.scRNAseq.expression.rds")
Idents(data)="anno_LPLCsubtype"
HepChol=subset(data,idents= c("Hep","Proliferating hep","LPLC_1","LPLC_2","Cholangiocyte"))

DefaultAssay(HepChol) <- "RNA"
#FindAllmarkers
annotated.markers <- FindAllMarkers(data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
annotated.markers$gene<-rownames(annotated.markers)
annotated.markers =filter(annotated.markers ,p_val_adj<=0.05)
annotated.markers2=arrange(annotated.markers ,cluster,desc(avg_log2FC))
annotated.markers2=filter(annotated.markers2,cluster %in% c("Hep","LPLC_1","LPLC_2","Cholangiocyte")) # Focusing on Hepatocyte differentiation.


library(monocle)
set.seed(42)
subobj <- subset(HepChol, downsample = 22000) #downsampling
table(subobj@active.ident)
table(subobj@active.ident) %>% sum

data2 <- as(as.matrix(subobj@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = subobj@meta.data)

fData <- data.frame(gene_short_name = row.names(data2), row.names = row.names(data2))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data2,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds
pData(monocle_cds)
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

express_genes <- unique(annotated.markers2$gene)
express_genes
express_genes %>% length

monocle_cds <- setOrderingFilter(monocle_cds, express_genes)
plot_ordering_genes(monocle_cds)
plot_pc_variance_explained(monocle_cds, return_all = F) # norm_method='log'
options(repr.plot.width = 5, repr.plot.height = 12)
monocle_cds <- reduceDimension(monocle_cds, num_dim = 25, reduction_method = 'tSNE')
monocle_cds <- clusterCells(monocle_cds, method = "louvain", num_clusters = 10)
plot_cell_clusters(monocle_cds, cell_size = 0.5) +
  labs(x = "tSNE1", y = "tSNE2")

monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                               method = 'DDRTree')

monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, root_state = 1) # take one of them as root.
p2=plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 1) 
p3=plot_cell_trajectory(monocle_cds, color_by = "anno_LPLCsubtype", cell_size = 1) 
options(repr.plot.width = 10, repr.plot.height = 5)
p2+p3
# because it's reported LPLC was reprogramming from hepatocytes during DDC injury, thus we set state 2 (hepatocyte population) as root.
monocle_cds_2 <- orderCells(monocle_cds, root_state = 2)
p1=plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 1) +
  scale_color_viridis_c()
p2=plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 1) 
p3=plot_cell_trajectory(monocle_cds, color_by = "anno_LPLCsubtype", cell_size = 1) 
options(repr.plot.width = 15, repr.plot.height = 5)
p1+p2+p3

# export the density of each population across pseudotime (differentiation into Cholangiocyte).
monometa=pData(monocle_cds)
ToChol=filter(monometa,State %in% c("3","1"))
library(ggsci)
library(ggplot2)
library(ggridges)
options(repr.plot.width = 5.5, repr.plot.height = 3.5)
p=ggplot(data=ToChol, aes(x=Pseudotime, y=anno_LPLCsubtype)) +
  geom_density_ridges_gradient(aes(fill = `anno_LPLCsubtype`), scale =12, size = 0.3) +  theme(legend.position = "none")+
  theme_bw(base_size = 16)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+ggtitle("Trajectory to cholangioycte")+
  scale_fill_manual(values=c('#E64B35FF','#4DBBD5FF','#00A087FF','#3C5488FF','#F39B7FFF'))
p




##--- Fig.4g
library("Seurat")
library("dplyr")
library("ggplot2")
library(patchwork)
# donwsampling
data <- readRDS("../DDC_ST_niche_5slides.rds") # using the dataset including representative sections (DY1_D0, FS1_D8, FS4_D17, FS5_R2 and FR8_R7).
table(data$domain)
data2=subset(data,idents=c("Hepa-domain","LPLC-domain1","LPLC-domain2","Chol-domain"))
set.seed(42)
subobj <- subset(data2, downsample = 8000) #downsampling

library(monocle)
data2 <- as(as.matrix(subobj@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = subobj@meta.data)
fData <- data.frame(gene_short_name = row.names(data2), row.names = row.names(data2))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
monocle_cds <- newCellDataSet(data2,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds
pData(monocle_cds)
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)

deg.cluster=read.csv("./findallmarker_annotated_domain_LPLCsubtyes.csv") #read the table in the supplementary table 3 (log2FC >=0.5), and remove the specific genes of Portal vein-area.
head(deg.cluster)
dim(deg.cluster)
express_genes <- unique(deg.cluster$gene)
express_genes
express_genes %>% length

monocle_cds <- setOrderingFilter(monocle_cds, express_genes)
plot_ordering_genes(monocle_cds)
plot_pc_variance_explained(monocle_cds, return_all = F) # norm_method='log'
options(repr.plot.width = 5, repr.plot.height = 12)
monocle_cds <- reduceDimension(monocle_cds, num_dim = 25, reduction_method = 'tSNE')
monocle_cds <- clusterCells(monocle_cds, method = "louvain", num_clusters = 10)
plot_cell_clusters(monocle_cds, cell_size = 0.5) +
  labs(x = "tSNE1", y = "tSNE2")

options(repr.plot.width = 5, repr.plot.height = 5)
p1=plot_cell_clusters(monocle_cds, cell_size = 0.5,color_by = "annotation") +
  labs(x = "tSNE1", y = "tSNE2")
p1

monocle_cds <- reduceDimension(monocle_cds, max_components = 2,
                               method = 'DDRTree')
monocle_cds <- orderCells(monocle_cds)
monocle_cds <- orderCells(monocle_cds, root_state = 1) # take state 1 (hepa-domain) as root.
p1=plot_cell_trajectory(monocle_cds, color_by = "Pseudotime", cell_size = 1) +
  scale_color_viridis_c()+theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p2=plot_cell_trajectory(monocle_cds, color_by = "State", cell_size = 1) +theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p3=plot_cell_trajectory(monocle_cds, color_by = "sub.cluster", cell_size = 1) +theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_color_manual(values=c('#9678AC','#739E30','#FFFF00','#DF7066'))
options(repr.plot.width = 18, repr.plot.height = 4)
P=p1+p2+p3
P

# export the density of each domain across pseudotime (differentiation into Cholangiocyte).
monometa=pData(monocle_cds)
monometa$sub.cluster=factor(monometa$sub.cluster,levels=c("Hepa-domain","LPLC-domain1","LPLC-domain2","Chol-domain"))

library(ggplot2)
library(ggridges)
options(repr.plot.width = 7, repr.plot.height = 3)
p1=ggplot(data=monometa, aes(x=Pseudotime, y=sub.cluster)) +
  geom_density_ridges_gradient(aes(fill = `sub.cluster`), scale = 3, size = 0.3) +  theme(legend.position = "none")+
  theme_bw(base_size = 18)+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_fill_manual(values=c('#9678AC','#739E30','#FFFF00','#DF7066'))+ggtitle("domain trajectory")
p1




## --- Fig.4h
library("Seurat")
library("dplyr")
library("ggplot2")
library(patchwork)
library(viridis)
data <- readRDS("Niche_22_slides_data_annotation_advanced_ed2.rds") # "Niche_22_slides_data_annotation_advanced_ed2.rds" this rds file contained domain annotation across 22 slides.
data_select2=subset(data,subset=annotation=="Portal vein-area",invert=T)
data_select2$annotation=factor(data_select2$annotation,levels=rev(c("Hepa-domain","LPLC-domain1","LPLC-domain2","Chol-domain")))

gene=c('Spink1','Mmp2','Mmp7','Timp1','Fbln1','Fbn1','Postn','Pf4','Sfrp1','S100a4','S100a9','Serpine2','Hspg2','Gas6','Cd24a')
options(repr.plot.width=20, repr.plot.height=4)
p1=VlnPlot(data_select2, features =gene,group.by="annotation",stack = T,
           log = TRUE,pt.size=0,fill.by='ident',cols=rev(c("#c77cff","#7CAE00","#ffff00","#fa776e")))+xlab("") + ylab("") + ggtitle("LPLC domain subtype specific ligands") 
p1
p2=p3+theme(legend.position = "none",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(),
            axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5))
p2




## --- Fig.4j, environment in Fig.4h.
ref.TGFbwiki=list(c("Bambi","Bmp4","Crebbp","Ctnnb1","Egf","Eng","Ep300","Fkbp1a","Fos","Foxh1","Fst","Hras1","Ifng","Inhba","Itgb6","Jak1","Jun","Lef1","Lif","Ltbp1","Mapk3","Mapk9","Nfkb1","Nog","Runx2","Runx3","Serpine1","Ski","Skil","Smad1","Smad2","Smad3","Smad4","Smad5","Smad6","Smad7","Smad9","Spp1","Stat1","Stat3","Tcfe3","Tgfb1","Tgfbr1","Tgfbr2","Tgfbr3","Tgif1","Thbs1","Tnf","Wnt1","Zeb2","Zfp423","Zfyve9"))
ref.Notch = list(c("Jag1","Jag2","Hey1","Hey2","Heyl","Hes1","Notch1","Notch2","Hes5"))
data <- AddModuleScore(object = data, features = get("ref.Notch"), name = gsub("ref.", "", "ref.Notch"))
data <- AddModuleScore(object = data, features = get("ref.TGFbwiki"), name = gsub("ref.", "", "ref.TGFbwiki"))

table(data$annotation)
# draw by niche
featureModule_niche <- function(timelist, feature_list) {
  meta <- data@meta.data
  meta <- meta[meta$time == timelist, c("annotation", feature_list)]
  
  res.df <- matrix(rep(0, 5*length(feature_list)), nrow = 5) %>% as.data.frame()
  rownames(res.df) <- c("Hepa-domain","LPLC-domain1","LPLC-domain2","Chol-domain","Portal vein-area")
  colnames(res.df) <- feature_list
  for (i in c("Hepa-domain","LPLC-domain1","LPLC-domain2","Chol-domain","Portal vein-area")) {
    meta.layer <- meta[meta$annotation == i,]
    res.df[i,] <- apply(meta.layer[,feature_list, drop=F], 2, mean)
  }
  res.df <- t(res.df) %>% as.data.frame()
  colnames(res.df) <- paste0(timelist, "_", colnames(res.df))
  rownames(res.df) <- names(feature_list)
  return(res.df)
}
pathways <- c('Notch1','TGFbwiki1')
names(pathways) <- pathways

for (i in c("D0","D8","D17")) {
  if (i == "D0") {
    pathways.df <- featureModule_niche(i, pathways)
  } else {
    t.df <- featureModule_niche(i, pathways)
    pathways.df <- cbind(pathways.df, t.df)
  }
}; rm(t.df)

info=pathways.df[2,c(1,4,6:9,11:14)] 
info
row.names(info)="TGFb signaling"
options(repr.plot.width = 6, repr.plot.height = 6)
ann_colors = list(time = c(D0="#440154FF",D8= "#414487FF",D17="#2A788EFF")) 
annotation_col=data.frame(time = factor(c(rep("D0",2),rep("D8",4),rep("D17",4)),
                                        levels=c("D0","D8","D17")), row.names = colnames(info))
p=pheatmap::pheatmap(info, show_colnames = T, scale = "row",
                     gaps_col = c(2,6),
                     labels_col = c("Hepa-domain","Chol-domain",rep(c("Hepa-domain","LPLC-domain1","LPLC-domain2","Chol-domain"),2)),  
                     cluster_rows =F,annotation_colors=ann_colors,
                     cluster_cols = F, annotation_col = annotation_col, color = colorRampPalette(c("navy",
                                                                                                   "white", "firebrick3"))(200),cellheight = 16,cellwidth = 20,fontsize = 10,main="Tgfb signaling in different domain")
p



## --- Fig.4k, environment in Fig.4h.
data_notvein=subset(data,subset=annotation=="Portal vein-area",invert=T)
gene=c("Tgfb2")
options(repr.plot.width=6, repr.plot.height=4)
p=DotPlot(data_notvein,features=gene,dot.scale=11,cols = c("Snow", "#228B22"),group.by="annotation",scale.max=15)+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+labs(x=NULL,y=NULL)+ggtitle("Tgfb2 in different domain")+theme(aspect.ratio=0.35)
plot(p2)