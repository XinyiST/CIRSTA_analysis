# --- Fig2

## --- Fig.2a
#"HeatmapSelectedGO.csv",included in the source data, was exported from metascape using the overlap of D17vsD0 DEG and Domain DEG in Stereo-seq (Supplementary table 3).
cts<- read.csv("./HeatmapSelectedGO.csv",row.names = 2)[,-1]
head(cts)
colnames(cts)=c("Hepa-domain","LPLC-domain","Chol-domain")
library(pheatmap)
b<-pheatmap(cts,scale = "row",show_rownames = T, main=paste0("domain heatmap"),cellwidth = 20, cellheight = 10,cluster_col = F,cluster_row =F, 
            color = colorRampPalette(c("navy", "white", "firebrick3"))(50),treeheight_col=10,treeheight_row=10,fontsize = 8)





## --- Fig.2b
library("Seurat")
library("dplyr")
library("ggplot2")
library(patchwork)
library(viridis)
data <- readRDS("Niche_22_slides_data_annotation_advanced_ed2.rds") # "Niche_22_slides_data_annotation_advanced_ed2.rds" this rds file contained domain annotation across 22 slides.

# module score calculation
ref.senescenceGO=list(c("Abi3","Abl1","Akt3","Arg2","Arntl","B2m","Bcl2l12","Bcl6","Bmpr1a","Brca2","Calr","Cav1","Cdk6","Cdkn1a","Cdkn2a","Cdkn2b","Cgas","Dnaja3","Ecrg4","Eef1e1","Fbxo5","Fzr1","Hmga1","Hmga1b","Hmga2","Hras","Htt","Icmt","Id2","Ilk","Ing2","Kat6a","Kras","Lmna","Map2k1","Map3k3","Mapk14","Mapkapk5","Mif","Mir15a","Mir15b","Mir16-1","Mir16-2","Mir19b-1","Mir19b-2","Mir20a","Mir22","Mir23a","Mir23b","Mir24-1","Mir24-2","Mir25","Mir27a","Mir30a","Mir31","Mir34a","Mir34b","Mir34c","Mir99a","Mir125a","Mir125b-1","Mir125b-2","Mir128-1","Mir128-2","Mir130a","Mir132","Mir134","Mir140","Mir183","Mir188","Mir191","Mir193b","Mir199a-1","Mir199a-2","Mir211","Mir290a","Mir296","Mir342","Mir370","Mir378a","Mir378b","Mir378d","Mir483","Mir671","Mir877","Mir1897","Mir5100","Mirlet7b","Mirlet7d","Mirlet7f-1","Mirlet7f-2","Mirlet7g","Mirlet7i","Mnt","Morc3","Mtor","Nampt","Nek4","Npm1","Nsmce2","Nuak1","Nup62","Opa1","Pawr","Pla2r1","Plk2","Pml","Pnpt1","Pot1b","Prelp","Prkcd","Prkdc","Prmt6","Pten","Rbl1","Rsl1d1","Sirt1","Sirt6","Slc30a10","Smc5","Smc6","Spi1","Srf","Suv39h1","Tbx2","Tbx3","Terf2","Tert","Trp53","Trp63","Twist1","Vash1","Wnt1","Wnt16","Wrn","Ybx1","Ypel3","Zfp217","Zfp277","Zmiz1","Zmpste24"))
ref.DAMPpaper=list(c("Bgn","Trem2","Trem1","Tlr3","Tlr2","Nlrp3","Cgas","Tlr4","Tlr9","Tlr7","Fpr1","Casr","Fpr2","Fpr3","Gprc6a","Aim2","Trpm2","Clec9a","Clec4e","Clec7a","Rigi","Ifih1","Ager","Fpr1","P2ry2","P2ry6","P2ry12","P2rx7"))
ref.ECM=list(c("2300002M23Rik","Abi3bp","Abl1","Abl2","Acan","Adamts1","Adamts2","Adamts3","Adamts4","Adamts5","Adamts6","Adamts7","Adamts8","Adamts9","Adamts10","Adamts12","Adamts13","Adamts14","Adamts15","Adamts16","Adamts17","Adamts18","Adamts19","Adamts20","Adamtsl1","Adamtsl2","Adamtsl3","Adamtsl4","Adtrp","Aebp1","Agt","Ambn","Angptl7","Antxr1","Anxa2","Apbb1","Apbb2","Aplp1","Aplp2","App","Atp7a","Atxn1l","Axin2","B4galt1","Bcl3","Bmp2","Carmil2","Cav1","Cav2","Ccdc80","Ccn1","Ccn2","Cflar","Chadl","Clasp1","Clasp2","Cma1","Col1a1","Col1a2","Col2a1","Col3a1","Col4a1","Col4a2","Col4a3","Col4a4","Col4a5","Col4a6","Col5a1","Col5a2","Col5a3","Col6a4","Col7a1","Col8a1","Col8a2","Col9a1","Col9a2","Col9a3","Col10a1","Col11a1","Col11a2","Col13a1","Col14a1","Col15a1","Col16a1","Col17a1","Col18a1","Col19a1","Col22a1","Col23a1","Col24a1","Col27a1","Col28a1","Colgalt1","Colq","Comp","Cpb2","Creb3l1","Crispld2","Crtap","Csgalnact1","Cst3","Ctsg","Ctss","Cyp1b1","Cyp2j6","Dag1","Ddr1","Ddr2","Dmp1","Dnajb6","Dpp4","Dpt","Dspp","Ecm2","Efemp2","Egfl6","Egflam","Elf3","Eln","Emilin1","Eng","Ercc2","Ero1a","Ero1b","Ets1","Exoc8","Ext1","Fap","Fbln1","Fbln2","Fbln5","Fermt1","Fgfr4","Fkbp10","Flot1","Flrt2","Fn1","Foxc1","Foxc2","Foxf1","Foxf2","Fscn1","Fshr","Gas2","Gfap","Gfod2","Gpm6b","Grem1","Hapln2","Has1","Has2","Has3","Hmcn1","Hpn","Hpse","Hpse2","Hsd17b12","Hspg2","Ibsp","Idua","Ier3ip1","Ihh","Il6","Impg1","Impg2","Itga8","Itgb1","Itgb3","Itih1","Kazald1","Kif9","Klk4","Klk5","Lama1","Lama2","Lamb1","Lamb2","Lamb3","Lamc1","Lcp1","Lemd3","Lgals3","Lmx1b","Lox","Loxl1","Loxl2","Loxl3","Loxl4","Ltbp4","Lum","Mad2l2","Mansc4","Marco","Matn1","Matn2","Matn3","Matn4","Meltf","Mfap4","Mia","Mia3","Mkx","Mmp1b","Mmp2","Mmp3","Mmp7","Mmp8","Mmp9","Mmp10","Mmp11","Mmp12","Mmp13","Mmp14","Mmp15","Mmp16","Mmp17","Mmp19","Mmp20","Mmp21","Mmp23","Mmp24","Mmp25","Mmp27","Mmp28","Mpv17","Mpzl3","Myf5","Myh11","Myo1e","Ndnf","Nepn","Nf1","Nfkb2","Nid1","Notch1","Nox1","Noxo1","Nphp3","Npnt","Nr2e1","Ntn4","Olfml2a","Olfml2b","Optc","Otol1","P3h1","P3h4","P4ha1","Papln","Pbxip1","Pdgfra","Pdpn","Phldb1","Phldb2","Plg","Plod3","Pmp22","Pomt1","Postn","Prdm5","Prdx4","Prickle1","Ptk2","Ptx3","Pxdn","Qsox1","Ramp2","Rb1","Reck","Rgcc","Ric1","Ric8a","Rxfp1","Scara3","Scx","Sema5a","Serac1","Serpinb5","Serpinf2","Serpinh1","Sfrp2","Sh3pxd2b","Slc2a10","Slc39a8","Smad3","Smarca4","Smoc1","Smoc2","Smpd3","Sox9","Spint1","Spint2","Spock2","Sulf1","Sulf2","Tcf15","Tfap2a","Tfip11","Tgfb1","Tgfb2","Tgfbi","Tgfbr1","Thsd4","Tie1","Tmem38b","Tnf","Tnfrsf1a","Tnfrsf1b","Tnfrsf11b","Tnr","Tnxb","Vhl","Vipas39","Vit","Vps33b","Vtn","Vwa1","Wdr72","Wnt3a","Wt1","Zfp469"))
ref.inflacytochem=list(c("Vcam1","Cxcl14","Ccl6","Ccl9","Ccl2","Ccl5","Cx3cl1","Cxcl1","Cxcl2","Cxcl5","Il1a","Il1b","Il1r1","Il18","Il2","Il4","Il7","Il9","Il13","Il15","Il3","Il5","Csf2","Il6","Il11","Csf3","Lif","Osm","Il10","Il20","Il16","Ifng","Ltb","Tnfsf9","Tnfsf13","Cd70","Tnfsf18","Tnfsf14","Tnfsf4","Tnfsf13b","Tnfsf10","Tnfsf12","Tnfsf11","Epo","Tpo","Kitl","Csf1r","Mst1"))
ref.DNAdamage=list(c("Abl1","Akt1","Apaf1","Atm","Atr","Atrip","Bax","Bbc3","Bid","Brca1","Casp3","Casp8","Casp9","Ccnb1","Ccnb2","Ccnb3","Ccnd1","Ccnd2","Ccnd3","Ccne1","Ccne2","Cdc25a","Cdc25c","Cdk2","Cdk4","Cdk5","Cdk6","Cdkn1a","Cdkn1b","Chek1","Chek2","Gm10053","Cycs","Ddb2","E2f1","Fancd2","Fas","Gadd45a","Gadd45b","Gadd45g","Hus1","Mdm2","Mre11a","Nbn","Pmaip1","Pml","Prkdc","Rad1","Rad17","Rad50","Rad51","Rad52","Rad9a","Rb1","Rpa2","Rrm2b","Sesn1","Sfn","Smc1a","Tlk1","Tlk2","Tnfrsf26","Tnfrsf22","Tnfrsf23","Trp53","Cdk1","Creb1","H2ax","Pidd1","Rfc1","Myc"))
ref.Oxidate=list(c("Nfe2l2","Gsr","Txnrd1","Gclc","Mt1","Mt2","Fth1","Hmox1","Txnip","Nqo1","Srxn1","Gsta1","Gsta2","Gstm2"))

for (i in ls(pattern = "ref.")) {
  data <- AddModuleScore(object = data, features = get(i), name = gsub("ref.", "", i))
}; rm(i)

table(data$annotation_macro)

# draw by domain
featureModule_niche <- function(timelist, feature_list) {
  meta <- data@meta.data
  meta <- meta[meta$time == timelist, c("annotation_macro", feature_list)]
  
  res.df <- matrix(rep(0, 4*length(feature_list)), nrow = 4) %>% as.data.frame()
  rownames(res.df) <- c("Hepa-domain","LPLC-domain","Chol-domain","Portal vein-area")
  colnames(res.df) <- feature_list
  for (i in c("Hepa-domain","LPLC-domain","Chol-domain","Portal vein-area")) {
    meta.layer <- meta[meta$annotation_macro == i,]
    res.df[i,] <- apply(meta.layer[,feature_list, drop=F], 2, mean)
  }
  res.df <- t(res.df) %>% as.data.frame()
  colnames(res.df) <- paste0(timelist, "_", colnames(res.df))
  rownames(res.df) <- names(feature_list)
  return(res.df)
}
pathways <- c('DAMPpaper1','Oxidate1','DNAdamage1','senescenceGO1','inflacytochem1','ECM1')
names(pathways) <- gsub("1$", "", pathways)

for (i in c("D0","D8","D17","R2","R7","R21")) {
  if (i == "D0") {
    pathways.df <- featureModule_niche(i, pathways)
  } else {
    t.df <- featureModule_niche(i, pathways)
    pathways.df <- cbind(pathways.df, t.df)
  }
}; rm(t.df)

# select the D0,D8 and D17 time points
info=pathways.df[,c(1,3,5,6,7,9,10,11,13:15,17:19,21:23)]
info

options(repr.plot.width = 6, repr.plot.height = 6)
ann_colors = list(time = c(D0="#440154FF",D8= "#414487FF",D17="#2A788EFF")) 

info2=info[,c(1:8)]
row.names(info2)=c("DAMP","Oxidative stress","DNA damage","Cellular senescence","Inflammation","ECM organization")

annotation_col=data.frame(time = factor(c(rep("D0",2),rep("D8",3),rep("D17",3)),
                                        levels=c("D0","D8","D17")), row.names = colnames(info2))
p=pheatmap::pheatmap(info2, show_colnames = T, scale = "row",
                      gaps_col = c(2,5),
                      labels_col = c("Hepa-domain","Chol-domain",rep(c("Hepa-domain","LPLC-domain","Chol-domain"),2)),  
                      cluster_rows =F,annotation_colors=ann_colors,
                      cluster_cols = F, annotation_col = annotation_col, color = colorRampPalette(c("navy",
                                                                                                    "white", "firebrick3"))(200),cellheight = 16,cellwidth = 20,fontsize = 10,main="Cellular damage signaling in different domain")
p




## --- Fig.2d
library("Seurat")
library("dplyr")
library("ggplot2")
library("ggpubr")
data <- readRDS("./allStage.scRNAseq.expression.rds")
# calculate the module score of DAMP and oxidation stress as Fig.2b.
Idents(data)="group"
#diffcellmean module score
diffcellmeanModule=function(pathways,x){
  meta = x@meta.data
  celltypes = unique(meta$anno_Hepmerge)
  res.df  = as.data.frame(matrix(nrow=length(celltypes),ncol=length(pathways)))
  rownames(res.df) <- celltypes
  colnames(res.df)=pathways
  for (type in celltypes) {
    cell.data = subset(x, subset = anno_Hepmerge == type)
    matrix = cell.data@meta.data[,pathways]
    result = apply(matrix, 2, mean)
    res.df[type,pathways] <- result
  }
  return(res.df)
}

D17=subset(data,idents='D17')
pathways_D17 <- c("DAMP1","Oxidate1")
info_D17=diffcellmeanModule(pathways_D17,D17)
table(data$anno_Hepmerge)
info2_D17=info_D17[c("Hep","LPLC","Cholangiocyte","PC_LVEC","PP_LVEC","LSEC","LyEC","Fibroblast","HSC","VSMC","Mesothelial cell","KC","LAM","Monocyte","T cell","ILC1","NK","Effector T","B cell","Neutrophil","pDC","cDC"),]

D0=subset(data,idents='D0')
pathways_D0 <- c("DAMP1","Oxidate1")
info_D0=diffcellmeanModule(pathways_D0,D0)
info2_D0=info_D0[c("Hep","LPLC","Cholangiocyte","PC_LVEC","PP_LVEC","LSEC","LyEC","Fibroblast","HSC","VSMC","Mesothelial cell","KC","LAM","Monocyte","T cell","ILC1","NK","Effector T","B cell","Neutrophil","pDC","cDC"),]

D8=subset(data,idents='D8')
pathways_D8 <- c("DAMP1","Oxidate1")
info_D8=diffcellmeanModule(pathways_D8,D8)
info2_D8=info_D8[c("Hep","LPLC","Cholangiocyte","PC_LVEC","PP_LVEC","LSEC","LyEC","Fibroblast","HSC","VSMC","Mesothelial cell","KC","LAM","Monocyte","T cell","ILC1","NK","Effector T","B cell","Neutrophil","pDC","cDC"),]

info_total=cbind(t(info2_D0),t(info2_D8),t(info2_D17))
info_total
colnames(info_total)

colnames(info_total)= c("Hep_D0","LPLC_D0","Cholangiocyte_D0","PC_LVEC_D0","PP_LVEC_D0","LSEC_D0","LyEC_D0","Fibroblast_D0","HSC_D0","VSMC_D0","Mesothelial cell_D0","KC_D0","LAM_D0","Monocyte_D0","T cell_D0","ILC1_D0","NK_D0","Effector T_D0","B cell_D0","Neutrophil_D0","pDC_D0","cDC_D0",
                        "Hep_D8","LPLC_D8","Cholangiocyte_D8","PC_LVEC_D8","PP_LVEC_D8","LSEC_D8","LyEC_D8","Fibroblast_D8","HSC_D8","VSMC_D8","Mesothelial cell_D8","KC_D8","LAM_D8","Monocyte_D8","T cell_D8","ILC1_D8","NK_D8","Effector T_D8","B cell_D8","Neutrophil_D8","pDC_D8","cDC_D8",
                        "Hep_D17","LPLC_D17","Cholangiocyte_D17","PC_LVEC_D17","PP_LVEC_D17","LSEC_D17","LyEC_D17","Fibroblast_D17","HSC_D17","VSMC_D17","Mesothelial cell_D17","KC_D17","LAM_D17","Monocyte_D17","T cell_D17","ILC1_D17","NK_D17","Effector T_D17","B cell_D17","Neutrophil_D17","pDC_D17","cDC_D17")

# sort the order by each cell type
info_total2=info_total[,c(1,22+1,44+1,2,2+22,2+44,3,3+22,3+44,4,4+22,4+44,5,5+22,5+44,
                          6,6+22,6+44,7,7+22,7+44,8,8+22,8+44,9,9+22,9+44,10,10+22,10+44,
                          11,11+22,11+44,12,12+22,12+44,13,13+22,13+44,14,14+22,14+44,20,20+22,20+44,
                          15,15+22,15+44,16,16+22,16+44,17,17+22,17+44,18,18+22,18+44,19,19+22,19+44,
                          21,21+22,21+44,22,22+22,22+44)]
info_total2

options(repr.plot.width = 25, repr.plot.height = 5)
annotation_col=data.frame(time = factor(rep(c("D0","D8","D17"), 22),
                                        levels=c("D0","D8","D17")), row.names = colnames(info_total2))
ann_colors = list(time = c(D0="#440154FF",D8= "#414487FF",D17="#2A788EFF"))  
map=pheatmap::pheatmap(info_total2, show_rownames = T, scale = "row", gaps_col = c(1:20)*3,
                       cluster_rows =F,border_color=NA,,annotation_col = annotation_col,annotation_colors=ann_colors,
                       cluster_cols = F,cellheight = 16,cellwidth = 20,fontsize = 10, color = colorRampPalette(c("navy",
                                                                                                                 "white", "firebrick3"))(200),main="Cell_stress Module scores in total")




## --- Fig.2e
library("Seurat")
library("dplyr")
library("ggplot2")
library(patchwork)
library(viridis)
data <- readRDS("Stereo-seq-DDCliver-22slides-integrated.rds")

ref.LSEC=list(c("Pecam1","Cdh5","Clec4g","Dnase1l3","Aqp1"))
ref.LyEC <- list(c("Mmrn1","Pard6g","Sned1","Col23a1","Meox1","Pdpn","Prox1","Tbx1","Slac45a3","Fxyd6"))
for (i in ls(pattern = "ref.")) {
  data <- AddModuleScore(object = data, features = get(i), name = gsub("ref.", "", i))
}; rm(i)


# draw by layers
featureModule <- function(timelist, feature_list) {
  meta <- meta #data@meta.data
  meta <- meta[meta$time == timelist, c("rank", feature_list)]
  
  res.df <- matrix(rep(0, 9*length(feature_list)), nrow = 9) %>% as.data.frame()
  rownames(res.df) <- c(1:9); colnames(res.df) <- feature_list
  for (i in c(1:9)) {
    meta.layer <- meta[meta$rank == i,]
    res.df[i,] <- apply(meta.layer[,feature_list, drop=F], 2, mean)
  }
  res.df <- t(res.df) %>% as.data.frame()
  colnames(res.df) <- paste0(timelist, "_layer", colnames(res.df))
  rownames(res.df) <- names(feature_list)
  return(res.df)
}
pathways <- c("LSEC1","LyEC1")
names(pathways) <- pathways
for (i in c("D0","D8","D17","R2","R7","R21")) {
  if (i == "D0") {
    pathways.df <- featureModule(i, pathways)
  } else {
    t.df <- featureModule(i, pathways)
    pathways.df <- cbind(pathways.df, t.df)
  }
}; rm(t.df)

pathways.dfa=pathways.df
pathways.dfa2=pathways.dfa[,c(1:27)] # select the D0,D8 and D17 time points

row.names(pathways.dfa2)=c("LSEC","LyEC")

options(repr.plot.width = 10, repr.plot.height = 5)
ann_colors = list(time = c(D0="#440154FF",D8= "#414487FF",D17="#2A788EFF"))
annotation_col=data.frame(time = factor(rep(c("D0","D8","D17"),each = 9),levels=c("D0","D8","D17")), row.names = colnames(pathways.dfa2))
map=pheatmap::pheatmap(pathways.dfa2, show_rownames = T, scale = "row", gaps_col = c(9,18),
                       labels_col = rep(c("CV", rep("", 7),"PV"), 3), annotation_colors=ann_colors,
                       cluster_rows =F,annotation_col = annotation_col, border_color=NA, 
                       cluster_cols = F,cellheight = 16,cellwidth = 8,fontsize = 10, color = colorRampPalette(c("navy",
                                                                                                                "white", "firebrick3"))(200),main="cell score in different zone")




## --- Fig.2g
library("Seurat")
library("dplyr")
library("ggplot2")
library("ggpubr")
data <- readRDS("./allStage.scRNAseq.expression.rds")

Idents(data)="group"
D17=subset(data,idents= "D17")
Idents(D17)="anno_Hepmerge"

gene=c("Mmrn1","Reln","Fxyd6","Thy1",
       "Pard6g","Col23a1",
       "Pdpn","Meox1",
       "Ehd2","Palm","Aph1b",
       "Eln",
       "Sema3a","Il7","Timp2","Lamb1")
options(repr.plot.width=10, repr.plot.height=10)
p=DotPlot(D17,features=rev(gene),dot.scale=11)+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p1=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=1)+ggtitle("LyEC DEGs at D17")
plot(p1)




## ---Fig.2h, environment in Fig.2g
Idents(D17)="anno_Hepmerge"
Immune=subset(D17,idents=c("KC","LAM","Monocyte","T cell","ILC1","NK","B cell","Neutrophil","pDC","cDC"))
feature=c("Nrp1","Plxna1")
options(repr.plot.width=8, repr.plot.height=5)
p=DotPlot(Immune,features=rev(feature),dot.scale=11,cols = c("Snow", "#228B22"))+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+labs(x=NULL,y=NULL)+ggtitle("Nrp1 and Plxna1")+theme(aspect.ratio=0.2)
plot(p2)
