# --- Fig5

## --- Fig.5a
library("Seurat")
library("dplyr")
library("ggplot2")
library(patchwork)
library(viridis)
data <- readRDS("Stereo-seq-DDCliver-22slides-integrated.rds")

ref.Monocyte <- list(c('Ly6c2','Vcan','Cd177','Chil3','S100a4','Serpinb2','Tppp3'))
ref.LAM <- list(c('Spp1','Gpnmb','Fabp5','Cd93','Cx3cr1','Ms4a7','Mmp14','Itgax','Cd63','Cd93','Trem2','Atf3'))
ref.Neutrophil=list(c("S100a9","Retnlg","Clec4d"))
ref.Cholangiocyte <- list(c('Krt19','Krt7','Epcam'))
ref.Fibroblast <- list(c("Cd34","Col15a1","Thy1","Fbln2","Gpx3","Clec3b","Mfap4","Dpt","Entpd2"))
for (i in ls(pattern = "ref.")) {
  data <- AddModuleScore(object = data, features = get(i), name = gsub("ref.", "", i))
}; rm(i)

# draw by layers
featureModule <- function(timelist, feature_list) {
  meta <- data@meta.data
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
pathways <- c("Monocyte1","LAM1","Neutrophil1","Cholangiocyte1","Fibroblast1")
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
row.names(pathways.dfa)=c("Monocyte","LAM","Neutrophil","Cholangiocyte","Fibroblast")
options(repr.plot.width = 10, repr.plot.height = 5)
ann_colors = list(time = c(D0="#440154FF",D8= "#414487FF",D17="#2A788EFF", R2="#22A884FF",R7="#7AD151FF", R21="#FDE725FF"))
annotation_col=data.frame(time = factor(rep(c("D0","D8","D17","R2","R7","R21"),each = 9),levels=c("D0","D8","D17","R2","R7","R21")), row.names = colnames(pathways.dfa))
map=pheatmap::pheatmap(pathways.dfa, show_rownames = T, scale = "row", gaps_col = c(9,18,27,36,45), gaps_row = c(3),
                       labels_col = rep(c("CV", rep("", 7),"PV"), 6), annotation_colors=ann_colors,
                       cluster_rows =F,annotation_col = annotation_col, border_color=NA, 
                       cluster_cols = F,cellheight = 16,cellwidth = 8,fontsize = 10, color = colorRampPalette(c("navy",
                                                                                                                "white", "firebrick3"))(200),main="cell score in different zone")




## --- Fig.5c
library("Seurat")
library("dplyr")
library("ggplot2")
library(patchwork)
library(viridis)
data <- readRDS("Stereo-seq-DDCliver-22slides-integrated.rds")

ref.senescenceGO=list(c("Abi3","Abl1","Akt3","Arg2","Arntl","B2m","Bcl2l12","Bcl6","Bmpr1a","Brca2","Calr","Cav1","Cdk6","Cdkn1a","Cdkn2a","Cdkn2b","Cgas","Dnaja3","Ecrg4","Eef1e1","Fbxo5","Fzr1","Hmga1","Hmga1b","Hmga2","Hras","Htt","Icmt","Id2","Ilk","Ing2","Kat6a","Kras","Lmna","Map2k1","Map3k3","Mapk14","Mapkapk5","Mif","Mir15a","Mir15b","Mir16-1","Mir16-2","Mir19b-1","Mir19b-2","Mir20a","Mir22","Mir23a","Mir23b","Mir24-1","Mir24-2","Mir25","Mir27a","Mir30a","Mir31","Mir34a","Mir34b","Mir34c","Mir99a","Mir125a","Mir125b-1","Mir125b-2","Mir128-1","Mir128-2","Mir130a","Mir132","Mir134","Mir140","Mir183","Mir188","Mir191","Mir193b","Mir199a-1","Mir199a-2","Mir211","Mir290a","Mir296","Mir342","Mir370","Mir378a","Mir378b","Mir378d","Mir483","Mir671","Mir877","Mir1897","Mir5100","Mirlet7b","Mirlet7d","Mirlet7f-1","Mirlet7f-2","Mirlet7g","Mirlet7i","Mnt","Morc3","Mtor","Nampt","Nek4","Npm1","Nsmce2","Nuak1","Nup62","Opa1","Pawr","Pla2r1","Plk2","Pml","Pnpt1","Pot1b","Prelp","Prkcd","Prkdc","Prmt6","Pten","Rbl1","Rsl1d1","Sirt1","Sirt6","Slc30a10","Smc5","Smc6","Spi1","Srf","Suv39h1","Tbx2","Tbx3","Terf2","Tert","Trp53","Trp63","Twist1","Vash1","Wnt1","Wnt16","Wrn","Ybx1","Ypel3","Zfp217","Zfp277","Zmiz1","Zmpste24"))
ref.DAMPpaper=list(c("Bgn","Trem2","Trem1","Tlr3","Tlr2","Nlrp3","Cgas","Tlr4","Tlr9","Tlr7","Fpr1","Casr","Fpr2","Fpr3","Gprc6a","Aim2","Trpm2","Clec9a","Clec4e","Clec7a","Rigi","Ifih1","Ager","Fpr1","P2ry2","P2ry6","P2ry12","P2rx7"))
ref.ECM=list(c("2300002M23Rik","Abi3bp","Abl1","Abl2","Acan","Adamts1","Adamts2","Adamts3","Adamts4","Adamts5","Adamts6","Adamts7","Adamts8","Adamts9","Adamts10","Adamts12","Adamts13","Adamts14","Adamts15","Adamts16","Adamts17","Adamts18","Adamts19","Adamts20","Adamtsl1","Adamtsl2","Adamtsl3","Adamtsl4","Adtrp","Aebp1","Agt","Ambn","Angptl7","Antxr1","Anxa2","Apbb1","Apbb2","Aplp1","Aplp2","App","Atp7a","Atxn1l","Axin2","B4galt1","Bcl3","Bmp2","Carmil2","Cav1","Cav2","Ccdc80","Ccn1","Ccn2","Cflar","Chadl","Clasp1","Clasp2","Cma1","Col1a1","Col1a2","Col2a1","Col3a1","Col4a1","Col4a2","Col4a3","Col4a4","Col4a5","Col4a6","Col5a1","Col5a2","Col5a3","Col6a4","Col7a1","Col8a1","Col8a2","Col9a1","Col9a2","Col9a3","Col10a1","Col11a1","Col11a2","Col13a1","Col14a1","Col15a1","Col16a1","Col17a1","Col18a1","Col19a1","Col22a1","Col23a1","Col24a1","Col27a1","Col28a1","Colgalt1","Colq","Comp","Cpb2","Creb3l1","Crispld2","Crtap","Csgalnact1","Cst3","Ctsg","Ctss","Cyp1b1","Cyp2j6","Dag1","Ddr1","Ddr2","Dmp1","Dnajb6","Dpp4","Dpt","Dspp","Ecm2","Efemp2","Egfl6","Egflam","Elf3","Eln","Emilin1","Eng","Ercc2","Ero1a","Ero1b","Ets1","Exoc8","Ext1","Fap","Fbln1","Fbln2","Fbln5","Fermt1","Fgfr4","Fkbp10","Flot1","Flrt2","Fn1","Foxc1","Foxc2","Foxf1","Foxf2","Fscn1","Fshr","Gas2","Gfap","Gfod2","Gpm6b","Grem1","Hapln2","Has1","Has2","Has3","Hmcn1","Hpn","Hpse","Hpse2","Hsd17b12","Hspg2","Ibsp","Idua","Ier3ip1","Ihh","Il6","Impg1","Impg2","Itga8","Itgb1","Itgb3","Itih1","Kazald1","Kif9","Klk4","Klk5","Lama1","Lama2","Lamb1","Lamb2","Lamb3","Lamc1","Lcp1","Lemd3","Lgals3","Lmx1b","Lox","Loxl1","Loxl2","Loxl3","Loxl4","Ltbp4","Lum","Mad2l2","Mansc4","Marco","Matn1","Matn2","Matn3","Matn4","Meltf","Mfap4","Mia","Mia3","Mkx","Mmp1b","Mmp2","Mmp3","Mmp7","Mmp8","Mmp9","Mmp10","Mmp11","Mmp12","Mmp13","Mmp14","Mmp15","Mmp16","Mmp17","Mmp19","Mmp20","Mmp21","Mmp23","Mmp24","Mmp25","Mmp27","Mmp28","Mpv17","Mpzl3","Myf5","Myh11","Myo1e","Ndnf","Nepn","Nf1","Nfkb2","Nid1","Notch1","Nox1","Noxo1","Nphp3","Npnt","Nr2e1","Ntn4","Olfml2a","Olfml2b","Optc","Otol1","P3h1","P3h4","P4ha1","Papln","Pbxip1","Pdgfra","Pdpn","Phldb1","Phldb2","Plg","Plod3","Pmp22","Pomt1","Postn","Prdm5","Prdx4","Prickle1","Ptk2","Ptx3","Pxdn","Qsox1","Ramp2","Rb1","Reck","Rgcc","Ric1","Ric8a","Rxfp1","Scara3","Scx","Sema5a","Serac1","Serpinb5","Serpinf2","Serpinh1","Sfrp2","Sh3pxd2b","Slc2a10","Slc39a8","Smad3","Smarca4","Smoc1","Smoc2","Smpd3","Sox9","Spint1","Spint2","Spock2","Sulf1","Sulf2","Tcf15","Tfap2a","Tfip11","Tgfb1","Tgfb2","Tgfbi","Tgfbr1","Thsd4","Tie1","Tmem38b","Tnf","Tnfrsf1a","Tnfrsf1b","Tnfrsf11b","Tnr","Tnxb","Vhl","Vipas39","Vit","Vps33b","Vtn","Vwa1","Wdr72","Wnt3a","Wt1","Zfp469"))
ref.inflacytochem=list(c("Vcam1","Cxcl14","Ccl6","Ccl9","Ccl2","Ccl5","Cx3cl1","Cxcl1","Cxcl2","Cxcl5","Il1a","Il1b","Il1r1","Il18","Il2","Il4","Il7","Il9","Il13","Il15","Il3","Il5","Csf2","Il6","Il11","Csf3","Lif","Osm","Il10","Il20","Il16","Ifng","Ltb","Tnfsf9","Tnfsf13","Cd70","Tnfsf18","Tnfsf14","Tnfsf4","Tnfsf13b","Tnfsf10","Tnfsf12","Tnfsf11","Epo","Tpo","Kitl","Csf1r","Mst1"))
ref.DNAdamage=list(c("Abl1","Akt1","Apaf1","Atm","Atr","Atrip","Bax","Bbc3","Bid","Brca1","Casp3","Casp8","Casp9","Ccnb1","Ccnb2","Ccnb3","Ccnd1","Ccnd2","Ccnd3","Ccne1","Ccne2","Cdc25a","Cdc25c","Cdk2","Cdk4","Cdk5","Cdk6","Cdkn1a","Cdkn1b","Chek1","Chek2","Gm10053","Cycs","Ddb2","E2f1","Fancd2","Fas","Gadd45a","Gadd45b","Gadd45g","Hus1","Mdm2","Mre11a","Nbn","Pmaip1","Pml","Prkdc","Rad1","Rad17","Rad50","Rad51","Rad52","Rad9a","Rb1","Rpa2","Rrm2b","Sesn1","Sfn","Smc1a","Tlk1","Tlk2","Tnfrsf26","Tnfrsf22","Tnfrsf23","Trp53","Cdk1","Creb1","H2ax","Pidd1","Rfc1","Myc"))
ref.Oxidate=list(c("Nfe2l2","Gsr","Txnrd1","Gclc","Mt1","Mt2","Fth1","Hmox1","Txnip","Nqo1","Srxn1","Gsta1","Gsta2","Gstm2"))

for (i in ls(pattern = "ref.")) {
  data <- AddModuleScore(object = data, features = get(i), name = gsub("ref.", "", i))
}; rm(i)


# draw by layers
featureModule <- function(timelist, feature_list) {
  meta <- data@meta.data
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
pathways <- c("DAMPpaper1","Oxidate1","DNAdamage1","senescenceGO1","inflacytochem1","ECM1")
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
row.names(pathways.dfa)=c("DAMP","Oxidative_stress","DNA_damage","Cellular_senescence","Inflammation","ECM_organization")
options(repr.plot.width = 10, repr.plot.height = 5)
ann_colors = list(time = c(D0="#440154FF",D8= "#414487FF",D17="#2A788EFF", R2="#22A884FF",R7="#7AD151FF", R21="#FDE725FF"))
annotation_col=data.frame(time = factor(rep(c("D0","D8","D17","R2","R7","R21"),each = 9),levels=c("D0","D8","D17","R2","R7","R21")), row.names = colnames(pathways.dfa))
map=pheatmap::pheatmap(pathways.dfa[c(2,3,5,1,4,6),], show_rownames = T, scale = "row", gaps_col = c(9,18,27,36,45), gaps_row=c(1,3),
                       labels_col = rep(c("CV", rep("", 7),"PV"), 6), annotation_colors=ann_colors,
                       cluster_rows =F,annotation_col = annotation_col, border_color=NA, 
                       cluster_cols = F,cellheight = 16,cellwidth = 8,fontsize = 10, color = colorRampPalette(c("navy",
                                                                                                                "white", "firebrick3"))(200),main="Injury pathway in different zone")



## --- Fig.5e
library("Seurat")
library("dplyr")
library("ggplot2")
library("ggpubr")
data <- readRDS("./allStage.scRNAseq.expression.rds")
Idents(data)='annotation'
Chol=subset(data,idents=c("Cholangiocyte"))
gene=c("Ccl2","Cxcl1","Cxcl2","Csf1","Cx3cl1","Cxcl5","Cxcl16")
options(repr.plot.width=7, repr.plot.height=6)
p=DotPlot(Chol,features=rev(WntL),dot.scale=11,cols = c("Snow", "#228B22"),group.by="time")+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+labs(x=NULL,y=NULL)+ggtitle("Cholangio-hub ligands in Cholangiocyte")
plot(p2)

ref.cholangioctyeimmune <- list(c("Ccl2","Cxcl1","Cxcl2","Csf1","Cx3cl1","Cxcl5","Cxcl16"))
for (i in ls(pattern = "ref.")) {
  Chol<- AddModuleScore(object = Chol, features = get(i), name = gsub("ref.", "", i))
}; rm(i)

# calculate p-value and plot.
Cholselect=Chol@meta.data[,c('time','annotation','cholangioctyeimmune1')]
colnames(Cholselect)=c('group','annotation',"score")
df=Cholselect
df$score=as.numeric(df$score)
df$score[df$score> quantile(df$score,0.995)[[1]]]=quantile(df$score,0.995)[[1]]
df$score[df$score< quantile(df$score,0.005)[[1]]]=quantile(df$score,0.005)[[1]]
my_comparisons=list(c("D17","R2"),c("D17","R7"),c("D17","R21"),c("R2","R7"))

options(repr.plot.width = 3, repr.plot.height = 5)
p=ggboxplot(df,x ="group",y= "score",
             bxp.errorbar=T,
             width=0.5,
             color="gray",size=0.2,outlier.shape = NA)+ggtitle("Cholangio-hub immune module in Cholangiocyte difftime")+
  stat_summary(fun.y=mean, colour="red", geom="line", aes(group = 1))+stat_compare_means(comparisons = my_comparisons)
p



## --- Fig.5f, environment in Fig.5e
Mac=subset(data,idents=c("KC","LAM","Monocyte"))
feature=c("Tnfsf12")
options(repr.plot.width=6, repr.plot.height=5)
p=DotPlot(Mac,features=rev(feature),dot.scale=11,cols = c("Snow", "#228B22"),group.by="time",scale.min=1)+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+labs(x=NULL,y=NULL)+ggtitle("Tnfsf12 in Mac(LAM,KC,Mono)")+theme(aspect.ratio=0.3)
p2




## --- Fig.5g, environment in Fig.5e
gene=unique(c("Fzd1","Fzd4","Fzd5","Erbb2","Egfr","Igf1r","Igf2r","Fgfr4","Tnfrsf12a"))
options(repr.plot.width=6, repr.plot.height=6)
p=DotPlot(Chol,features=rev(gene),dot.scale=11,cols = c("Snow", "#228B22"),group.by="time",scale.min=0)+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+labs(x=NULL,y=NULL)+ggtitle("receptors in Cholangiocyte")+theme(aspect.ratio=1.8)
plot(p2)



