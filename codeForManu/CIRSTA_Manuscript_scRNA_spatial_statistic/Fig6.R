# --- Fig6

## --- Fig.6a
library("Seurat")
library("dplyr")
library("ggplot2")
library(patchwork)
library(viridis)
data <- readRDS("Stereo-seq-DDCliver-22slides-integrated.rds")
gene_1=read.csv("expression_Allspatialgenes_eachzone_expressionratio0.001.csv",row.names=1) # this table calculate the average expression of each genes in each layer at different time points, and keep the genes whose expressed ratio > 1бы bins in at least one timepoint.

info=t(gene_1["Mki67",])
info=as.data.frame(info)
colnames(info)="Expression"
info$time=rep(c("D0","D8","D17","R2","R7","R21"),each = 9)
info$time=factor(info$time,levels=c("D0","D8","D17","R2","R7","R21"))
info$layer=rep(c("1","2","3","4","5","6","7","8","9"),6)
info$layer = as.numeric(info$layer)
info
theme_set(theme_bw(base_size = 20))
options(repr.plot.width = 5.5, repr.plot.height = 6)
p <- ggplot(info, aes(x = layer, y = Expression, color = time))+
  ylim(0,0.28)+
  annotate("rect",xmin=1,xmax=2,ymin=0.198,ymax=0.215,alpha = .2)+ 
  annotate("rect",xmin=3,xmax=6,ymin=0.170,ymax=0.187,alpha = .2)+
  annotate("rect",xmin=7,xmax=9,ymin=0.158,ymax=0.164,alpha = .2)+
  geom_line(size = 1) +
  scale_x_continuous(breaks = c(1:9)) + labs(title = "Mki67") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+  theme(legend.title = element_blank(), 
                                                    legend.key=element_rect(fill='white'), 
                                                    legend.text = element_text(size=20), 
                                                    legend.key.size=unit(1,'cm') ) +  
  guides(color = guide_legend(override.aes = list(size=5)))
p2= p +scale_color_npg()
p2

# calculate p-value
data_R2=subset(data,subset=time == "R2")
data_R2$region='0'
data_R2$region[data_R2$rank %in% c("1","2")]="central"
data_R2$region[data_R2$rank %in% c("3","4","5","6")]="mid"
data_R2$region[data_R2$rank %in% c("7","8","9")]="portal"
table(data_R2$region)
library(ggpubr)
library(patchwork)
set.seed(123)

df_R2=data_R2@meta.data[,c("region"),drop=F]
expression=data_R2@assays$RNA@data["Mki67",row.names(data_R2@meta.data)]
df_R2_2=merge(df_R2,expression,by=0)
row.names(df_R2_2)=df_R2_2[,1]
df_R2_2=df_R2_2[,-1]
colnames(df_R2_2)=c("region","expression")
df_R2_2$expression[df_R2_2$expression> quantile(df_R2_2$expression,0.995)[[1]]]=quantile(df_R2_2$expression,0.995)[[1]]
df_R2_2$expression[df_R2_2$expression< quantile(df_R2_2$expression,0.005)[[1]]]=quantile(df_R2_2$expression,0.005)[[1]]

my_comparisons=list(c("mid","central"),c("mid","portal"),c("central","portal"))
df_R2_2$expression=as.numeric(df_R2_2$expression)
options(repr.plot.width = 4.5, repr.plot.height = 6)
p1=ggboxplot(df_R2_2, x = "region", y = "expression",
             bxp.errorbar=T,
             width = 0.5,
             color = "region", 
             palette="npg")+ stat_compare_means(comparisons=my_comparisons)+ggtitle("R2 boxplot pvalue_69zone")
p1




## --- Fig.6b, environment in Fig.6a.
PlotLR=function(ligand,receptor){
  Ligandinfo=gene_1[ligand,]
  Receptorinfo=gene_1[receptor,]
  a=t(Ligandinfo)
  b=t(Receptorinfo)
  form<-data.frame(cbind(a,b))
  form$c <- form[,1] * form[,2]
  d=t(form[,3])
  colnames(d)=rownames(form)
  rownames(d)=paste(ligand,"_",receptor)
  return(d)
}

LRlist=c('Postn_Egfr','Igf2_Igf2r','Sfrp1_Fzd2','Wnt5a_Ror2','Tgfb2_Tgfbr1')
a=1
for (i in LRlist){ 
  if (a==1){
    LR=strsplit(i, '_')
    infoLRall=PlotLR(LR[[1]][1],LR[[1]][2])
    a=2
  }else{
    LR=strsplit(i, '_')
    infoLR=PlotLR(LR[[1]][1],LR[[1]][2])
    infoLRall=rbind(infoLRall,infoLR)
  }
}
infoLRall=na.omit(infoLRall)


x=1
for (i in row.names(infoLRall)){ 
  LRscore = as.data.frame(matrix(nrow=54,ncol=0))
  LRscore$score=infoLRall[i,]
  row.names(LRscore)=colnames(infoLRall)
  LRscore$time=rep(c("D0","D8","D17","R2","R7","R21"),each = 9)
  LRscore$time=factor(LRscore$time,levels=c("D0","D8","D17","R2","R7","R21"))
  LRscore$layer=rep(c("1","2","3","4","5","6","7","8","9"),6)
  LRscore$layer = as.numeric(LRscore$layer)
  theme_set(theme_bw(base_size = 20))
  options(repr.plot.width = 6, repr.plot.height = 5)
  p <- ggplot(LRscore, aes(x = layer, y = score, color = time)) + geom_line(size = 1) +
    scale_x_continuous(breaks = c(1:9)) + labs(title = i) + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(plot.margin=unit(rep(1,4),'cm'))+theme(aspect.ratio=1)+scale_color_npg()
  plot(p)
  assign(paste("p",sep = "",x),p)
  x=x+1
}
print(x)



## --- Fig.6d
library("Seurat")
library("dplyr")
library("ggplot2")
library("ggpubr")
data <- readRDS("./allStage.scRNAseq.expression.rds") 
Idents(data)="time"
R2=subset(data,idents= "R2")
Idents(R2)="anno_Hepzone"
gene=c("Tgfb2")
options(repr.plot.width=10.5, repr.plot.height=5)
p=DotPlot(R2,features=rev(gene),dot.scale=11)+scale_y_discrete("")+scale_x_discrete("")+coord_flip()
p2=p+theme(panel.grid=element_blank())+theme_bw(base_size = 18)+RotatedAxis()+
  scale_color_viridis_c()+labs(x=NULL,y=NULL)+theme(aspect.ratio=0.1)
plot(p2)



## --- Fig.6f, environment in Fig.6a.
library(ggsci)
info=as.matrix(gene_1)
LRlist=c("Atoh8")
for (i in LRlist){ 
  #print(i)
  LRscore = as.data.frame(matrix(nrow=54,ncol=0))
  LRscore$Expression=info[i,]
  row.names(LRscore)=colnames(info)
  LRscore$time=rep(c("D0","D8","D17","R2","R7","R21"),each = 9)
  LRscore$time=factor(LRscore$time,levels=c("D0","D8","D17","R2","R7","R21"))
  LRscore$layer=rep(c("1","2","3","4","5","6","7","8","9"),6)
  #print(LRscore)
  LRscore$layer = as.numeric(LRscore$layer)
  theme_set(theme_bw(base_size = 20))
  options(repr.plot.width = 6, repr.plot.height = 4.8)
  p <- ggplot(LRscore, aes(x = layer, y = Expression, color = time)) + geom_line(size = 1) +
    scale_x_continuous(breaks = c(1:9)) + labs(title = i) + theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(plot.margin=unit(rep(1,4),'cm'))+scale_color_npg()+  theme(legend.title = element_blank(), 
                                                                     legend.key=element_rect(fill='white'), 
                                                                     legend.text = element_text(size=20),
                                                                     legend.key.size=unit(1,'cm') ) +  
    guides(color = guide_legend(override.aes = list(size=2.5)))
  plot(p)
  
}

options(repr.plot.width = 6, repr.plot.height = 4.8)
p2=p+ylim(0.13, 0.8)+
  annotate("rect",xmin=1,xmax=2,ymin=0.611,ymax=0.632,alpha = .2)+
  annotate("rect",xmin=7,xmax=9,ymin=0.476,ymax=0.516,alpha = .2)+
  annotate("rect",xmin=1,xmax=2,ymin=0.345,ymax=0.370,alpha = .2)+
  annotate("rect",xmin=7,xmax=9,ymin=0.326,ymax=0.346,alpha = .2)
p2

# calculate p-value
data_D17=subset(data,subset=time== "D17")
data_D17$region='0'
data_D17$region[data_D17$rank %in% c("1","2")]="central"
data_D17$region[data_D17$rank %in% c("3","4","5","6")]="mid"
data_D17$region[data_D17$rank %in% c("7","8","9")]="portal"
df_D17=data_D17@meta.data[,c("region"),drop=F]
expression=data_D17@assays$RNA@data["Atoh8",row.names(data_D17@meta.data)]
df_D17_2=merge(df_D17,expression,by=0)
row.names(df_D17_2)=df_D17_2[,1]
df_D17_2=df_D17_2[,-1]
colnames(df_D17_2)=c("region","expression")
df_D17_2$expression[df_D17_2$expression> quantile(df_D17_2$expression,0.995)[[1]]]=quantile(df_D17_2$expression,0.995)[[1]]
df_D17_2$expression[df_D17_2$expression< quantile(df_D17_2$expression,0.005)[[1]]]=quantile(df_D17_2$expression,0.005)[[1]]
my_comparisons=list(c("mid","central"),c("mid","portal"),c("central","portal"))
df_D17_2$expression=as.numeric(df_D17_2$expression)
options(repr.plot.width = 4.5, repr.plot.height = 6)
p1=ggboxplot(df_D17_2, x = "region", y = "expression",
             bxp.errorbar=T,
             width = 0.5,
             color = "region", 
             palette="npg")+ stat_compare_means(comparisons=my_comparisons)+ggtitle("D17 Atoh8 pvalue")
p1

data_D8=subset(data,subset=time == "D8")
data_D8$region='0'
data_D8$region[data_D8$rank %in% c("1","2")]="central"
data_D8$region[data_D8$rank %in% c("3","4","5","6")]="mid"
data_D8$region[data_D8$rank %in% c("7","8","9")]="portal"
df_D8=data_D8@meta.data[,c("region"),drop=F]
expression=data_D8@assays$RNA@data["Atoh8",row.names(data_D8@meta.data)]
df_D8_2=merge(df_D8,expression,by=0)
row.names(df_D8_2)=df_D8_2[,1]
df_D8_2=df_D8_2[,-1]
colnames(df_D8_2)=c("region","expression")
df_D8_2$expression[df_D8_2$expression> quantile(df_D8_2$expression,0.995)[[1]]]=quantile(df_D8_2$expression,0.995)[[1]]
df_D8_2$expression[df_D8_2$expression< quantile(df_D8_2$expression,0.005)[[1]]]=quantile(df_D8_2$expression,0.005)[[1]]
my_comparisons=list(c("mid","central"),c("mid","portal"),c("central","portal"))
df_D8_2$expression=as.numeric(df_D8_2$expression)
library(ggpubr)
library(patchwork)
set.seed(123)
options(repr.plot.width = 4.5, repr.plot.height = 6)
p2=ggboxplot(df_D8_2, x = "region", y = "expression",
             bxp.errorbar=T,
             width = 0.5,
             color = "region", 
             palette="npg")+ stat_compare_means(comparisons=my_comparisons)+ggtitle("D8 Atoh8 pvalue")
p2



## --- Fig.6h, this figure was plotted using GraphPad.

