library(ReactomePA)
library("org.Mm.eg.db")
library("DOSE")
library("clusterProfiler")
library(ggraph)
library(ggplot2)
library(edgeR)
library(limma)
library(AnnotationDbi)
library(IRanges)
library(S4Vectors)
library(Biobase)
library(BiocGenerics)
library(parallel)
library(stats4)
setwd("D:/")
mouse.mat<-read.csv(file="heart RNAseq.csv",header = TRUE,row.names = 1)
mouse.mat<-mouse.mat[!duplicated(mouse.mat$gene_symbol),]
row.names(mouse.mat)<-unique(mouse.mat$gene_symbol)
mouse.mat<-mouse.mat[,3:26]
name<-colnames(mouse.mat)
hc=hclust(dist(t(mouse.mat)))
plot(hc)
name<-colnames(mouse.mat)
ID<-rep(c("P7hLPS4h heart","P7hLP14h heart","P7h heart","control heart"),each=6)
meta.data1<-cbind(name,ID)
colnames(mouse.mat)=rownames(meta.data1)
meta.data1<-as.data.frame(meta.data1)
sum(mouse.mat$exp_P7hLP14h heart4)
#Quality control
boxplot(mouse.mat,ylim=c(0,50))
mouse.mat["Gapdh",]
#PCA and Heatmap
library("pheatmap")
setwd("D:/")
mouse.mat<-read.csv(file="heart RNAseq.csv",header = TRUE,row.names = 1)
mouse.mat<-mouse.mat[!duplicated(mouse.mat$gene_symbol),]
row.names(mouse.mat)<-unique(mouse.mat$gene_symbol)
mouse.mat<-mouse.mat[,3:26]
a<-na.omit(mouse.mat)
df=as.data.frame(t(a))
df.pca<-PCA(df,graph =F)
meta.data1<-cbind(name,ID)
meta.data1<-as.data.frame(meta.data1)
group<-meta.data1$ID
fviz_pca_ind(df.pca,
             geom.ind = "point", # show points only (nbut not "text")
             col.ind = group, # color by groups
             palette = c("#00AFBB", "#E7B800","#FF6B6B", "#4ECDC4"),
             addEllipses = TRUE, # Concentration ellipses
             legend.title = "Groups")
a<-log(edgeR::cpm(a)+1)
choose_gene=names(tail(sort(apply(a,1,mad)),18045))
choose_matrix=a[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix,show_rownames = F)
mouse.mat<-mouse.mat[,-4]
meta.data1<-meta.data1[-4,]
library("edgeR")
library("limma")
meta.data1$ID<-factor(meta.data1$ID,levels = c("P7hLPS4h heart","P7hLP14h heart","P7h heart","control heart"))
mouse.mat<-na.omit(mouse.mat)
y<-DGEList(counts=mouse.mat,group=meta.data1$ID)
mouse.mat1<-log2(mouse.mat)
group<-meta.data1$ID
design<-model.matrix(~0+factor(group))
colnames(design)<-levels(factor(group))
rownames(design)<-colnames(mouse.mat)
design
#DEG
mouse.mat1<-na.omit(mouse.mat1)
fit<-lmFit(mouse.mat1,design)
cont.matrix<-makeContrasts(contrasts = c("P7h heart-control heart","P7hLPS14h heart-P7h heart","P7hLPS4h heart-P7hLPS14h heart"),levels = design)
cont.matrix
fit2=contrasts.fit(fit,cont.matrix)
fit2<-eBayes(fit2)
tempOutput=topTable(fit2,coef = 1,n=Inf)
nrDEG=na.omit(tempOutput)
#Volcano 
DEG=nrDEG
logFC_cutoff<-with(DEG,mean(abs(logFC))+2*sd(abs(logFC)))
DEG$change=as.factor(ifelse(DEG$P.Value<0.05&abs(DEG$logFC)>logFC_cutoff,ifelse(DEG$logFC>logFC_cutoff,"UP","DOWN"),"NOT"))
this_tile<-paste0("Cutoff for logFC is",round(logFC_cutoff,3),"\nThe number of up gene is",nrow(DEG[DEG$change=="UP",]),"\nThe number of down gene is",nrow(DEG[DEG$change=="DOWN",]))
g=ggplot(data=DEG,aes(x=logFC,y=-log10(P.Value),color=change))+
  geom_point(alpha=0.4,size=0.75) +
  theme_set(theme_set(theme_bw(base_size = 10)))+
  xlab("log2 fold change")+ylab("-log10,p=adjst")+
  ggtitle(this_tile)+scale_color_manual(values = c("blue","black","red"))
print(g)
#DEG selection for GSEA
DEG.up<-DEG[which(DEG$logFC>=1.5&DEG$P.Value<0.05),]
DEG.down<-DEG[which(DEG$logFC<=1.5&DEG$P.Value<0.05),]
write.csv(DEG.up,file = "D:/DEG.up.csv",row.names = TRUE)
write.csv(DEG.down,file = "D:/DEG.down.csv",row.names = TRUE)
