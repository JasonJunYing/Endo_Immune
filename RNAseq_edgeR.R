setwd('../Intermediate Files/')

#########EdgeR#########
library(edgeR)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(pheatmap)

samp='Endo'
type='Imputed'
version='2'
filter=F
test="LRT"
ThreUSE = "Pval" #Pval/FDR
comp1="ACA"
comp2="WCA"
coef=2
toAnalysis="Group" 
model="Group"
inter="~"

#Import count matrix
counts<-read.csv(paste0(samp,'.',type,'counts',version,'.csv'),row.names = 1,stringsAsFactors = F)

# #Define group design
meta<-read.csv(paste0('meta_',samp,'.csv'),row.names = 1)
# meta$Sample = strsplit2(meta$Sample,'\\.')[,2]
# meta["Endo.L3A","Group"]='0W'
# meta$Treatment <- substr(meta$Sample,8,8)
# meta["Endo.L3A","Treatment"]='Base'
# meta$Geno <- as.factor(substr(meta$Sample,6,6))
# meta$Geno <- relevel(meta$Geno,ref="W")
# if(filter){
#   meta<-meta[(meta$Sample!="Endo.L3CA"),]
#   meta <- meta[meta$Treatment=="C",]
#   meta<-meta[meta$Geno!='L',]
#   counts<-counts[,meta$Sample]
# }
# meta$Group <- relevel(as.factor(meta$Group),ref="WCA")

if(model=="Group"){
  design<-model.matrix(~0+Group,data = meta)
  pair<-makeContrasts(paste0(model,comp1,"-",model,comp2),
                      levels = design)
}else{
  design<-model.matrix(~Group,data = meta)
}

study="Endo"
thre=50
comp=paste0(comp1,"v",comp2)

dgelist<-DGEList(counts = counts,group = meta[,toAnalysis])
keep<-filterByExpr(dgelist,group = meta[,toAnalysis])
dgelist<-dgelist[keep,,keep.lib.sizes=F]

# PCA/MDS/Cor
# pca = prcomp(t(dgelist$counts),center = T,scale. = T)
# df1 = data.frame(pca$x)
# summ1 <- summary(pca)
# xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
# ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
# p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = meta$Group))+
#   stat_ellipse(aes(fill = meta$Group),
#                type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
#   geom_point(size = 3.5)+
#   labs(x = xlab1,y = ylab1,color = "Group",title = "PCA")+
#   guides(fill = "none")+
#   theme_bw()+
#   geom_text_repel(aes(label=meta$Sample))+
#   #scale_fill_manual(values = c("purple","orange","pink"))+
#   #scale_colour_manual(values = c("purple","orange","pink"))+
#   theme(plot.title = element_text(hjust = 0.5,size = 15),
#         axis.text = element_text(size = 11),axis.title = element_text(size = 13),
#         legend.text = element_text(size = 11),legend.title = element_text(size = 13),
#         plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
# pdf(paste(study,'pca','pdf',sep = '.'),width = 7,height = 6)
# print(p.pca1)
# dev.off()
# 
# 
# pdf(paste(study,'mds','pdf',sep = '.'),width = 7,height = 6)
# colors <- rainbow(length(meta$Group),v=0.8)
# plotMDS(dgelist,var.explained = T,
#         labels = meta$Sample,
#         col = colors[as.numeric(factor(meta$Group))]
#         )
# dev.off()
# 
# pdf(paste(study,'corr_filtered','pdf',sep = '.'),width = 7,height = 6)
# pheatmap(cor(dgelist$counts))
# dev.off()

# Normalization and fit model
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
dge <- estimateGLMCommonDisp(dgelist_norm, design)
dge <- estimateGLMTrendedDisp(dge,design)
dge <- estimateGLMTagwiseDisp(dge,design)
#dge <- estimateDisp(dgelist,design) # This is NOT equal to the above
fit <- glmFit(dge, design, 
                robust = TRUE
                )

# Make contrast
if(model!="Group"){
  fitted <- glmLRT(fit,coef = coef)
}else{
  fitted <- glmLRT(fit,contrast = pair)
}

lrt <- topTags(fitted, n = nrow(dgelist$counts))
write.table(lrt, paste(study,paste0(type,version),test,comp,'txt',sep = '.'), sep = '\t', col.names = NA, quote = FALSE)

###Volcano###
library(ggplot2)
library(ggrepel)
lrt<-read.delim(paste(study,paste0(type,version),test,comp,'txt',sep = '.'),row.names = 1)
data<-lrt

  # Select DEG
thres_p = 0.05
thres_FDR = 0.1
if(ThreUSE == 'Pval'){
  up<- data[(data$logFC>=1)&
              (data$PValue<=thres_p),]
  down<- data[(data$logFC<=(-1)&
                 (data$PValue<=thres_p)),]
}else{
  up<- data[(data$logFC>=1)&
            (data$FDR<=thres_FDR),]
  down<- data[(data$logFC<=(-1)&
                 (data$FDR<=thres_FDR)),]
}
up<-up[order(up$PValue),]
down<-down[order(down$PValue),]

  # Plot
DEGs<- data.frame(Group=rep("normal",times= nrow(data)),row.names = rownames(data))
DEGs[rownames(down),"Group"]<-"down"
DEGs[rownames(up),"Group"]<-"up"
#DEGs<- factor(DEGs,levels = c("up","down","normal"))
data$symbol = rownames(data)
data$logp = (-log10(data$PValue))
data$logFDR = (-log10(data$FDR))
data$DEG = DEGs$Group

  ## MA_PLOT
ma<-ggplot(data,aes(x=logCPM,y=logFC,color=DEG))+
  scale_color_manual(values = c("up"="red","down"="green","normal"="grey"))+
  geom_point()+
  ggtitle(paste(study,type,'maplot',test,comp,version,sep = '.'))+
  xlab("logCPM")+
  ylab("logFC")
pdf(paste(study,type,'maplot',test,comp,version,'pdf',sep = '.'),height = 5,width = 6 )
print(ma)
dev.off()

data = data[data$logFDR<20,]
# lab<-data[c("Lgals3","Apoe"),]
if(ThreUSE=="Pval"){
  lab<-data[(data$PValue<=thres_p)|(rownames(data)=="Lgals3"),]
}else{
  lab<-data[data$FDR<=thres_FDR,]
}
xline=c(-2,2)
yline=ifelse(ThreUSE == 'FDR',-log(thres_FDR,10),-log(thres_p,10))

if(ThreUSE == 'Pval'){
  p<- ggplot(data=data,aes(x=logFC,y=logp,color=DEG))
}else{
  p<- ggplot(data=data,aes(x=logFC,y=logFDR,color=DEG))
}
p <- p+
  xlab("logFC")+
  ylab(ifelse(ThreUSE=="Pval","-lg(p)","-lg(FDR)"))+
  geom_point(data=data)+
  geom_text_repel(data=lab,aes(label=rownames(lab)),
                  color='black',max.overlaps = 50,size=3)+
  scale_color_manual(values = c("up"="red","down"="green","normal"="black"))+
  geom_vline(xintercept = xline, lty=2,size=I(0.2),col="grey11")+
  geom_hline(yintercept = yline, lty=2,size=I(0.2),col="grey11")+
  ggtitle(paste(comp,inter,model))+
  theme_bw()+theme(panel.background = element_rect(colour = "black",size = 1,fill = "white"),panel.grid = element_blank())
pdf(paste(study,type,'volcano',test,comp,version,'pdf',sep = '.'),height = 5,width = 6 )
print(p)
dev.off()

  # Extract counts for DEG
CPMcounts <- cpm(counts,log = T)
upcnts<-CPMcounts[rownames(up),]
downcnts<-CPMcounts[rownames(down),]
cnts<-rbind(head(upcnts,thre),
            head(downcnts,thre))

###Heatmap###
library(pheatmap)
comps1 = meta[meta[,toAnalysis]==comp1,]$Sample
comps2 = meta[meta[,toAnalysis]==comp2,]$Sample
comps = c(comps2,comps1)
cnt_matrix<-cnts[,comps] 
p<-pheatmap(cnt_matrix,
         scale = "row",
         color = colorRampPalette(c("blue", "white","red"))(100),
         annotation = meta[,c("Geno","Label")],
         #annotation_colors = annotation_col,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = TRUE, 
         height = 12,width = 6,
         filename = paste(study,type,'heatmap',test,comp,version,'pdf',sep = '.')
         ) 
p

