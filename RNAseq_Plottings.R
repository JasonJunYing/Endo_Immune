setwd("D:\\Researches/Final Reports/Endo-immune Crosstalk and Atherosclerosis/Intermediate Files/")

library(stringr)
library(dplyr)
library(rstatix)
library(ggpubr)
library(ggplot2)
library(patchwork)
library(edgeR)
library(impute)
library(xlsx)

# Select genes to plot
  ##Check surface markers
surfs <- read.xlsx('D:\\Researches/Sc_common_GeneLists/Mouse_Surface_Maker_List_expanded.xlsx',sheetIndex = 1)
surfLS <- surfs$Full_name_from_nomenclature_authority
surfLS <- surfLS[surfLS!='-']

samp = 'Endo'
comp='ACAvWCA'
type='Imputed'
version='2'
thres=paste0('p_',thres_p)
compLRT = read.table(paste0(samp,'.',type,version,'.LRT.',comp,'.txt'))

data = compLRT
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

surfLS_down = intersect(surfLS,rownames(down))
surfLS_up = intersect(surfLS,rownames(up))

chklist = c(surfLS_down,surfLS_up)
#chklist = c("Cd34","Pecam1","Itga5","Itgb1","Cxcl12","Ccl5")
group = "Selected"
chknum = 8
tochk <- na.omit(chklist[1:chknum])

# Figure parameters
ncol = 3
height = max(length(tochk),3)
width = 8

# Prepare for matrix and meta
counts<-read.csv(paste0(samp,'.',type,'counts',version,'.csv'),row.names = 1,stringsAsFactors = F)
meta<-read.csv(paste0('meta_',samp,'.csv'),row.names = 1,stringsAsFactors = F)
counts <- counts[filterByExpr(counts,group = meta$Group),]
CPMcounts = cpm(counts,log = T)
cntTochk = CPMcounts[intersect(tochk,rownames(CPMcounts)),]
#write.csv(cntTochk,'cntTochk_org.csv')

  ## Imputation
  # counts <- counts[filterByExpr(counts,group = meta$Group),]
  # counts[counts==0] <- NA
  # counts <- as.matrix(counts)
  # 
  #   ### Test rowmax, k
  # imputer <- impute.knn(counts)
  # imputer2 <- impute.knn(counts,rowmax = 0.8,k=5)
  # 
  # counts <- imputer$data
  # write.csv(counts,paste0(samp,'.Imputedcounts.csv'))
  # 
  # counts2 <- imputer2$data
  # write.csv(counts2,paste0(samp,'.Imputedcounts2.csv'))
  # CPMcounts = cpm(counts,log = T)
  # cntTochk = CPMcounts[intersect(tochk,rownames(CPMcounts)),]
  # write.csv(cntTochk,'cntTochk_imputed.csv')

meta$Label = factor(meta$Label,
               levels = c("WT(Base)","Apoe(Base)","Ldlr(Base)","WT(Con)","WT(Rux)","Apoe(Con)","Apoe(Rux)","Ldlr(Con)","Ldlr(Rux)")
)

# Plot and combine
for(x in rownames(cntTochk)){
  pldata = data.frame(logCPM=cntTochk[x,],
                       Sample=colnames(cntTochk),
                       Geno = meta$Geno,
                       Treatment = meta$Treatment,
                       Label = meta$Label
                       )
  wt_med = median(c(pldata["W1A","logCPM"],pldata["W2A","logCPM"],pldata["W3A","logCPM"]))
  
  p<-ggplot(pldata,aes(x=Label,y=logCPM))+
    geom_boxplot(aes(color=Treatment))+
    geom_hline(yintercept = wt_med,
               linetype=2,color='grey')+
    xlab("")+
    ylab(x)+
    theme_classic()+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle=45,hjust=1))
  
  # if(x==rownames(cntTochk)[2]){p<-p+ggtitle(title)}
  # if(ncol==2){p<-p+ggtitle(title)}
  if(x==rownames(cntTochk)[1]){q<-p}else{q<-q+p}
}

title = paste(samp,paste0(type,version),thres,comp,sep='_')
q<-q+ plot_layout(ncol = ncol,heights = 0.7)+ plot_annotation(title=title)

pdf(paste('Boxplot',samp,type,version,group,thres,comp,'pdf',sep = '.'),
    height = height,width = width)
print(q)
dev.off()

