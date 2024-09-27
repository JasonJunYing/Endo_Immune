library(limma)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)
library(rstatix)

# Preprocess

data <- read.csv('New Plex_result_analytes.csv',
                 row.names = NULL,header = F)
data <- data[-1,]
k=1
for(i in 1:nrow(data)){
  if(data$V2[i]=="        "){
    sub<-data[(i+12):(i+66),]
    colnames(sub)<-data[2,]
    write.csv(sub,paste0('A',k,'.csv'))
    k=k+1
  }
}

analytes <- c("IL-6","IL-10","MCP-1","IFN-γ","TNF","IL12p70")

# Plotting Boxplot

for(i in 1:6){
  # Separately annotate (could be optimized)
  result <- read.csv(paste0('A',i,'.csv'))
  result$FinalCC <- result$Final.CC.........*2
  resultW3 <- result[!startsWith(result$Results.File,'W0'),]
  resultW3$Time <- "W3"
  resultW3$SampleID <- strsplit2(resultW3$Results.File,'-')[,1]
  resultW3$Group<-substr(resultW3$Results.File,1,2)
  resultW3$Geno <- substr(resultW3$Group,1,1)
  resultW3$Treatment<-substr(resultW3$Group,2,2)
  resultW3$Treatment <- recode(resultW3$Treatment,
                               C='Con',
                               R='Rux'
  )
  resultW3$Geno <- recode(resultW3$Geno,
                          L='Ldlr',
                          W='WT',
                          A='Apoe'
  )
  resultW3$Label <- paste0(resultW3$Geno,'(',resultW3$Treatment,')')
  resultW3$Label = factor(resultW3$Label,
                          levels = c("WT(Con)","WT(Rux)","Apoe(Con)","Apoe(Rux)","Ldlr(Con)","Ldlr(Rux)")
  )
  
  
  resultW0 <- result[startsWith(result$Results.File,'W0'),]
  resultW0$Time <- "W0"
  resultW0$SampleID <- substr(resultW0$Results.File,4,5)
  resultW0$Group<-paste0(substr(resultW0$Results.File,4,4),"C")
  resultW0 <- resultW0[resultW0$SampleID!='L2',]
  resultW0$Treatment <- "Con"
  
  resultW0$Geno <- substr(resultW0$Group,1,1)
  resultW0$Geno <- recode(resultW0$Geno,
                          L='Ldlr',
                          W='WT',
                          A='Apoe'
  )
  resultW0$Label <- paste0(resultW0$Geno,'(',resultW0$Treatment,')')
  resultW0$Label = factor(resultW0$Label,
                          levels = c("WT(Con)","WT(Rux)","Apoe(Con)","Apoe(Rux)","Ldlr(Con)","Ldlr(Rux)")
  )
  
  result2 <- rbind(resultW0,resultW3)
  result2$Time <- recode(result2$Time,W0="Base(Week 0)",W3="Week 3")
  
  # Independent T-test
  stats <- result2 %>%
    group_by(Time) %>%
    pairwise_t_test( # pairwise, not paired
      FinalCC ~ Label, paired = F,
      p.adjust.method = "bonferroni"
    )
  stats <- stats %>% add_xy_position(x = "Label",
                                             step.increase = 0.005)
  stats$y.position <- stats$y.position-100
  stats$G1 = strsplit2(stats$group1,'\\(')[,1]
  stats$T1 = strsplit2(stats$group1,'\\(')[,2]
  stats$G2 = strsplit2(stats$group2,'\\(')[,1]
  stats$T2 = strsplit2(stats$group2,'\\(')[,2]
  stats <- stats[(stats$G1==stats$G2)|(stats$T1==stats$T2),]
  
  stats$cytokine <- analytes[i]
  if(i==1)statAll<-stats
  else statAll <- rbind(statAll,stats)
  
  #Plotting
  m <- ggplot(result2,aes(x=Label,y=FinalCC))+
    geom_boxplot()+
    xlab(analytes[i])+
    ylab("FinalCC (pg/ml)")+
    ylab("")+
    theme_bw()+
    theme(legend.position = "none")+
    rotate_x_text(45)+
    facet_wrap(~Time,)+
    stat_pvalue_manual(
      stats, label = "p.signif",
      hide.ns = 'p',
      step.increase = 0.05,tip.length = 0
    )
  if(i==1)n<-m
  else n <- n+m
}

n<-n+ plot_layout(ncol = 2,widths = 0.1,heights = 0.1)
ggsave('Cytokine_All.pdf',n,width = 7,height = 7)

