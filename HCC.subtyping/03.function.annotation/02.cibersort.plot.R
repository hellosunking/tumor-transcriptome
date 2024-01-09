library(GSVA)
library(GSEABase)
library(clusterProfiler)
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(tibble)
library(limma)
library(ggplotify)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(readxl)
###########

data <- read.csv("infiltration_estimation_for_tcga.csv")
data1 <- data[grep(c("cell_type|_CIBERSORT.ABS$"),colnames(data))]
rownames(data1) <- data1$cell_type
sub <- read.table("HCC.cluster.clinical.include.normal.txt",header = T,sep="\t",row.names = 1)
sub <- subset(sub,sub$Subtype!="Normal")

rownames(sub) <- gsub("-01A$","-01",rownames(sub))
data2 <- merge(sub, data1,by="row.names")
data2 <- data2[,c(8:ncol(data2))]
colnames(data2) <- gsub("_CIBERSORT.ABS$","",colnames(data2))
a2 <- tidyr::gather(data2, key = ssgsea, value = expression, -c(Subtype,cell_type))

a2 <- a2[c(460:1530,1:459,1531:1836,1837:2448,2449:2754,2755:3060,3061:3213,3214:3366),]
order <- factor(a2$ssgsea,levels = c("T.cell.CD4..naive","T.cell.CD4..memory.resting",
                                     "T.cell.CD4..memory.activated","T.cell.CD8.","T.cell.follicular.helper",
                                     "T.cell.gamma.delta","T.cell.regulatory..Tregs.","B.cell.plasma",
                                     "B.cell.naive","B.cell.memory","NK.cell.resting","NK.cell.activated",
                                     "Monocyte","Macrophage.M0","Macrophage.M1","Macrophage.M2",
                                     "Myeloid.dendritic.cell.resting","Myeloid.dendritic.cell.activated",
                                     "Mast.cell.resting","Mast.cell.activated","Eosinophil","Neutrophil"
                                     ))
  #CIBERSORT.ABS

library(ggplot2)
p <- ggplot(data=a2, aes(x=order,y=as.numeric(expression),color= Subtype))+
  geom_boxplot()+
  #geom_beeswarm()+
  scale_color_manual(values = c("#65b5f1","#b34df3","#e94849"))+
  theme(axis.title=element_text(face="bold", size=12,colour = 'black'), #
                       axis.text.x =element_text(size=10,colour = 'black',angle = 60,hjust = 1, vjust = 1), #
                       plot.title = element_text(hjust = 0.5,face="bold"), ###
                       panel.border = element_blank(),
                       panel.grid = element_blank(),
                       panel.background = element_blank(), #
                       axis.line = element_line(colour = "black"),
                       legend.key = element_blank() ,#
                       legend.key.size = unit(0.3,"inches"),
        strip.background =  element_rect(fill = NA, colour = NA), 
        strip.text = element_text(size = 10,face = "bold")
  )+
  labs(x="",y="expression",title = "")+stat_compare_means(label = "p.format") #label = "p.signif"
p
ggsave("ssgsea.immune.CIBERSORT.ABS.sort2.pdf",width = 12,height = 5)


