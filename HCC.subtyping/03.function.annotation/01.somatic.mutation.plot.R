library(maftools)
library(ggpubr)
library(rstatix)
mytheme <- theme(axis.title=element_text(face="bold", size=18,colour = 'black'), #
                 axis.text=element_text(size=18,colour = 'black'), #
                 plot.title = element_text(hjust = 0.5,face="bold"), ###
                 panel.border = element_blank(),
                 panel.grid = element_blank(),
                 panel.background = element_blank(), #
                 axis.line = element_line(colour = "black"),
                 legend.key = element_blank() ,#
                 legend.key.size = unit(0.3,"inches")
)
laml.maf <- read.maf("TCGA.LIHC.mutect.a630f0a0-39b3-4aab-8181-89c1dde8d3e2.DR-10.0.somatic.maf", clinicalData ="HCC.cluster.clinical.txt", isTCGA = T)
data <- read.table("HCC.cluster.clinical.txt",header=T,row.names = 1)
sid <- rownames(data)
laml.maf <- subsetMaf(maf=laml.maf,tsb=sid)

#write.mafSummary(maf = laml.maf, basename = 'laml.maf')
pdf("mutect.mutation.pdf",width = 8,height = 6)
plotmafSummary(maf = laml.maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)


col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')
fabcolors = c("#65b5f1","#b34df3","#e94849")
names(fabcolors) = c("C0","C1","C2")
fabcolors = list(Subtype = fabcolors)
colsAnn <- list('Subtype' = c("C0"="#65b5f1","C1"="#b34df3","C2"="#e94849"),
                'Gender' = c('female' = '#ff7f00', 'male' = '#1f78b4'),
                'Stage' = c('Stage_I' = '#00AFBB', 'Stage_II' = '#E7B800','Stage_III' = 'purple','Stage_IV' = '#FC4E07')
                )

oncoplot(maf = laml.maf, colors=col, top=20,clinicalFeatures =c('Subtype','Gender','Stage'),fontSize = 0.8,
         annotationOrder=c("C0","C1","C2"),keepGeneOrder=T,anno_height=1.5,bgCol = "white",
         sortByAnnotation = T,annotationColor = colsAnn,removeNonMutated=F)

dev.off()
#############
pdf("mutation.stat.boxplot.pdf",width = 8,height = 6)
info <- getSampleSummary(laml.maf)
info <- data.frame(info,row.names = info$Tumor_Sample_Barcode)
sample <- merge(data,info,by="row.names")
my_comparisons <- list(c("C1","C0"),c("C1","C2"),c("C0","C2"))
sample$Subtype <- factor(sample$Subtype,levels = c("C0","C1","C2"))
p1 <- ggplot(data=sample,aes(x=Subtype,y=total,color=Subtype))+
  geom_boxplot()+
  geom_beeswarm()+
  theme_classic()+
  labs(x="Subtype",y="Total mutation number")+
  #mytheme+
  #guides(color=T)+
 # scale_color_manual(values = c("#194FE7","purple","#e31a1c"))+
  scale_color_manual(values = c("#65b5f1","#b34df3","#e94849"))+
  stat_compare_means(method = "kruskal.test")+
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",aes(label=..p.signif..),hide.ns = T)+
  theme(axis.text = element_text(size=30,color = "black"),
       axis.title = element_text(size=30),
       legend.text = element_text(size=30))
p1
dev.off()
kruskal.test(sample$total~sample$Subtype)
dunn_test(total~Subtype,data=sample,p.adjust.method = "bonferroni")

################ differences of genes among subtypes ##################
C0 <- subset(data,data$Subtype=="C0")
C1 <- subset(data,data$Subtype=="C1")
C2 <- subset(data,data$Subtype=="C2")
C0 <- rownames(C0)
C1 <- rownames(C1)
C2 <- rownames(C2)
laml.c0 <- subsetMaf(maf=laml.maf,tsb=C0)
laml.c1 <- subsetMaf(maf=laml.maf,tsb=C1)
laml.c2 <- subsetMaf(maf=laml.maf,tsb=C2)
info.c0 <- getGeneSummary(laml.c0)
info.c1 <- getGeneSummary(laml.c1)
info.c2 <- getGeneSummary(laml.c2)

gene <- c("TP53","TTN","CTNNB1","MUC16","MUC4","PCLO","OBSCN","AXIN1","RYR1","ALB","CSMD3","FLG","RB1","ABCA13","APOB","HMCN1","RYR2","CACNA1E","ARID1A","CCDC168")
info.c0 <- info.c0[info.c0$Hugo_Symbol %in% gene,c("Hugo_Symbol","total")]
info.c1 <- info.c1[info.c1$Hugo_Symbol %in% gene,c("Hugo_Symbol","total")]
info.c2 <- info.c2[info.c2$Hugo_Symbol %in% gene,c("Hugo_Symbol","total")] ##info.c2 ## 缺少 AXIN1和RB1
add <- data.frame(Hugo_Symbol=c("AXIN1","RB1"),total=c(0,0))
info.c2 <- rbind(info.c2,add)

info.c0$Subtype <- "C0"
info.c1$Subtype <- "C1"
info.c2$Subtype <- "C2" 
info.c0$Sample <- 76
info.c1$Sample <- 56
info.c2$Sample <- 21

dat <- rbind(info.c0,info.c1,info.c2)
dat$Hugo_Symbol <- factor(dat$Hugo_Symbol,levels = gene)
dat$prop <- dat$total/dat$Sample*100
#write.table(dat,file = "Gene.mutation.prop.txt",quote=F,sep="\t",row.names = F)
ggplot(dat,aes(x=Hugo_Symbol,y=prop,fill=Subtype))+
  geom_bar(position=position_dodge(preserve = "single"), stat="identity",width = 0.7,alpha=1)+
  #geom_col(position = position_dodge(preserve = "single"),width = 0.7)+
  mytheme+
  scale_y_continuous(expand = c(0,0,0,5))+
  #ylim(0,85)+
  labs(x="",y="Mutation frequency (80%)")+
  theme(legend.position = c(0.9,0.7),
        #legend.direction = "horizontal",
        legend.key.size = unit(0.5,"inches"),
        legend.text = element_text(size=20),
        legend.title = element_text(size=24),
        axis.text=element_text(size=20,colour = 'black'),
        axis.title = element_text(size=24,colour = 'black'),
        #panel.background = theme_rect(fill = "transparent",colour = NA)
        plot.background = element_rect(fill = "transparent",colour = NA)
        )+
  scale_fill_manual(values = c("#65b5f1","#b34df3","#e94849"))


ggsave("top20.gene.mutation.stat.pdf",width = 24,height = 5)
#############
dat$no_mut <- dat$Sample-dat$total
p <- data.frame(list("Gene","compare","pvalue"))
for (i in gene){
print(i)
  n <- which(dat$Hugo_Symbol==i)
t <- dat[dat$Hugo_Symbol %in% i,c(3,2,6)]
c0.c1 <- t[1:2,c(2,3)]
c0.c2 <- t[c(1,3),c(2,3)]
c1.c2 <- t[c(2,3),c(2,3)]
c0.c1<- chisq.test(c0.c1,correct = F)
c0.c2 <- chisq.test(c0.c2,correct = F)
c1.c2 <- chisq.test(c1.c2,correct = F)

p <-  rbind(p,list(i,"c0.vs.c1",c0.c1$p.value))
p <-  rbind(p,list(i,"c0.vs.c2",c0.c2$p.value))
p <-  rbind(p,list(i,"c1.vs.c2",c1.c2$p.value))
}
write.csv(p,file = "top20.gene.mutation.frequency.csv")
####################################
