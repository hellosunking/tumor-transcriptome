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
###########
gene.set = read.table("cellreports.txt", stringsAsFactors = F, sep = "\t", header = F, check.names = F)
gene.set <- gene.set %>%
  column_to_rownames("V1")%>%t()
a <- gene.set
a <- a[1:nrow(a),]
set <- colnames(a)
l <- list()
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}

exp <- read.table("exp.file",header = T,sep="\t",row.names = 1)
colnames(exp) = gsub("[.]", "-", colnames(exp))
exp <- as.matrix(exp)

exp <- exp[,-c(1:6)]
ssgsea<- gsva(log2(exp+1),l, method='ssgsea', kcdf='Gaussian',abs.ranking=F) 

################# plot the results ###############
sub <- read.table("HCC.cluster.clinical.include.normal.txt",header = T,sep="\t",row.names = 1)
sub <- subset(sub,sub$Subtype!="Normal")

############## boxplot 统计##########
a <- ssgsea %>% t() %>% as.data.frame()
sub <- read.table("HCC.cluster.clinical.include.normal.txt",header = T,sep="\t",row.names = 1)
sub <- subset(sub,sub$Subtype!="Normal")
a1 <- merge(a,sub,by="row.names")
a1 <- a1[,-c(30:35)]
a2 <- tidyr::gather(a1, key = ssgsea, value = expression, -c(Subtype,Row.names))
order <- factor(a2$ssgsea,levels = unique(a2$ssgsea))

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
  labs(x="",y="ssGSEA score")+stat_compare_means(label = "p.signif")
p
ggsave("ssgsea.immune.pdf",width = 12,height = 4.5)

################### kegg ssgsea################
gene.set <- clusterProfiler::read.gmt("c2.cp.kegg.v2023.1.Hs.symbols.gmt")
 gene.set = split(gene.set$gene, gene.set$term)
 str(head(gene.set))
 ssgsea.2<- gsva(log2(exp+1),gene.set, method='gsva', kcdf='Gaussian',abs.ranking=F) 
 kegg.all <- ssgsea.2[,row.names(sub)]
 ha = HeatmapAnnotation(Subtype=sub$Subtype, Gender=sub$Gender,Stage=sub$Stage,
                        col = list(Subtype = c("Normal" = "#CCCCCC", "C0" = "#65b5f1","C1"="#b34df3","C2"="#e94849"),
                                   Gender = c('female' = '#ff7f00', 'male' = '#1f78b4'),
                                   Stage = c('Stage_I' = '#00AFBB', 'Stage_II' = '#E7B800','Stage_III' = 'purple','Stage_IV' = '#FC4E07')
                        ))
 p2 <- Heatmap(kegg.all, name = "ssGSEA score",show_column_names = F,show_row_dend = F,
              row_names_gp = grid::gpar(fontsize = 5), top_annotation = ha,cluster_columns = F,cluster_rows = T)
p2

#### limma difference analysis  ####
##
c0 <- subset(sub,sub$Subtype=="C0")
c0$Subtype2 <- "C0"
c1 <- subset(sub,sub$Subtype=="C1")
c1$Subtype2 <- "C1"
c2 <- subset(sub,sub$Subtype=="C2")
c2$Subtype2 <- "C2"
all <- rbind(c0,c1,c2)

c1_c2 <- subset(sub,sub$Subtype=="C1" | sub$Subtype=="C2")
c1_c2$Subtype2 <- "C1_C2"
c0_c2 <- subset(sub,sub$Subtype=="C0" | sub$Subtype=="C2")
c0_c2$Subtype2 <- "C0_C2"
c0_c1 <- subset(sub,sub$Subtype=="C0" | sub$Subtype=="C1")
c0_c1$Subtype2 <- "C0_C1"

################################################# C0 vs other ####################
subtype <- rbind(c0,c1_c2)
group_list <- factor(c(c0$Subtype,c1_c2$Subtype2))
design <- model.matrix(~ group_list+0)
rownames(design) <- c(row.names(c0),row.names(c1_c2))
colnames(ssgsea.2) = as.character( gsub("\\.", "-", colnames(ssgsea.2)) );
ssgsea.2 <- ssgsea.2[,row.names(design)]
comparE <- makeContrasts(group_listC0 - group_listC1_C2, levels=design)
fiT <- lmFit(ssgsea.2, design)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
keggDiff <- topTable(fiT3, coef=1, number=200)
keggDiff_c0 <- keggDiff
padj_cutoff=0.05
log2FC_cutoff=0.5
keep_c0 <- rownames(keggDiff[keggDiff$adj.P.Val < padj_cutoff & abs(keggDiff$logFC)>log2FC_cutoff, ])
length(keep_c0)
###################################### C1 vs other #################################
subtype <- rbind(c1,c0_c2)
group_list <- factor(c(c1$Subtype,c0_c2$Subtype2))
design <- model.matrix(~ group_list+0)
rownames(design) <- c(row.names(c1),row.names(c0_c2))
colnames(ssgsea.2) = as.character( gsub("\\.", "-", colnames(ssgsea.2)) );
ssgsea.2 <- ssgsea.2[,row.names(design)]
comparE <- makeContrasts(group_listC1 - group_listC0_C2, levels=design)
fiT <- lmFit(ssgsea.2, design)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
keggDiff <- topTable(fiT3, coef=1, number=200)
keggDiff_c1 <- keggDiff
padj_cutoff=0.05
log2FC_cutoff=0.5
keep_c1 <- rownames(keggDiff[keggDiff$adj.P.Val < padj_cutoff & abs(keggDiff$logFC)>log2FC_cutoff, ])
length(keep_c1)
#######################################C2 vs other #################################
subtype <- rbind(c2,c0_c1)
group_list <- factor(c(c2$Subtype,c0_c1$Subtype2))
design <- model.matrix(~ group_list+0)
rownames(design) <- c(row.names(c2),row.names(c0_c1))
colnames(ssgsea.2) = as.character( gsub("\\.", "-", colnames(ssgsea.2)) );
ssgsea.2 <- ssgsea.2[,row.names(design)]
comparE <- makeContrasts(group_listC2 - group_listC0_C1, levels=design)
fiT <- lmFit(ssgsea.2, design)
fiT2 <- contrasts.fit(fiT, comparE)
fiT3 <- eBayes(fiT2)
keggDiff <- topTable(fiT3, coef=1, number=200)
keggDiff_c2 <- keggDiff
padj_cutoff=0.05
log2FC_cutoff=0.5
keep_c2 <- rownames(keggDiff[keggDiff$adj.P.Val < padj_cutoff & abs(keggDiff$logFC)>log2FC_cutoff, ])
length(keep_c2)
######################################################################################

keep <- c(keep_c0,keep_c1,keep_c2)
keep <- unique(keep)
kegg <- kegg.all[keep,]

ha = HeatmapAnnotation(Subtype=all$Subtype, Gender=all$Gender,Stage=all$Stage,
	col = list(Subtype = c("Normal" = "#CCCCCC", "C0" = "#65b5f1","C1"="#b34df3","C2"="#e94849"),
	Gender = c('female' = '#ff7f00', 'male' = '#1f78b4'),
	Stage = c('Stage_I' = '#00AFBB', 'Stage_II' = '#E7B800','Stage_III' = 'purple','Stage_IV' = '#FC4E07') ))
mycols <- colorRamp2(breaks = c(-0.9, 0, 0.9),colors = c("#349FD4", "white", "red"))
rownames(kegg) <- gsub("KEGG_","",row.names(kegg))
p3 <- Heatmap(kegg, name = "GSVA score",show_column_names = F,show_row_dend = F,col=mycols,
              row_names_gp = grid::gpar(fontsize = 8), top_annotation = ha,cluster_columns = F,cluster_rows = T)
p3 <- as.ggplot(p3)
p3
ggsave("KEGG.GSVA.pdf",width = 9,height =6 )

