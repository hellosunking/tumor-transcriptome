#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

suppressPackageStartupMessages({
	library(Seurat);library(plyr);library(dplyr);library(cowplot);
	library(ggplot2);library(gridExtra);library(ggthemes);library(patchwork);
	library(gridExtra);library(forcats);library(stringr);
});
options( stringsAsFactors=F );

ginfo = read.table( "ENSEMBL.v101.info" );
RP = ginfo[ grep("^RP[SL]", ginfo$V2), ];
MT = ginfo[ grep("^MT-", ginfo$V2), ];

metadata=read.table("HCC.clinical.race.info",header=T,row.names=1);
rownames(metadata)=gsub("-",".",row.names(metadata));
a = read.table( "asian.sample.list.unstrand", head=T, row.names=1 );
rownames(a) = as.character( gsub("\\.\\d+", "", rownames(a)) );
sc.raw = CreateSeuratObject(counts=a, project='HCC',meta.data = metadata);
all.gene = rownames(a);
keep_genes = all.gene[!(all.gene %in% c(RP$V1, MT$V1))];  

sce = subset(sc.raw, features=keep_genes);
sce = NormalizeData( sce, normalization.method='LogNormalize', scale.factor=10000 );
sce = FindVariableFeatures(sce, selection.method='vst', nfeatures=2000);
#sce = FindVariableFeatures(sce, selection.method='mean.var.plot');
sce = ScaleData(sce);
sce = RunPCA(sce, features=VariableFeatures(sce));
sce = FindNeighbors(sce, dims=1:30);
sce = FindClusters(sce, resolution = 0.7);

sce = RunTSNE(sce, dims=1:30);
sce = RunUMAP(sce, dims=1:30);

sce=SetIdent(sce,value='RNA_snn_res.0.7');
sce$seurat_clusters=sce$RNA_snn_res.0.7;
write.table(Idents(sce), "Cluster.txt", sep="\t",col.names=F, quote=F);

cluster_markers = FindAllMarkers(sce, logfc.threshold=0.2, min.pct=0.25);
cluster_markers = cluster_markers[with(cluster_markers,order(cluster, p_val_adj, -avg_log2FC)),];
cluster_markers_list = split(cluster_markers, cluster_markers$cluster)
all.markers = cluster_markers %>% group_by(cluster);
write.table(all.markers, file="All.marker.Genes.txt", sep="\t", quote=F, row.names=F);

save( sce, file="HCC.seurat.rds" );

## Plots
outfileName = 'TCGA.HCC.asian.meta.pdf';
pdf( outfileName, width=6, height=4 );
UMAPPlot(sce, pt.size=2, label=T);
TSNEPlot(sce, pt.size=2, label=T);
UMAPPlot(sce, pt.size=2, label=F, group.by = "race")
TSNEPlot(sce, pt.size=2, label=F, group.by = "race")
UMAPPlot(sce, pt.size=2, label=F, group.by = "ajcc_pathologic_stage")
TSNEPlot(sce, pt.size=2, label=F, group.by = "ajcc_pathologic_stage")
UMAPPlot(sce, pt.size=2, label=F, group.by = "vital_status")
TSNEPlot(sce, pt.size=2, label=F, group.by = "vital_status")
#Gender
UMAPPlot(sce, pt.size=2, label=F, group.by = "gender")
TSNEPlot(sce, pt.size=2, label=F, group.by = "gender")
