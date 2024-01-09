#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) != 4 ) {
	print( 'usage: R --slave --args <in.matrix> <out.prefix> <N1> <N2> < plot.R' );
	q();
}

suppressPackageStartupMessages({ library( DESeq2 ); });

RNAseq = read.table( argv[1], head=T, row.names=1 );
all.gene = rownames( RNAseq );

## regress out Ribosomal RNA and Mitochonrial genes
RP.MT = read.table( "ENSEMBL.v101.RP.MT" );
keep.gene = all.gene[!(all.gene %in% c(RP.MT$V1, RP.MT$V2))];	## both ENS and Symbol
RNAseq = RNAseq[ keep.gene, ];

group_list = c( rep("0.Ctrl", as.numeric(argv[3]) ), rep("1.Treatment", as.numeric(argv[4])) );
colData = data.frame( row.names=colnames(RNAseq), group_list=group_list );

dds = DESeqDataSetFromMatrix( countData=RNAseq, colData=colData, design = ~ group_list );
# FILTER NONE-EXPRESSED GENES TO REDUCE FALSE-NEGATIVE
# It seems that DESEQ2 will do this filter itself
#dds = dds[rowSums(counts(dds)) > 1, ];
dds = estimateSizeFactors( dds );
dsq = DESeq( dds );
res = results( dsq );
write.table( res, paste0(argv[2], ".DESeq2.tsv"), sep="\t", quote=F, col.names=NA);

normalized.data = counts( dds, normalized=T );
write.table( normalized.data, file=paste0(argv[2], ".DESeq2.norm.exp.tsv.gz"), sep="\t", quote=F);

