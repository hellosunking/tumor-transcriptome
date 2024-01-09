#!/bin/bash


export cancer=$1
export subtype=$2
#export prop=$3
pd=`pwd`;

wkdir=/lustre/home/dxhu/project/project4_desurv/paper_figure/all.program/01.SAHR
dir=/lustre/home/dxhu/project/project4_desurv/paper_figure/all.program/results/SAHR
TCGA=/lustre/home/dxhu/data/TCGA/RNA-seq.new.pipeline.2022

#mkdir $subtype
cd $dir/$subtype
#####DEG analysis #######

perl $wkdir/01.data.pre.to.deg.pl $wkdir/data/$subtype.clinical.tsv.clean $cancer $TCGA tum.nor.number >all.sample.list
#
sh $wkdir/02.deg.analysis.sh $dir/$subtype
echo "DEG analysis finished !!"
########
####sample select######
perl $wkdir/03.rand.clinical.info.pl  $wkdir/data/$subtype.clinical.tsv.clean clinical.info $prop
##########
perl $wkdir/04.data.tumor.pl clinical.info.training $cancer $TCGA >tumor.matrix
perl $wkdir/05.data.normal.pl clinical.info.training $cancer $TCGA >normal.matrix

sh $wkdir/07.survival.analysis.sh $dir/$subtype/
 
echo "survival anlysis finished!!!"
###########
perl $wkdir/09.anno.survival.pl surv.summary $subtype >surv.summary.surv.sig.anno


######select the top genes ##########
for Type in Oncogene TumorSuppressor
do
    awk '{print $1"\t"$4"\t"$5"\t"$6}' surv.summary.surv.sig.anno | grep -w $Type | perl -lane 'next unless $F[1]>0 && $F[3]<20 && 1/$F[3]<20; $F[3]=1/$F[3] if $F[3]<1; print join("\t", @F)' | sort -k2,2g | head -n 20
done > surv.summary.top60
for Type in UpSaver DownSaver
do
    awk '{print $1"\t"$4"\t"$5"\t"$6}' surv.summary.surv.sig.anno | grep -w $Type | perl -lane 'next unless $F[1]>0&& $F[3]<20 && 1/$F[3]<20; $F[3]=1/$F[3] if $F[3]<1; print join("\t", @F)' | sort -k2,2g | head -n 10
done >>surv.summary.top60

############# calculate SAHR value ##############
echo "Start to calculate SAHR value:"
perl $wkdir/10.prepare.datatoSAHR.pl surv.summary.top60 surv.all top60.train.matrix
perl $wkdir/11.pre.test.SAHR.pl surv.summary.top60 Gene.exp.cutoff clinical.info.testing top60.test.matrix $cancer $TCGA
perl $wkdir/12.calc.SAHR.pl surv.summary.top60 top60.train.matrix clinical.info.training top60.train.matrix.SAHR
perl $wkdir/12.calc.SAHR.pl surv.summary.top60 top60.test.matrix clinical.info.testing top60.test.matrix.SAHR
R --slave --args top60.train.matrix.SAHR top60.test.matrix.SAHR ${subtype}<$wkdir/14.plot.survival.R
#R --slave <13.plot.paper.figure.R
