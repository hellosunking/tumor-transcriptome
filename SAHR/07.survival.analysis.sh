#!/bin/bash


export dir=$1
pd=`pwd`;
cd $dir;
wkdir=/lustre/home/dxhu/project/project4_desurv/paper_figure/all.program/01.SAHR
TCGA=/lustre/home/dxhu/data/TCGA/RNA-seq.new.pipeline.2022

#sh survival.analysis.ALL.DEG.hr.sh dir
####################################################################
perl $wkdir/06.smy.DEG.pl $dir/all.sample.list.DESeq2.tsv >$dir/All.DEGs.list.symbol

timestamp=`date +%y%m%d%H%M%S`
[ -s Gene.exp.cutoff ] && mv Gene.exp.cutoff Gene.exp.cutoff.$timestamp
mkdir -p surv.all
cd surv.all

[ -s surv.summary ] && mv surv.summary surv.summary.$timestamp


cnt=0
#Gene   Symbol  Biotype DEG     fc      p
while read Gene Symbol Biotype DEG fc p
do
        if [[ $Gene =~ ^ENSG ]]
        then
                let cnt=$cnt+1
                echo -en "\r$cnt: $Gene"
                perl $wkdir/07.sub.prepare.surv.pl $Gene $subtype ../ > $Gene.$Symbol.$DEG
		/lustre/home/ksun/bin/R --slave --args $Gene.$Symbol.$DEG surv.summary $subtype $p < $wkdir/08.do.surv.analysis.R 2>>err.log
                mv $Gene.$Symbol.$DEG $Gene

        fi
done < ../All.DEGs.list.symbol

echo -e "\rDone: $cnt genes analyzed.";
cp surv.summary ../ 
cd ../
echo "survival analysis finished for DEGS !!!!"
########################################################################

sed -i 's/.Up/\tUp/' surv.summary
sed -i 's/.Down/\tDown/' surv.summary

sed -i '1i Gene\tDEG\tpval\tobs1\tobs2\texp1\texp2\tHR' surv.summary

echo "pro1 finished !!"
