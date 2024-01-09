#!/bin/sh
#  Author: 
#  Date: 2023-12-27
#  Description: 
#

wkdir=/lustre/home/dxhu/project/project4_desurv/paper_figure/all.program/01.SAHR

export dir=$1;
while read Normal  Tumor
do
   /lustre/home/ksun/bin/R --slave --args $dir/all.sample.list $dir/all.sample.list $Normal $Tumor <$wkdir/02.run.DESeq2.no.spike.in.R
done <$dir/tum.nor.number
echo "DEG analysis finished !!"

