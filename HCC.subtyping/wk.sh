#!/bin/bash
perl data.pl HCC.clinical.asian asian.tumor.sample >asian.sample.list.unstrand
R4 --slave <do.subtyping.R
R --slave <replot.old.R
perl attach.subtype.pl HCC.clinical.race.info Cluster.txt >HCC.cluster.clinical.txt
R --slave <do.surv.R

##C1.vs.C0_C2## DEG analyses##
perl data.C1.pl HCC.cluster.clinical.txt C1.vs.C0_C2.txt >C1.vs.C0_C2.unstrand.list
while read C0_C2  C1
do
        R --slave --args C1.vs.C0_C2.unstrand.list C1.vs.C0_C2 $C0_C2 $C1 <2.run.DESeq2.no.spike.in.R
done <C1.vs.C0_C2.txt

perl 3.filter.pl C1.vs.C0_C2.DESeq2.tsv >C1.vs.C0_C2.DESeq2.tsv.filter
perl 4.attach.symbol.pl C1.vs.C0_C2.DESeq2.tsv.filter >C1.vs.C0_C2.DESeq2.tsv.filter.symbol
######## C0 VS. C1_C2 #############
perl data.C0.pl HCC.cluster.clinical.txt C0.vs.C1_C2.txt >C0.vs.C1_C2.unstrand.list
while read C1_C2 C0
do
   R --slave --args C0.vs.C1_C2.unstrand.list C0.vs.C1_C2 $C1_C2 $C0 <2.run.DESeq2.no.spike.in.R
done <C0.vs.C1_C2.txt
perl 3.filter.pl C0.vs.C1_C2.DESeq2.tsv >C0.vs.C1_C2.DESeq2.tsv.filter
perl 4.attach.symbol.pl C0.vs.C1_C2.DESeq2.tsv.filter >C0.vs.C1_C2.DESeq2.tsv.filter.symbol

############## C2 VS. C0_C1 #############
perl data.C2.pl HCC.cluster.clinical.txt C2.vs.C0_C1.txt >C2.vs.C0_C1.unstrand.list
while read C0_C1 C2
do
   R --slave --args C2.vs.C0_C1.unstrand.list C2.vs.C0_C1 $C0_C1 $C2 <2.run.DESeq2.no.spike.in.R
done <C2.vs.C0_C1.txt
perl 3.filter.pl C2.vs.C0_C1.DESeq2.tsv >C2.vs.C0_C1.DESeq2.tsv.filter
perl 4.attach.symbol.pl C2.vs.C0_C1.DESeq2.tsv.filter >C2.vs.C0_C1.DESeq2.tsv.filter.symbol









