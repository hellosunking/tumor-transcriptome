#
# Author: Kun Sun @SZBL (sunkun@szbl.ac.cn)
# Date  :
#
# R script for 
#

argv = commandArgs(T);
if( length(argv) < 3 ) {
	print( 'usage: R --slave --args <in.surv> <surv.result> <cancer.type> [p=1] < do.surv.analysis.R' );
	q();
}

library(survival);
library(beeswarm);

p=1;
if( length(argv) > 3 ) {
	p = as.numeric(argv[4]);

}

#outfileName = paste0(argv[1], ".pdf");
#pdf( outfileName );

dat=read.table( argv[1], head=T);
dat$FPKM  = log2( dat$FPKM + 1 );

cancer=argv[3]
#
dat.diseaes = subset(dat, Type==cancer);
## Kruskal-Wallis rank sum test
#kw=kruskal.test(dat.diseaes$FPKM, dat.diseaes$Stage);
#boxplot(dat$FPKM ~ dat$Stage, outline=F,
#		col=c("grey", rainbow(length(levels(dat$Stage)))),
#		main=paste0("P(DEG)=", sprintf("%.4g",p), ", P(stage)=", sprintf("%.4f", kw$p.value) ),
#		ylab="Gene expression (log2-FPKM)", xlab="Adj, or AJCC stage");
#beeswarm( dat$FPKM ~ dat$Stage, pch=19, add=T, corral="wrap" );


## survival
n.high = nrow( subset(dat.diseaes, Category=="H") );
n.low  = nrow( subset(dat.diseaes, Category=="L") );
attach( dat.diseaes );
surv = Surv(Days, Vital);
fit  = survfit( surv ~ Category );
diff = survdiff( surv ~ Category );
pval = 1 - pchisq(diff$chisq, length(diff$n) - 1);

fit.coxph = coxph( surv~Category, data=dat.diseaes );
hr = exp( coef(fit.coxph) );	## hazard ratio
info = paste( argv[1], pval,diff$obs[1], diff$obs[2], diff$exp[1], diff$exp[2], hr, sep="\t");
write( info, file=argv[2], append=T );

#plot( fit, col=c('red', 'blue'), xlab="Time (days)", ylab="Survival Probability", lwd=4, mark.time=T,
#		main=paste0(cancer," P=", sprintf("%.4g", pval), ", HR=", sprintf("%.4f", hr)) );
#legend('bottomleft', paste0( c("Higher expressed", "Lower expressed"), " (n=", c(n.high, n.low), ")" ),
#		col=c('red', 'blue'), lty=c(1,1), bty='n', cex=1.1);

#dev.off();
