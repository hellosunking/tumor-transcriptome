# Author: ahfyth
#
# R script for 
#

library( plyr     );
library( survival );
library( beeswarm );

argv = commandArgs(T);
if( length(argv) != 3 )
{
	print( 'usage: R --slave --args <in.train.SAHR> <in.test.SAHR> <out.pdf> < plot.R' );
	q();
}

outfileName=argv[3];
pdf( outfileName ,width=16,height=8);
require(stringr)
#layout(matrix(c()))
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,4))

## plot train
raw=read.table( argv[1], head=T );
name = str_split(argv[1],'/')
name = name[[1]][length(name[[1]])]
#raw=subset(raw, raw$Age<70*365 );	## the age issue, discard the too old guys

dat=raw[ order(raw$SAHR), ];

cancer = gsub(".matrix.SAHR","",argv[1])
plot(density(dat$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer))


## Figure 1: positive SAHR vs negative SAHR values
pos=subset(dat, dat$SAHR >= 0);
PCode=rep('Pos', nrow(pos));
pos.labeled=cbind( pos, PCode );
colnames( pos.labeled ) = c( colnames(raw), 'Code' );

neg=subset(dat, dat$SAHR <  0);
NCode=rep('Neg', nrow(neg));
neg.labeled=cbind( neg, NCode );
colnames( neg.labeled ) = c( colnames(raw), 'Code' );

SAHR.labeled=rbind( neg.labeled, pos.labeled );
#colnames( SAHR.labeled ) = c( 'Sample', 'Vital', 'Days', 'SAHR', 'Code' );
attach( SAHR.labeled );

surv = Surv(Days, Vital);
fit  = survfit(  surv~Code, data=SAHR.labeled );
diff = survdiff( surv~Code, data=SAHR.labeled );
pval = 1 - pchisq(diff$chisq, length(diff$n)-1 );

fit.coxph = coxph( surv~Code, data=SAHR.labeled );
hr=exp(coef(fit.coxph));	## hazard ratio

plot( fit, col=c('blue', 'red'), xlab="Time (days) of Train", ylab="Survival Probability", lwd=2, mark.time=T,
		main=paste("P=", sprintf("%.1e", pval), ", HR=", sprintf("%.2f", hr), sep="") );
legend('bottomleft', c(paste("SAHR- (N=", nrow(neg), ")", sep=""),
					   paste("SAHR+ (N=", nrow(pos), ")", sep="")),
	   col=c('blue', 'red'), lty=c(1,1), bty='n', cex=1.1);

detach( SAHR.labeled );

## Figure 2: Survival vs stage
stage = subset( raw, raw$Stage != "NA" );
stage = subset( raw, raw$Stage != "Stage_X" );
if( nrow(stage) > 0 && nrow( count(stage, 'Stage') ) >= 2 )
## discard those with only one stage or no stage info (for UCEC)
{

	stage.cat = count(stage, 'Stage');
	stage.cat.num = nrow( stage.cat );
	color.scheme=c("blue", "purple","red","red4","black");
	cs = head( color.scheme, stage.cat.num );

	attach( stage );
	surv3 = Surv(Days, Vital);
	fit3  = survfit(  surv3~Stage, data=stage );
	diff3 = survdiff( surv3~Stage, data=stage );
	pval3 = 1 - pchisq(diff3$chisq, length(diff3$n)-1 );

	fit.coxph3 = coxph( surv3~Stage, data=stage );
	hr3=exp(coef(fit.coxph3));	## hazard ratio

	hr3.combined=paste(sprintf("%.2f", hr3), collapse=";");
	plot( fit3, col=cs, xlab="Time (days) of Train", ylab="Survival Probability", lwd=2, mark.time=T,
			main=paste("P=", sprintf("%.1e", pval3), ", HR=", sprintf("%.2f", hr3), sep="") );;
	legend('bottomleft', paste("", stage.cat[,1], " (n=", stage.cat[,2], ")", sep=""),
			col=cs, lty=rep(1, stage.cat.num), bty='n', cex=1.1);
	detach( stage );
}


## Figure 3: SAHR value vs stage
stage = subset( raw, raw$Stage != "Stage_X" );
if( nrow(stage) > 0 && nrow( count(stage, 'Stage') ) >= 2 )

{
    boxplot( stage$SAHR~stage$Stage, outline=F, ylim=range(stage$SAHR),
            ylab="SAHR values", xlab="Pathological stage",
            main=paste("SAHR values vs Stage in ", name, sep=""));
    beeswarm( stage$SAHR~stage$Stage, add=T, pch=19, corral="wrap", col=cs );
}





## plot testing
raw=read.table( argv[2], head=T );
name = str_split(argv[2],'/')
name = name[[1]][length(name[[1]])]
#raw=subset(raw, raw$Age<70*365 );	## the age issue, discard the too old guys

dat=raw[ order(raw$SAHR), ];

cancer = gsub(".matrix.SAHR","",argv[2])
plot(density(dat$SAHR),xlab="Value", ylab="Density", lwd=2, main=paste(cancer))

## Figure 1: positive SAHR vs negative SAHR values
pos=subset(dat, dat$SAHR >= 0);
PCode=rep('Pos', nrow(pos));
pos.labeled=cbind( pos, PCode );
colnames( pos.labeled ) = c( colnames(raw), 'Code' );

neg=subset(dat, dat$SAHR <  0);
NCode=rep('Neg', nrow(neg));
neg.labeled=cbind( neg, NCode );
colnames( neg.labeled ) = c( colnames(raw), 'Code' );

SAHR.labeled=rbind( neg.labeled, pos.labeled );
#colnames( SAHR.labeled ) = c( 'Sample', 'Vital', 'Days', 'SAHR', 'Code' );
attach( SAHR.labeled );

surv = Surv(Days, Vital);
fit  = survfit(  surv~Code, data=SAHR.labeled );
diff = survdiff( surv~Code, data=SAHR.labeled );
pval = 1 - pchisq(diff$chisq, length(diff$n)-1 );

fit.coxph = coxph( surv~Code, data=SAHR.labeled );
hr=exp(coef(fit.coxph));	## hazard ratio

plot( fit, col=c('blue', 'red'), xlab="Time (days) of test", ylab="Survival Probability", lwd=2, mark.time=T,
	main=paste("P=", sprintf("%.1e", pval), ", HR=", sprintf("%.2f", hr), sep="") );
legend('bottomleft', c(paste("SAHR lower than median (N=", nrow(neg), ")", sep=""),
	paste("SAHR hihger than median (N=", nrow(pos), ")", sep="")),
	col=c('blue', 'red'), lty=c(1,1), bty='n', cex=1.1);

detach( SAHR.labeled );


## Figure 2: Survival vs stage
stage = subset( raw, raw$Stage != "NA" );
stage = subset( raw, raw$Stage != "Stage_X" );
if( nrow(stage) > 0 && nrow( count(stage, 'Stage') ) >= 2 )
## discard those with only one stage or no stage info (for UCEC)
{

	stage.cat = count(stage, 'Stage');
	stage.cat.num = nrow( stage.cat );
	color.scheme=c("blue", "purple","red","red4","black");
	cs = head( color.scheme, stage.cat.num );

	attach( stage );
	surv3 = Surv(Days, Vital);
	fit3  = survfit(  surv3~Stage, data=stage );
	diff3 = survdiff( surv3~Stage, data=stage );
	pval3 = 1 - pchisq(diff3$chisq, length(diff3$n)-1 );

	fit.coxph3 = coxph( surv3~Stage, data=stage );
	hr3=exp(coef(fit.coxph3));	## hazard ratio

	hr3.combined=paste(sprintf("%.2f", hr3), collapse=";");
	plot( fit3, col=cs, xlab="Time (days) of test", ylab="Survival Probability", lwd=2, mark.time=T,
		main=paste("P=", sprintf("%.1e", pval3), ", HR=", sprintf("%.2f", hr3), sep="") );;
	legend('bottomleft', paste("", stage.cat[,1], " (n=", stage.cat[,2], ")", sep=""),
			col=cs, lty=rep(1, stage.cat.num), bty='n', cex=1.1);
	detach( stage );
}

## Figure 3: SAHR value vs stage
stage = subset( raw, raw$Stage != "Stage_X" );
if( nrow(stage) > 0 && nrow( count(stage, 'Stage') ) >= 2 )
{
    boxplot( stage$SAHR~stage$Stage, outline=F, ylim=range(stage$SAHR),
            ylab="SAHR values", xlab="Pathological stage",
            main=paste("SAHR values vs Stage in ", name, sep=""));
    beeswarm( stage$SAHR~stage$Stage, add=T, pch=19, corral="wrap", col=cs );
}

par(opar)
