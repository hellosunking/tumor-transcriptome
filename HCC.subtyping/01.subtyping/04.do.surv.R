# R script for 
#

library(survival);
library(beeswarm);
library(ggplot2);
library( plyr );

outfileName = 'TCGA.HCC.subtype.clinical.paper.figure.pdf';
pdf( outfileName,width=6,height=5 );

dat = read.table( "HCC.cluster.clinical.txt" ,head =T );
table( dat$Subtype, dat$Stage );
#
## 3 subtypes
attach( dat );
surv = Surv(Days, Vital);
fit  = survfit( surv ~ Subtype );
diff = survdiff( surv ~ Subtype );
pval = 1 - pchisq(diff$chisq, length(diff$n) - 1);
fit.coxph = coxph( surv~Subtype, data=dat );
hr = exp( coef(fit.coxph) );
subtype.cat = count(dat, 'Subtype');
subtype.cat.num = nrow( subtype.cat );
color.scheme=c("blue", "purple","red");
cs = head( color.scheme, subtype.cat.num );

plot( fit, col=cs, ylim=c(0.4, 1),xlab="Time (days)", ylab="Survival Probability", lwd=3, mark.time=T,
		main=paste0("P=", sprintf("%.4g", pval), ", HR=", sprintf("%.4f", hr)) );
legend('bottomleft',paste("", subtype.cat[,1], " (n=", subtype.cat[,2], ")", sep=""),
                        col=cs, lty=rep(1, subtype.cat.num), bty='n', cex=1.1);

## ANNOVA shows age is correlated with subtype
boxplot( dat$Age ~ dat$Subtype, outline=F, ylim=range(dat$Age),
			names=c("Subtype 0","Subtype 1", "Subtype 2"),
			ylab="Age at diagnosis (days)",
			main="Age vs Subtype" );
beeswarm( dat$Age ~ dat$Subtype, pch=19, add=T, col=c("blue","purple", "red") );

aov = aov( Age ~ Subtype, dat );
summary(aov)


tukey = TukeyHSD( aov );
tukey = as.data.frame( tukey$Subtype );
tukey$pair = rownames( tukey );

ggplot(tukey, aes(colour=cut(`p adj`, c(0, 0.01, 0.05, 1),
		label=c("p<0.01", "p<0.05", "Non-Sig")))) +
theme_bw(base_size = 16) +
geom_hline(yintercept=0, lty="11", colour="grey30", size = 1) +
geom_errorbar(aes(pair, ymin=lwr, ymax=upr), width=0.2, size = 1) +
geom_point(aes(pair, diff), size = 2) +
labs(colour="")+
theme(axis.text.x = element_text(size = 14));

# pool B+C as 1 subtype as these two are similar
dat.A = subset( dat, Subtype == 'C0' );
dat.x = subset( dat, Subtype != 'C0' );
dat.x$Subtype = 'C1_C2';
dat2 = rbind( dat.A, dat.x );
#table( dat$Subtype, dat$Stage );
detach( dat );
attach( dat2 );
surv2 = Surv( Days, Vital );
fit2  = survfit(  surv2 ~ Subtype );
diff2 = survdiff( surv2 ~ Subtype );
pval2 = 1 - pchisq(diff2$chisq, length(diff2$n) - 1);
fit.coxph2 = coxph( surv2~Subtype, data=dat2 );
hr2 = exp( coef(fit.coxph2) )['SubtypeC1_C2'];

subtype.cat = count(dat2, 'Subtype');
subtype.cat.num = nrow( subtype.cat );

plot( fit2, col=c("blue", "red"), ylim=c(0.4, 1),
		xlab="Time (days)", ylab="Survival Probability", lwd=4, mark.time=T,
		main=paste0("P=", sprintf("%.4g", pval2), ", HR=", sprintf("%.4f", hr2)) );
legend('bottomleft',paste(subtype.cat[,1], " (n=", subtype.cat[,2], ")", sep=""),
                        col=c("blue","red"), lty=rep(1,2), bty='n', cex=1.1);
#####################C2. vs C0_1
dat.A = subset( dat, Subtype == 'C2' );
dat.x = subset( dat, Subtype != 'C2' );
dat.x$Subtype = 'C0_C1';
dat3 = rbind( dat.A, dat.x );
#table( dat$Subtype, dat$Stage );
detach( dat2 );
attach( dat3 );
surv3 = Surv( Days, Vital );
fit3  = survfit(  surv3 ~ Subtype );
diff3 = survdiff( surv3 ~ Subtype );
pval3 = 1 - pchisq(diff3$chisq, length(diff3$n) - 1);
fit.coxph3 = coxph( surv3~Subtype, data=dat3 );
hr3 = exp( coef(fit.coxph3) )['SubtypeC2'];

subtype.cat = count(dat3, 'Subtype');
subtype.cat.num = nrow( subtype.cat );

plot( fit3, col=c("blue", "red"), ylim=c(0.4, 1),
                xlab="Time (days)", ylab="Survival Probability", lwd=4, mark.time=T,
                main=paste0("P=", sprintf("%.4g", pval3), ", HR=", sprintf("%.4f", hr3)) );
legend('bottomleft',paste(subtype.cat[,1], " (n=", subtype.cat[,2], ")", sep=""),
                        col=c("blue","red"), lty=rep(1,2), bty='n', cex=1.1);

#########C1 vs. C0_2
dat.A = subset( dat, Subtype == 'C1' );
dat.x = subset( dat, Subtype != 'C1' );
dat.x$Subtype = 'C0_C2';
dat4 = rbind( dat.A, dat.x );
#table( dat$Subtype, dat$Stage );
detach( dat3 );
attach( dat4 );
surv4 = Surv( Days, Vital );
fit4  = survfit(  surv4 ~ Subtype );
diff4 = survdiff( surv4 ~ Subtype );
pval4 = 1 - pchisq(diff4$chisq, length(diff4$n) - 1);
fit.coxph4 = coxph( surv4~Subtype, data=dat4 );
hr4 = exp( coef(fit.coxph4) )['SubtypeC1'];

subtype.cat = count(dat4, 'Subtype');
subtype.cat.num = nrow( subtype.cat );

plot( fit4, col=c("blue", "red"), ylim=c(0.4, 1),
                xlab="Time (days)", ylab="Survival Probability", lwd=4, mark.time=T,
                main=paste0("P=", sprintf("%.4g", pval4), ", HR=", sprintf("%.4f", hr4)) );
legend('bottomleft',paste(subtype.cat[,1], " (n=", subtype.cat[,2], ")", sep=""),
                        col=c("blue","red"), lty=rep(1,2), bty='n', cex=1.1);

