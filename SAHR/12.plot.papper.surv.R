library( plyr     );
library( survival );
library( ggbeeswarm );

mytheme <- theme(axis.title=element_text(face="bold", size=14,colour = 'black'), #
                 #axis.title.y = element_text(margin = margin(r=-2,unit = "cm")),
                 plot.title = element_text(face="bold", size=14,colour = 'black',hjust = 0.5),
                 axis.text=element_text(size=12,colour = 'black'), #
                 axis.line = element_line(size=0.5, colour = 'black'), #
                 panel.background = element_blank(),
                 #plot.background = element_blank(),
                 legend.key = element_blank() #关闭图例边框
                 
)

###############################################################################
survplot <- function(raw,i,m){
  #raw=read.table(paste0("HCC",".surv.summary.top60.train.matrix.SAHR"), head=T );
  raw=data.frame(raw)
  dat=raw[ order(raw$SAHR), ];
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
hr=exp(coef(fit.coxph));        ## hazard ratio
plot( fit,col=c('blue', 'red'), xlab="Time (days)", ylab="Survival Probability", 
      lwd=2, mark.time=T, bty="l", font.lab=2,cex.lab=1.3,#cex.lab=1.6,cex.main=1.8,
      main=paste0(i,"\n","(",m," dataset)")
      );
legend('bottomleft', c(paste("SAHR- (N=", nrow(neg), ")", sep=""),
                       paste("SAHR+ (N=", nrow(pos), ")", sep="")),
       col=c('blue', 'red'), lty=c(1,1), bty='n', cex=1.1);
  text(max(fit$time),0.9,adj=1,c(paste("P=", ifelse(pval < 1e-10, "<1e-10", ifelse(pval < 0.01, sprintf("%.1e", pval),sprintf("%.2g", pval))), ", HR=", ifelse(hr<100000,sprintf("%.2f", hr),sprintf("%.1e", hr) ),sep="")));
detach( SAHR.labeled );
}

survplot_stage <- function(raw,i,m){
  s1.neg=subset(raw, (raw$Stage=='Stage_I' | raw$Stage=='Stage_II') & raw$SAHR<0);
  s1.pos=subset(raw, (raw$Stage=='Stage_I' | raw$Stage=='Stage_II') & raw$SAHR>=0);
  s2.neg=subset(raw, (raw$Stage=='Stage_III' | raw$Stage=='Stage_IV') & raw$SAHR<0);
  s2.pos=subset(raw, (raw$Stage=='Stage_III' | raw$Stage=='Stage_IV') & raw$SAHR>=0);
  
  code.1=rep("C1", nrow(s1.neg));
  code.2=rep("C2", nrow(s1.pos));
  code.3=rep("C3", nrow(s2.neg));
  code.4=rep("C4", nrow(s2.pos));
  
  s1.neg.code=cbind( s1.neg, code.1);
  colnames( s1.neg.code )=c( colnames(raw), "Code" );
  s1.pos.code=cbind( s1.pos, code.2 );
  colnames( s1.pos.code )=c( colnames(raw), "Code" );
  
  s2.neg.code=cbind( s2.neg, code.3 );
  colnames( s2.neg.code )=c( colnames(raw), "Code" );
  s2.pos.code=cbind( s2.pos, code.4 );
  colnames( s2.pos.code )=c( colnames(raw), "Code" );
  ## Stagge I
  SAHR.labeled=rbind(s1.neg.code, s1.pos.code);
  attach( SAHR.labeled );
  
  surv = Surv(Days, Vital);
  fit  = survfit(  surv~Code, data=SAHR.labeled );
  diff = survdiff( surv~Code, data=SAHR.labeled );
  pval = 1 - pchisq(diff$chisq, length(diff$n)-1 );
  
  fit.coxph = coxph( surv~Code, data=SAHR.labeled );
  hr=exp(coef(fit.coxph));        ## hazard ratio(s)
  plot( fit, col=c('blue', 'purple'), xlab="Time (days)", ylab="Survival Probability",
        lwd=2, mark.time=T,bty="l",font.lab=2,cex.lab=1.3,#cex.lab=1.6,cex.main=1.8,
        main=paste(i, "\n(early Stage)"));
  legend( 'bottomleft', c(paste("SAHR- (N=", nrow(s1.neg), ")", sep=""),
                          paste("SAHR+ (N=", nrow(s1.pos), ")", sep="") ),
          col=c('blue', 'purple'), lty=rep(1,2), bty='n', cex=1.1);
  text(max(fit$time),0.9,adj=1,c(paste("P=", ifelse(pval < 1e-10, "<1e-10", ifelse(pval < 0.01, sprintf("%.1e", pval),sprintf("%.2g", pval))), ", HR=",ifelse(hr<100000,sprintf("%.2f", hr),sprintf("%.1e", hr)), sep="")));
  detach( SAHR.labeled );
}

survplot_stage_2 <- function(raw,i,m){
  s1.neg=subset(raw, (raw$Stage=='Stage_I' | raw$Stage=='Stage_II') & raw$SAHR<0);
  s1.pos=subset(raw, (raw$Stage=='Stage_I' | raw$Stage=='Stage_II') & raw$SAHR>=0);
  s2.neg=subset(raw, (raw$Stage=='Stage_III' | raw$Stage=='Stage_IV') & raw$SAHR<0);
  s2.pos=subset(raw, (raw$Stage=='Stage_III' | raw$Stage=='Stage_IV') & raw$SAHR>=0);
  
  code.1=rep("C1", nrow(s1.neg));
  code.2=rep("C2", nrow(s1.pos));
  code.3=rep("C3", nrow(s2.neg));
  code.4=rep("C4", nrow(s2.pos));
  
  s1.neg.code=cbind( s1.neg, code.1);
  colnames( s1.neg.code )=c( colnames(raw), "Code" );
  s1.pos.code=cbind( s1.pos, code.2 );
  colnames( s1.pos.code )=c( colnames(raw), "Code" );
  
  s2.neg.code=cbind( s2.neg, code.3 );
  colnames( s2.neg.code )=c( colnames(raw), "Code" );
  s2.pos.code=cbind( s2.pos, code.4 );
  colnames( s2.pos.code )=c( colnames(raw), "Code" );
  ## Stage II
  SAHR.labeled=rbind(s2.neg.code, s2.pos.code);
  attach( SAHR.labeled );
  
  surv = Surv(Days, Vital);
  fit  = survfit(  surv~Code, data=SAHR.labeled );
  diff = survdiff( surv~Code, data=SAHR.labeled );
  pval = 1 - pchisq(diff$chisq, length(diff$n)-1 );
  
  fit.coxph = coxph( surv~Code, data=SAHR.labeled );
  hr=exp(coef(fit.coxph));        ## hazard ratio(s)
  
  plot( fit, col=c('red', 'red4'), xlab="Time (days)", ylab="Survival Probability",
        lwd=2, mark.time=T,bty="l",font.lab=2,cex.lab=1.3,#cex.lab=1.6,cex.main=1.8,
        main=paste(i, "\n(late Stage)"));
  legend( 'bottomleft', c(paste("SAHR- (N=", nrow(s2.neg), ")", sep=""),
                          paste("SAHR+ (N=", nrow(s2.pos), ")", sep="") ),
          col=c('red', 'red4'), lty=rep(1,2), bty='n', cex=1.1);
  text(max(fit$time),0.9,adj=1,c(paste("P=", ifelse(pval < 1e-10, "<1e-10", ifelse(pval < 0.01, sprintf("%.1e", pval),sprintf("%.2g", pval))), ", HR=", ifelse(hr<100000,sprintf("%.2f", hr),sprintf("%.1e", hr)), sep="")));
  detach( SAHR.labeled );
 
}
