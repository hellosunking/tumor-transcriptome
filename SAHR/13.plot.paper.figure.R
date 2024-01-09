
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
##################################################
source("./plot.papper.surv.R")
pdf("SAHR.all.cancer.old.model.pdf",width = 16,height = 8)
for (cancer in c("BLCA","BRCA","COAD","HCC","HNSC.Larynx","HNSC.Tongue","KIRC","KIRP","LUAD","LUSC","STAD",
                 "THCA")){
train=read.table(paste0(cancer,".top60.train.matrix.SAHR"), head=T );
test=read.table(paste0(cancer,".top60.test.matrix.SAHR"), head=T );
#filename=paste0(cancer,".survival.paper.figure.pdf")
#pdf(filename,width = 11,height = 7)
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,4))
plot(density(train$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer," in training"))
survplot(train,cancer,"training")
survplot_stage(train,cancer,"train")
survplot_stage_2(train,cancer,"train")

plot(density(test$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer," in testing"))
survplot(test,cancer,"testing")
survplot_stage(test,cancer,"testing")
survplot_stage_2(test,cancer,"testing")

par(opar)
}

for (cancer in c("UCEC")){
  train=read.table(paste0(cancer,".Endometrioid.top60.train.matrix.SAHR"), head=T );
  test=read.table(paste0(cancer,".Endometrioid.top60.test.matrix.SAHR"), head=T );
  #filename=paste0(cancer,".survival.paper.figure.pdf")
  #pdf(filename,width = 11,height = 7)
  opar <- par(no.readonly=TRUE)
  par(mfrow=c(2,4))
  plot(density(train$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer," in training"))
  survplot(train,cancer,"training")

  plot(density(test$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer," in testing"))
  survplot(test,cancer,"testing")
  par(opar)
}

dev.off()
