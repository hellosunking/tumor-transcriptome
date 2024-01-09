argv = commandArgs(T);

source("/lustre/home/dxhu/project/project4_desurv/paper_figure/all.program/01.SAHR/plot.papper.surv.R")
pdf("SAHR.pdf",width = 16,height = 8)

train=read.table(argv[1], head=T );
test=read.table(argv[2], head=T );
cancer = argv[3]
opar <- par(no.readonly=TRUE)
par(mfrow=c(2,4))
if (cancer != "UCEC"){
plot(density(train$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer," in training"))
survplot(train,cancer,"training")
survplot_stage(train,cancer,"training")
survplot_stage_2(train,cancer,"training")

plot(density(test$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer," in testing"))
survplot(test,cancer,"testing")
survplot_stage(test,cancer,"testing")
survplot_stage_2(test,cancer,"testing")
}
if(cancer == "UCEC"){
plot(density(train$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer," in training"))
survplot(train,cancer,"training")
plot(density(test$SAHR),xlab="Value", ylab="Density", lwd=2,main=paste(cancer," in testing"))
survplot(test,cancer,"testing")
}

par(opar)

dev.off()
