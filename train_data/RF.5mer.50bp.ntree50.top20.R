#! /usr/bin/env Rscript
#RF.10fold<-function(){

library(randomForest)
x.5mer.50bp<-as.matrix(read.table(file="p53.1130sgRNA_50bp_5mer.txt"))
x<-x.5mer.50bp
(load("sgRNA.1130.Rdata"))
sgRNA.1130<-sgRNA.1130.10fold
y<-sgRNA.1130$Z_score
names(y)<-rownames(sgRNA.1130.10fold)
sgRNA.top20.id<-rownames(sgRNA.1130[1:20,])
PCC.RF.5mer.50bp.top20<-matrix(,1,50)
PCC.RF.5mer.50bp.all<-matrix(,1,50)
for(m in 1:50){
	colnames(PCC.RF.5mer.50bp.top20)[m]=paste("mtry=",10+10*(m-1),sep="")
	colnames(PCC.RF.5mer.50bp.all)[m]=paste("mtry=",10+10*(m-1),sep="")
}

ntree<-50#ntree = 50,100,...,450,500
y.test.pred.total<-list()
	for(m in 1:50){y.test.pred.total[[m]]<-c(1)}

#for(j in 1:11){
#	ntree<-500+50*(j-1)
#	cat("ntree=",ntree)
for(m in 1:50){
	mtry<-10+10*(m-1)
	cat("mtry=",mtry)
	y.test.total<-c()
	tmp<-c()
	for(i in 1:10){
		cat("\r",i);
		test.group.id<-which(sgRNA.1130[1:980,5]==i)
		tmp<-c(tmp,rownames(sgRNA.1130)[test.group.id])
		train.group.id<-which(sgRNA.1130[,5]!=i)
		x.train<-x[train.group.id,]
		x.test<-x[test.group.id,]
		y.train<-y[train.group.id]
		y.test<-y[test.group.id]
		RF.cv<-randomForest(x.train,y.train,ntree=ntree,mtry=mtry)
		y.test.pred.RF<-predict(RF.cv,x.test)
		y.test.pred.total[[m]]<-c(y.test.pred.total[[m]],y.test.pred.RF)
		y.test.total<-c(y.test.total,y.test)
	}
	y.test.pred.total[[m]]<-y.test.pred.total[[m]][-1]
	names(y.test.total)<-tmp
	names(y.test.pred.total[[m]])<-tmp
	PCC.RF.5mer.50bp.top20[1,m]<-cor(y.test.total[sgRNA.top20.id],y.test.pred.total[[m]][sgRNA.top20.id])
	PCC.RF.5mer.50bp.all[1,m]<-cor(y.test.total,y.test.pred.total[[m]])
}
save(PCC.RF.5mer.50bp.top20,PCC.RF.5mer.50bp.all,y.test.total,y.test.pred.total,file="PCC.RF.5mer.50bp.50ntree.Rdata")
#}