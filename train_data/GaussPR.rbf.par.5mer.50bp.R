#! /usr/bin/env Rscript
#LASSO.10fold<-function(){
library(kernlab)

x.5mer.500bp<-as.matrix(read.table(file="1130sgRNA_500bp_5mer.txt"))
x.5mer.200bp<-as.matrix(read.table(file="1130sgRNA_200bp_5mer.txt"))
x.5mer.100bp<-as.matrix(read.table(file="1130sgRNA_100bp_5mer.txt"))
x.5mer.50bp<-as.matrix(read.table(file="1130sgRNA_50bp_5mer.txt"))
x.5mer.25bp<-as.matrix(read.table(file="1130sgRNA_25bp_5mer.txt"))
x.list<-list()
x.list[[1]]<-x.5mer.500bp
x.list[[2]]<-x.5mer.200bp
x.list[[3]]<-x.5mer.100bp
x.list[[4]]<-x.5mer.50bp
x.list[[5]]<-x.5mer.25bp

(load("sgRNA.1130.Rdata"))
sgRNA.1130<-sgRNA.1130.10fold
y<-sgRNA.1130$Z_score
sgRNA.top20.id<-rownames(sgRNA.1130[1:20,])

#PCC.GaussPR<-matrix(,5,7)
#rownames(PCC.GaussPR)<-c("500bp","200bp","100bp","50bp","25bp")
#colnames(PCC.GaussPR)<-c("rbfdot","polydot","vanilladot","tanhdot","laplacedot","besseldot","anovadot","splinedot")
#colnames(PCC.GaussPR)<-c("rbfdot","polydot","vanilladot","tanhdot","laplacedot","besseldot","anovadot")
rbf.par=c(seq(0.00001,0.01,0.0005),seq(0.02,0.08,0.02),seq(0.1,1,0.2), seq(5,95,5),seq(100,1000,50),seq(1100,5000,200))
PCC.GaussPR.50bp.RBF<-matrix(,2,length(rbf.par));
PCC.GaussPR.50bp.RBF[1,]<-rbf.par;
PCC.GaussPR.50bp.RBF.top20<-matrix(,2,length(rbf.par));
PCC.GaussPR.50bp.RBF.top20[1,]<-rbf.par;
#j=4,50bp
for(j in 4:4){
	x<-x.list[[j]]
	cat("j=",j)
	for(m in 1:length(rbf.par)){
		cat("m=",m)
		y.test.total <- c()
		y.test.pred.total<-c()
		tmp2<-c()
		GaussPR.kernel="rbfdot"
		rbf.kpar=rbf.par[m]
		for(i in 1:10){
		#	cat("\r",i);
			test.group.id<-which(sgRNA.1130[1:980,5]==i)
			tmp2<-c(tmp2,rownames(sgRNA.1130)[test.group.id])
			train.group.id<-which(sgRNA.1130[,5]!=i)
			x.train<-x[train.group.id,];
			x.test<-x[test.group.id,];
			y.train<-y[train.group.id];
			y.test<-y[test.group.id];
			train.model<-gausspr(x.train,y.train,kernel=GaussPR.kernel,kpar=list(sigma=rbf.kpar))
			y.test.pred.gausspr<-predict(train.model,x.test);
			y.test.pred.total<-c(y.test.pred.total,y.test.pred.gausspr);
			y.test.total<-c(y.test.total,y.test);
		}
		names(y.test.total)<-tmp2
		names(y.test.pred.total)<-tmp2
		PCC.GaussPR.50bp.RBF[2,m]<-cor(y.test.total,y.test.pred.total)
		PCC.GaussPR.50bp.RBF.top20[2,m] <- cor(y.test.total[sgRNA.top20.id],y.test.pred.total[sgRNA.top20.id])
	}
}
save(PCC.GaussPR.50bp.RBF,PCC.GaussPR.50bp.RBF.top20,y.test.total,y.test.pred.total,file="GaussPR_RBF_5mer_50bp_10fold_v4.Rdata")
#}