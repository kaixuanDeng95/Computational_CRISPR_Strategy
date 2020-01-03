##Lasso
library(glmnet)
x<-as.matrix(read.table(file="1130sgRNA_100bp_7mer.txt"))
(load("sgRNA.1130.Rdata"))
sgRNA.1130<-sgRNA.1130.10fold
y<-sgRNA.1130$Z_score
sgRNA.top20.id<-rownames(sgRNA.1130[1:20,])
sum_y_top20<-rep(0,20)
sum_y_total<-rep(0,980)
for(m in 1:10){
	cat("m=",m,'\t')
	y.test.total <- c()
	y.test.pred.total<-c()
	tmp2<-c()
	for(i in 1:10){
		test.group.id<-which(sgRNA.1130[1:980,5]==i)
		tmp2<-c(tmp2,rownames(sgRNA.1130)[test.group.id])
		train.group.id<-which(sgRNA.1130[,5]!=i)
		x.train<-x[train.group.id,]
		x.test<-x[test.group.id,]
		y.train<-y[train.group.id]
		y.test<-y[test.group.id]
		Lasso.cv<-cv.glmnet(x.train,y.train,alpha=1)
		bestlambda<-Lasso.cv$lambda.min
		train.model<-Lasso.cv$glmnet.fit
		y.test.pred.Lasso<-predict(train.model,newx=x.test,s=bestlambda)
		y.test.pred.total<-c(y.test.pred.total,y.test.pred.Lasso)
		y.test.total<-c(y.test.total,y.test)
		}
	names(y.test.total)<-tmp2
	names(y.test.pred.total)<-tmp2
	sum_y_top20 <- y.test.pred.total[sgRNA.top20.id] + sum_y_top20
	sum_y_total <- y.test.pred.total+sum_y_total
	}
y.test.pred.top20 <- sum_y_top20/10
y.test.pred.total_mean <- sum_y_total/10
PCC.Lasso.7mer.top20<-cor(y.test.total[sgRNA.top20.id],y.test.pred.top20)
PCC.Lasso.7mer.all <- cor(y.test.total,y.test.pred.total_mean)
save(y.test.total,y.test.pred.total_mean,PCC.Lasso.7mer.top20,PCC.Lasso.7mer.all,file="PCC.Lasso.7mer.Rdata")

##Ridge Regression
sum_y_top20<-rep(0,20)
sum_y_total<-rep(0,980)
for(m in 1:10){
	cat("m=",m,'\t')
	y.test.total <- c()
	y.test.pred.total<-c()
	tmp2<-c()
	for(i in 1:10){
		test.group.id<-which(sgRNA.1130[1:980,5]==i)
		tmp2<-c(tmp2,rownames(sgRNA.1130)[test.group.id])
		train.group.id<-which(sgRNA.1130[,5]!=i)
		x.train<-x[train.group.id,]
		x.test<-x[test.group.id,]
		y.train<-y[train.group.id]
		y.test<-y[test.group.id]
		Ridge.cv<-cv.glmnet(x.train,y.train,alpha=0)
		bestlambda<-Ridge.cv$lambda.min
		train.model<-Ridge.cv$glmnet.fit
		y.test.pred.Ridge<-predict(train.model,newx=x.test,s=bestlambda)
		y.test.pred.total<-c(y.test.pred.total,y.test.pred.Ridge)
		y.test.total<-c(y.test.total,y.test)
		}
	names(y.test.total)<-tmp2
	names(y.test.pred.total)<-tmp2
	sum_y_top20 <- y.test.pred.total[sgRNA.top20.id] + sum_y_top20
	sum_y_total <- y.test.pred.total+sum_y_total
	}
y.test.pred.top20 <- sum_y_top20/10
y.test.pred.total_mean <- sum_y_total/10
PCC.Ridge.7mer.top20<-cor(y.test.total[sgRNA.top20.id],y.test.pred.top20)
PCC.Ridge.7mer.all <- cor(y.test.total,y.test.pred.total_mean)
save(y.test.total,y.test.pred.total_mean,PCC.Ridge.7mer.top20,PCC.Ridge.7mer.all,file="PCC.Ridge.7mer.Rdata")

##RandomForest
library(randomForest)
names(y)<-rownames(sgRNA.1130.10fold)
PCC.RF.7mer.100bp.top20<-matrix(,1,50)
PCC.RF.7mer.100bp.all<-matrix(,1,50)
colnames(PCC.RF.7mer.100bp.top20)<-rep("ntree",50)
colnames(PCC.RF.7mer.100bp.all)<-rep("ntree",50)
for(m in 1:50){
	colnames(PCC.RF.7mer.100bp.top20)[m]=paste("mtry=",10+10*(m-1),sep="")
	colnames(PCC.RF.7mer.100bp.all)[m]=paste("mtry=",10+10*(m-1),sep="")
}
ntree<-50#range from 50 to 500
y.test.pred.total<-list()
for(m in 1:50){y.test.pred.total[[m]]<-c(1)}
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
	PCC.RF.7mer.100bp.top20[1,m]<-cor(y.test.total[sgRNA.top20.id],y.test.pred.total[[m]][sgRNA.top20.id])
	PCC.RF.7mer.100bp.all[1,m]<-cor(y.test.total,y.test.pred.total[[m]])
}
save(PCC.RF.7mer.100bp.top20,PCC.RF.7mer.100bp.all,y.test.total,y.test.pred.total,file="PCC.RF.7mer.100bp.50ntree.Rdata")

##SVR
library(e1071)
rbf.par=c( c(0.00001,0.00002,0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.1,0.2,0.5,1), seq(5,100,20),seq(100,1000,200))
PCC.SVR.RBF.top20<-matrix(,2,length(rbf.par))
PCC.SVR.RBF.top20[1,]<-rbf.par
PCC.SVR.RBF.all<-matrix(,2,length(rbf.par))
PCC.SVR.RBF.all[1,]<-rbf.par
y.test.pred.total_list<-list()
for(m in 1:length(rbf.par)){
	cat("m=",m)
	y.test.total <- c()
	y.test.pred.total<-c()
	tmp2<-c()
	SVR.kernel="radial"
	rbf.kpar=rbf.par[m]
	for(i in 1:10){
		test.group.id<-which(sgRNA.1130[1:980,5]==i)
		tmp2<-c(tmp2,rownames(sgRNA.1130)[test.group.id])
		train.group.id<-which(sgRNA.1130[,5]!=i)
		x.train<-x[train.group.id,];
		x.test<-x[test.group.id,];
		y.train<-y[train.group.id];
		y.test<-y[test.group.id];
		train.model<-svm(x.train,y.train,kernel=SVR.kernel,gamma=rbf.kpar);
		y.test.pred.svr<-predict(train.model,x.test);
		y.test.pred.total<-c(y.test.pred.total,y.test.pred.svr);
		y.test.total<-c(y.test.total,y.test);
	}
	names(y.test.total)<-tmp2
	names(y.test.pred.total)<-tmp2
	PCC.SVR.RBF.top20[2,m]<-cor(y.test.total[sgRNA.top20.id],y.test.pred.total[sgRNA.top20.id])
	PCC.SVR.RBF.all[2,m] <- cor(y.test.total,y.test.pred.total)
	y.test.pred.total_list[[m]] <- y.test.pred.total
	y.test.real.total <- y.test.total
}
save(PCC.SVR.RBF.top20,PCC.SVR.RBF.all,y.test.pred.total_list,y.test.real.total,file="SVR_RBF_PCC_7mer_100bp.Rdata")

##GPR
library(kernlab)
x<-as.matrix(read.table(file="1130sgRNA_100bp_7mer.txt"))
(load("sgRNA.1130.Rdata"))
sgRNA.1130<-sgRNA.1130.10fold
y<-sgRNA.1130$Z_score
sgRNA.top20.id<-rownames(sgRNA.1130[1:20,])
rbf.par=c(seq(1000,5000,200), seq(5000,16000,500));
PCC.GaussPR.100bp.RBF<-matrix(,2,length(rbf.par))
PCC.GaussPR.100bp.RBF[1,]<-rbf.par
PCC.GaussPR.100bp.RBF.top20<-matrix(,2,length(rbf.par))
PCC.GaussPR.100bp.RBF.top20[1,]<-rbf.par
for(m in 1:length(rbf.par)){
	cat("m=",m)
	y.test.total <- c()
	y.test.pred.total<-c()
	tmp2<-c()
	GaussPR.kernel="rbfdot"
	rbf.kpar=rbf.par[m]
	for(i in 1:10){
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
	PCC.GaussPR.100bp.RBF[2,m]<-cor(y.test.total,y.test.pred.total)
	PCC.GaussPR.100bp.RBF.top20[2,m] <- cor(y.test.total[sgRNA.top20.id],y.test.pred.total[sgRNA.top20.id])
}
save(PCC.GaussPR.100bp.RBF,PCC.GaussPR.100bp.RBF.top20,y.test.total,y.test.pred.total,file="GaussPR_RBF_PCC_7mer_100bp.Rdata")
