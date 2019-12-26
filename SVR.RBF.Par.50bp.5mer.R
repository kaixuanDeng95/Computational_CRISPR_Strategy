#! /usr/bin/env Rscript

#LASSO.10fold<-function(){
library(e1071)

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
names(y)<-rownames(sgRNA.1130.10fold)
sgRNA.top20.id<-rownames(sgRNA.1130[1:20,])

#PCC.SVR<-matrix(,5,4)
#rownames(PCC.SVR)<-c("500bp","200bp","100bp","50bp","25bp")
#colnames(PCC.GaussPR)<-c("rbfdot","polydot","vanilladot","tanhdot","laplacedot","besseldot","anovadot","splinedot")
#colnames(PCC.SVR)<-c("linear","polynomial","radial","sigmoid")
rbf.par=c( c(0.00001,0.00002,0.00005,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.1,0.2,0.5,1), seq(5,100,20),seq(100,1000,200));
#50bp
for(j in 4:4){
    PCC.SVR.RBF.top20<-matrix(,2,length(rbf.par))
    PCC.SVR.RBF.top20[1,]<-rbf.par
    PCC.SVR.RBF.all<-matrix(,2,length(rbf.par))
    PCC.SVR.RBF.all[1,]<-rbf.par
    y.test.pred.total_list<-list()
    x<-x.list[[j]]
    cat("j=",j)
    for(m in 1:length(rbf.par)){
    #for(m in 1:2){
        cat("m=",m)
        y.test.total <- c()
        y.test.pred.total<-c()
        tmp2<-c()
        SVR.kernel="radial"
        rbf.kpar=rbf.par[m]
        for(i in 1:10){
#          cat("\r",i)
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
}
save(PCC.SVR.RBF.top20,PCC.SVR.RBF.all,y.test.pred.total_list,y.test.real.total,file="SVR_RBF_PCC_5mer_50bp.Rdata")
#}