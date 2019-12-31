#! /usr/bin/env Rscript

#Ridge.10fold<-function(){
library(glmnet)

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
PCC.Ridge.5mer.top20<-c()
PCC.Ridge.5mer.all <- c()
y.test.pred.top20<-list()#存y.pred的top20,5个元素
y.test.pred.total_mean<-list()
for(j in 1:5){
	x<-x.list[[j]]
	cat("j=",j,'\n')
	sum_y_top20<-rep(0,20)#存储每个j对应的10个y_top20值的和
    sum_y_total<-rep(0,980)
	for(m in 1:10){
		cat("m=",m,'\t')
		y.test.total <- c()
		y.test.pred.total<-c()
		tmp2<-c()
		for(i in 1:10){
#           cat("\r",i);
			test.group.id<-which(sgRNA.1130[1:980,5]==i)
			tmp2<-c(tmp2,rownames(sgRNA.1130)[test.group.id])
			train.group.id<-which(sgRNA.1130[,5]!=i)
			x.train<-x[train.group.id,];
			x.test<-x[test.group.id,];
			y.train<-y[train.group.id];
			y.test<-y[test.group.id];
			Ridge.cv<-cv.glmnet(x.train,y.train,alpha=0);
			bestlambda<-Ridge.cv$lambda.min;
			train.model<-Ridge.cv$glmnet.fit;
			y.test.pred.Ridge<-predict(train.model,newx=x.test,s=bestlambda);
			y.test.pred.total<-c(y.test.pred.total,y.test.pred.Ridge);
			y.test.total<-c(y.test.total,y.test);
			}
		names(y.test.total)<-tmp2
		names(y.test.pred.total)<-tmp2
		sum_y_top20 <- y.test.pred.total[sgRNA.top20.id] + sum_y_top20
        sum_y_total <- y.test.pred.total+sum_y_total
		}
		y.test.pred.top20[[j]] <- sum_y_top20/10
        y.test.pred.total_mean[[j]] <- sum_y_total/10
		PCC.Ridge.5mer.top20[j]<-cor(y.test.total[sgRNA.top20.id],y.test.pred.top20[[j]])
        PCC.Ridge.5mer.all[j] <- cor(y.test.total,y.test.pred.total_mean[[j]])
	}
names(y.test.pred.total_mean) <- c("500bp","200bp","100bp","50bp","25bp")
names(y.test.pred.top20) <- c("500bp","200bp","100bp","50bp","25bp")
names(PCC.Ridge.5mer.top20) <-c ("500bp","200bp","100bp","50bp","25bp")
names(PCC.Ridge.5mer.all) <-c ("500bp","200bp","100bp","50bp","25bp")

save(y.test.total,y.test.pred.total_mean,PCC.Ridge.5mer.top20,PCC.Ridge.5mer.all,file="PCC.Ridge.5mer.Rdata")
#}