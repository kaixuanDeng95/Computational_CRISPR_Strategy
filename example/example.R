x=as.matrix(read.table("example_7mer.txt"))
library(randomForest)
y_pred=predict(RF.model,x)
#> y_pred
#> 4.7719043 0.1287047