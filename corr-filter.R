library(dplyr)

corrs<-read.csv("/home/dilyana/dssc-c1/data_out/correlation.csv", skip=1, header=FALSE)
cells=as.vector(corrs[,1])
corrs<-corrs[,-1]
rownames(corrs)<-cells
colnames(corrs)<-cells

y = function(x) if (abs(x)>=0.8) x else 0

corrs10 = apply(corrs,1:2,y)

corrs10.df<- cbind(cell = cells, corrs10)
write.csv(corrs10.df, "correlation-new.csv",row.names = FALSE,quote=FALSE)

