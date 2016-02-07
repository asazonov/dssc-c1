file_loc = "~/Dropbox/Work/hackathon/"
out_loc = "~/Dropbox/Work/hackathon/dssc-c1/data_out/"
setwd(out_loc)

source("https://bioconductor.org/biocLite.R")

library(useful)
library(devtools)
library(DESeq)
library(dynamicTreeCut)
library(Rtsne)
library(diffusionMap)
#library(BASiCS)
library(topGO)
#library(scde)
library(dplyr)
library(ctc)



cellcycle = read.csv(paste0(file_loc, "dssc-c1/cellcycle.csv"), row.names = 1)
zebra = read.table(sep = "\t", file = paste0(file_loc, "dssc-c1/zebrafish.txt"))
ear = read.table(sep = '\t', file = paste0(file_loc, "dssc-c1/ear.txt"), header = T, row.names = 1)
mouse = read.table(sep = "\t", file = paste0(file_loc, "mouse.txt")) 
brain = read.csv(file = paste0(file_loc, "brain.csv"), header = T, row.names = 1) 


data = round(ear) #### MAKE THIS POINT TO THE UPLOADED FILE
data[1:5, 1:5] # as a check - put this on the page?
#step 1: size factors via DESeq
counts = newCountDataSet(data, conditions = names(data))
counts = estimateSizeFactors(counts)
counts.matrix = counts(counts)
#Put the counts on the common scale if it made sizefactors:
if(!is.na(sizeFactors(counts)[1])){
  counts.adj = t(t(counts.matrix)/sizeFactors(counts))
  } else {counts.adj = t(t(counts.matrix))}

#step 2: high var genes (dimension reduction)
coef_var2 = function(vec_gene){
  return(var(vec_gene)/(mean(vec_gene)^2))
}

counts.cv2 = apply(counts.adj, MARGIN = 1, FUN = coef_var2)
counts.avg = apply(counts.adj, MARGIN = 1, FUN = mean)
counts.var = apply(counts.adj, MARGIN = 1, FUN = var)

#in the paper they removed these low ones from the regression
lower_limit = quantile(counts.avg)[2]

cv2.fit = counts.cv2[which(counts.avg > max(0.5, lower_limit))]
avg.fit = counts.avg[which(counts.avg > max(0.5, lower_limit))]


if(is.na(sizeFactors(counts)[1])){size_factors = rep(1, dim(counts.matrix)[2])} else {size_factors = sizeFactors(counts)}

#data fit using gamma family:
model = glm(cv2.fit ~ I(1/avg.fit), family = Gamma(link = 'identity'))# native (slower) implementation
grad = model$coefficients[2]; a1 = grad
int = model$coefficients[1]; a0 = int
#formalised estimation of genes (no spike)
df = dim(counts.matrix)[2]-1
psia1theta = mean(1/size_factors) + a1
minBiolDisp = .5 ^ 2
cv2th = a0 + minBiolDisp + a0 * minBiolDisp
testDenom = ( counts.avg * psia1theta + counts.avg^2 * cv2th ) / ( 1 + cv2th/df )
p = 1 - pchisq( counts.var * (df-1) / testDenom, df-1 )
padj = p.adjust(p, "BH")
sig = padj<0.1
high.var = names(sig)[which(sig)]
x_vals = 10^seq(from = -3, to = 5, length.out=1000)

pdf("highvarplot.pdf", width = 10, height = 10)
plot(x= counts.avg, y = counts.cv2, pch='.', log='xy', col = ifelse(names(counts.avg)%in%high.var , 'red', 'grey'))
lines(x_vals, grad/x_vals + int, col = 'blue', lwd=1)
lines(x_vals, ((grad)/x_vals + int) * qchisq(.975, df)/df, col = 'red', lwd=1, lty='dashed')
lines(x_vals, ((grad)/x_vals + int) * qchisq(.025, df)/df, col = 'red', lwd=1, lty='dashed')
dev.off()
counts.sig = counts.adj[which(rownames(counts.matrix)%in%high.var), ]




#step 3: clustering and identification
spearman.base = cor(counts.sig); spearman = (-spearman.base+1)/2; spearman[which(is.na(spearman))]=0.25 #JANKY FIX!!
cluster = hclust(as.dist(spearman), method = 'average')
clust.labels = cutreeDynamic(cluster, method = "hybrid", minClusterSize = 10, deepSplit = 0, distM = spearman)
summary(clust.labels)
rownames(clust.labels) = names(counts.sig)
hist(clust.labels, breaks = length(unique(clust.labels)))
plot(cluster, labels = clust.labels)

corrs = spearman.base
cells=rownames(corrs)

top_only = function(x) if (abs(x)>=0.8) x else 0

corrs10 = apply(corrs,1:2,top_only)

corrs10.df<- cbind(cell = cells, corrs10)
write.csv(corrs10.df, "correlation-new.csv",row.names = FALSE,quote=FALSE)



#step 4: try dimension reduction steps: PCA/t-SNE
tsne = Rtsne(X=t(counts.sig), verbose = T, initial_dims = dim(counts.sig)[1])
plot(tsne$Y, col = clust.labels)



pca = prcomp(t(counts.sig))
plot(x=pca$x[,1], y = pca$x[,2], col = clust.labels)
plot(x=pca$x[,1], y = pca$x[,3], col = clust.labels)
plot(x=pca$x[,2], y = pca$x[,3], col = clust.labels)
plot(x=pca$x[,2], y = pca$x[,4], col = clust.labels)


dm = diffuse(dist(t(counts.sig)))
plot(dm, color = clust.labels)



##### 
#cluster analysis

#scde time (see vignette)




#####
#data saving
comps = pca$x
remove_outliers = function(vector){
  limits = quantile(vector, c(.95,.05))
  vector[which(vector<limits[2])] = NA
  vector[which(vector>limits[1])] = NA
  return(vector)
}
comps.outliersgone = apply(comps, 2, remove_outliers)
write.csv(comps.outliersgone, "princomps_outliersgone.csv")
write.csv(comps, "princomps_full.csv")

vargene_record = rbind(counts.cv2, counts.avg, as.integer(sig))
rownames(vargene_record)[3] = "signif_var"

vargene_record=cbind(gene = rownames(vargene_record), vargene_record)
vargene_record = vargene_record[, -which(is.na(vargene_record[3,]))]
write.table(x = t(vargene_record), file="vargenes.csv",quote=FALSE, col.names=FALSE,sep=",")

write.csv(spearman.base, 'correlation.csv')

clust.df = data.frame(tsne_x=tsne$Y[,1], tsne_y= tsne$Y[,2], row.names = colnames(counts.matrix), cluster = clust.labels)
write.csv(clust.df, "clust.csv")

write.table(hc2Newick(cluster), file="newick.tab",row.names=FALSE,col.names=FALSE)
