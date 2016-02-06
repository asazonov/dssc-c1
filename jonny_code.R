source("https://bioconductor.org/biocLite.R")
library(useful)
library(devtools)
library(DESeq)
library(dynamicTreeCut)
library(Rtsne)


bertie = read.csv("~/Dropbox/Work/hackathon/dssc-c1/bertie.csv")
zebra = read.table(sep = "\t", file = "~/Dropbox/Work/hackathon/dssc-c1/zebrafish.txt")
ear = read.table(sep = '\t', file = "~/Dropbox/Work/hackathon/dssc-c1/ear.txt", header = T); rownames(ear) = ear[,1]; ear = ear[, -1]
#mouse = read.table(sep = "\t", file = "~/Dropbox/Work/hackathon/dssc-c1/mouse.txt")

data = round(ear)
#step 1: size factors
counts = newCountDataSet(data, conditions = names(data))
counts = estimateSizeFactors(counts)
counts.matrix = counts(counts)
#Put the counts on the common scale:
counts.adj = t(t(counts.matrix)/sizeFactors(counts))
#sanity check
counts.matrix[1:5,1:5]/counts.adj[1:5, 1:5]

#step 2: high var genes (dimension reduction)
coef_var2 = function(vec_gene){
  return(var(vec_gene)/(mean(vec_gene)^2))
}

counts.cv2 = apply(counts.adj, MARGIN = 1, FUN = coef_var2)
counts.avg = apply(counts.adj, MARGIN = 1, FUN = mean)
counts.var = apply(counts.adj, MARGIN = 1, FUN = var)

#in the paper they removed these low ones from the regression
cv2.fit = counts.cv2[which(counts.avg > 10)]
avg.fit = counts.avg[which(counts.avg > 10)]


#data fit using gamma family:
model = glm(cv2.fit ~ I(1/avg.fit), family = Gamma(link = 'identity'))# native (slower) implementation
grad = model$coefficients[2]; a1 = grad
int = model$coefficients[1]; a0 = int
#formalised estimation of genes (no spike)
df = dim(counts.matrix)[2]-1
psia1theta = mean(1/sizeFactors(counts)) + a1
minBiolDisp = .5 ^ 2
cv2th = a0 + minBiolDisp + a0 * minBiolDisp
testDenom = ( counts.avg * psia1theta + counts.avg^2 * cv2th ) / ( 1 + cv2th/df )
p = 1 - pchisq( counts.var * (df-1) / testDenom, df-1 )
padj = p.adjust(p, "BH")
sig = padj<0.1
high.var = names(sig)[which(sig)]
plot(x= counts.avg, y = counts.cv2, pch='.', log='xy', col = ifelse(names(counts.avg)%in%high.var , 'red', 'grey'))
lines(x_vals, grad/x_vals + int, col = 'blue', lwd=1)
counts.sig = counts.adj[which(rownames(counts.matrix)%in%high.var), ]



#step 3: clustering and identification
spearman = cor(counts.sig); spearman = (-spearman+1)/2
cluster = hclust(as.dist(spearman), method = 'average')
clust.labels = cutreeDynamic(cluster, method = "hybrid", minClusterSize = 10, deepSplit = 2, distM = spearman)
summary(clust.labels)
hist(clust.labels, breaks = length(unique(clust.labels)))
plot(cluster, labels = clust.labels)

#step 4: try dimension reduction steps: PCA/t-SNE
tsne = Rtsne(X=t(counts.sig), verbose = T, initial_dims = dim(counts.sig)[1])
plot(tsne$Y, col = clust.labels)

pca = prcomp(t(counts.sig))
plot(x=pca$x[,1], y = pca$x[,2], col = clust.labels)
plot(x=pca$x[,1], y = pca$x[,3], col = clust.labels)
plot(x=pca$x[,2], y = pca$x[,3], col = clust.labels)

