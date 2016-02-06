source("https://bioconductor.org/biocLite.R")
library(DESeq)
library(statmod)
library(biomaRt)
library(fpc)
library(dynamicTreeCut)
###### 
#DATA IN
vargenes.reference = read.table("~/Dropbox/Work/Cambridge/Marioni/high_var_genes_all.tsv")
vargenes.reference = as.vector(var.genes[,1], mode = "character")
meta = read.table("~/Dropbox/Work/Cambridge/Marioni/meta", header = T)
counts.in = read.table("~/Dropbox/Work/Cambridge/Marioni/counts.txt", header = T)
ensembl = useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), mart = ensembl, values=rownames(counts.in)[1:5], filters = 'ensembl_gene_id')
counts = newCountDataSet(counts.in, conditions = names(counts.in))
counts = estimateSizeFactors(counts)
counts.matrix = counts(counts)
#Put the counts on the common scale:
counts.adj = t(t(counts.matrix)/sizeFactors(counts))
#sanity check
counts.matrix[1:5,1:3]/counts.adj[1:5, 1:3]
######
#GET HIGHLY VARIABLE GENES
#simply calculates squared coef of variation
coef_var2 = function(vec_gene){
  return(var(vec_gene)/(mean(vec_gene)^2))
}
#this is how many of names1 are in names2
fraction_shared = function(names1, names2){
  test = names1 %in% names2
  return(sum(test)/length(test))
}

counts.cv2 = apply(counts.adj, MARGIN = 1, FUN = coef_var2)
counts.avg = apply(counts.adj, MARGIN = 1, FUN = mean)
counts.var = apply(counts.adj, MARGIN = 1, FUN = var)

#in the paper they removed these low ones from the regression
cv2.fit = counts.cv2[which(counts.avg > 10)]
avg.fit = counts.avg[which(counts.avg > 10)]


#data fit using gamma family:
model = glmgam.fit(cbind(a0=1, a1=1/avg.fit), cv2.fit)
model = glm(cv2.fit ~ I(1/avg.fit), family = Gamma(link = 'identity'))# native (slower) implementation
grad = model$coefficients[2]; a1 = grad
int = model$coefficients[1]; a0 = int

#plot of fit
plot(x= counts.avg, y = counts.cv2, pch='.', log='xy')
x_vals = 10^seq(from = -2, to = 5, length.out=1000)
lines(x_vals, grad/x_vals + int, col = 'red', lwd=3)
df=dim(counts.matrix)[2]-1
xi = mean(1/sizeFactors(counts))
lines(x_vals, ((xi+grad)/x_vals + int) * qchisq(.975, df)/df, col = 'red', lwd=2, lty='dashed')
lines(x_vals, ((xi+grad)/x_vals + int) * qchisq(.025, df)/df, col = 'red', lwd=2, lty='dashed')

#my janky count
residuals = counts.cv2-((grad/counts.avg + int))
residuals[which(residuals<0)]=NA
resid.test = residuals[!is.na(residuals)]
#want to test if these are >0, errors distributed as chisq/df
tests = p.adjust(pchisq(q=residuals, df = dim(counts.matrix)[2]-1, lower.tail = T), method = "BH")/(dim(counts.matrix)[2]-1)
minimum = qchisq(p=0.1, df = (dim(counts.matrix)[2]-1))
candidates = counts.cv2
resid = residuals(model)
high.var = names(test[which(test==T)])


fraction_shared(vargenes.reference, high.var)
plot(x= counts.avg, y = counts.cv2, pch='.', log='xy', col = ifelse(test, 'red', 'grey'))


#trying to formalise as per Brenneke: (This isn't quite working but these genes will do for now)
psia1theta = mean(1/sizeFactors(counts)) + a1
minBiolDisp = .5 ^ 2
cv2th = a0 + minBiolDisp + a0 * minBiolDisp
testDenom = ( counts.avg * psia1theta + counts.avg^2 * cv2th ) / ( 1 + cv2th/df )
p = 1 - pchisq( counts.var * (df-1) / testDenom, df-1 )
padj = p.adjust(p, "BH")
sig = padj<0.1
table(sig)

high.var = names(sig)[which(sig)]
length(high.var)
fraction_shared(high.var, vargenes.reference)
plot(x= counts.avg, y = counts.cv2, pch='.', log='xy', col = ifelse(names(counts.avg)%in%high.var , 'red', 'grey'))
lines(x_vals, grad/x_vals + int, col = 'blue', lwd=1)

#get the counts matrix for the significant genes
counts.sig = counts.adj[which(rownames(counts.matrix)%in%vargenes.reference), ]

#random mortgage inefficiency calcualtor
mortgage = function(amount, rate, payment){
  amount_paid = 0
  total = amount
  years = 0
  while(total>0){
    if(total*rate - total > payment){return("You'll never win")}
    years = years+1
    print(paste("you owe", total))
    else if(total>payment){
      total = total - payment
      total = total*rate
      amount_paid = amount_paid + payment
    }
    else{amount_paid = amount_paid + total; total = 0}
  }
  print(paste("It took you", years, "years, and you paid this fraction of the house value:"))
  return(amount_paid/amount)
}
###### 
#Find the clusters
# distance = dist(counts.sig)
# plot(hclust(distance), labels = F)
spearman = cor(counts.sig); spearman = (-spearman+1)/2
cluster = hclust(as.dist(spearman), method = 'average')
clust.labels = cutreeDynamic(cluster, method = "hybrid", minClusterSize = 10, deepSplit = 2, distM = spearman, cutHeight = 0.02)
summary(clust.labels)
hist(clust.labels, breaks = length(unique(clust.labels)))
plot(cluster, labels = clust.labels)



###### 
#Seurat approach
library(Seurat)

seurat.matrix = log(counts.in + 1)
nbt=new("seurat",raw.data=seurat.matrix)
# Take all genes in > 3 cells, all cells with > 1k genes, use an expression threshold of 1
# Cell type is encoded in the second _ field, will be stored in nbt@ident and also placed in the "orig.ident" field of object@data.info 
nbt=setup(nbt,project="NBT",min.cells = 3,names.field = 3,names.delim = "_",min.genes = 1000,is.expr=1)
par(mfrow=c(2,2))
cellPlot(nbt,nbt@cell.names[1],nbt@cell.names[2],do.ident = FALSE)
cellPlot(nbt,nbt@cell.names[3],nbt@cell.names[4],do.ident = FALSE)
par(mfrow = c(1,1))
nbt=mean.var.plot(nbt,y.cutoff = 2,x.low.cutoff = 2,fxn.x = expMean,fxn.y = logVarDivMean)
length(nbt@var.genes)

# Run a PCA using the variable genes as input (to change the input gene set, use the pc.genes argument)
# For example, try running PCA on all genes - i.e. , nbt=pca(nbt,pc.genes=rownames(nbt@data)) - which does not perform as well
nbt=pca(nbt,do.print=FALSE)
pca.plot(nbt,1,2,pt.size = 2)
raw = nbt@pca.x
