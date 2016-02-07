library(Seurat)

#Read in log-space expression matrix. Has been pre-computed and normalized (see manuscript for exact details)
#The data was run in three batches (zf1, zf2, zf3), which is denoted in the column name

zfish.data=read.table("/home/dilyana/Desktop/ucl-data-science/seurat_files/seurat_files/zdata.matrix.txt",sep="\t",header=TRUE)
zfish.data=read.table("/home/dilyana/dssc-c1/bertie.csv",sep=",",header = TRUE)
z<- zfish.data[,2:dim(zfish.data)[2]]


rownames(z)<- as.character(as.vector(zfish.data[,1]))
colnames(z)<- names(zfish.data)[2:length(names(zfish.data))]
z<-tFrame(z)
#Due to the high dynamic range in RNA-seq, transform data to log-space. Not required for Seurat, but highly recommended
zfish.log=log(z+1)

#Create and setup the Seurat object. Include all genes detected in > 3 cells (expression >0.01), and all cells with > 2k genes
#Cells will be initially assigned to an identity class (grouping) based on the first field (underscore-delimited)
zf.all=new("seurat",raw.data=z) 
zf.all=setup(zf.all,project="zfish",min.cells = 3,min.genes = 2000,is.expr=0.01,names.field = 1,names.delim = "_") 

zf.all

#View the expression of genes and basic metrics across samples
vlnPlot(zf.all,c("Cdh1","Ubc"))

#Plot two cells against each other
#Set do.ident=TRUE to identify outliers by clicking on points (ESC to exit after)

par(mfrow=c(2,2))
cellPlot(zf.all,zf.all@cell.names[1],zf.all@cell.names[2],do.ident = FALSE)
cellPlot(zf.all,zf.all@cell.names[3],zf.all@cell.names[4],do.ident = FALSE)

#Plot two genes against each other, can do this in limited groups of cells
genePlot(zf.all,"MIXL1","OSR1",cex.use = 1)
genePlot(zf.all,"MIXL1","SOX3",cell.ids = which.cells(zf.all,"zf1"),cex.use = 1)

#Genes placed into 20 bins based on X-axis (average expression). Y-axis is within-bin z-score of log(Variance/mean). 
zf.all=mean.var.plot(zf.all,y.cutoff = 2,do.plot=TRUE,x.low.cutoff=0.25,x.high.cutoff=7,fxn.x = expMean,fxn.y=logVarDivMean)

#remove markers which primarily define batch
markers.remove=batch.gene(zf.all,idents.use = c("zf1","zf2","zf3"),genes.use=zf.all@var.genes,auc.cutoff = 0.7)
zf.all@var.genes=zf.all@var.genes[!(zf.all@var.genes%in%markers.remove)]
length(zf.all@var.genes)

#run pca and examine results

#To demonstrate, run PCA on all genes (note that the result isn't great, PC1/2 contain little known biology)
zf.all=pca(zf.all,pc.genes = rownames(zf.all@data),pcs.print = 3,genes.print = 5)

#Instead, run PCA on zf.all@var.genes (default value for pc.genes).
zf.all=pca(zf.all,pcs.print = 3,genes.print = 5)

#If desired, project PCA structure across entire dataset (so all genes are included)
zf.all=project.pca(zf.all, do.print = FALSE) 


#visualize top genes for each PC, use.full selects the projected PCA. See also print.pca
#Explicit thanks to Olga Botvinnik, who showed me this visualization in her Flotilla package
viz.pca(zf.all,pcs.use = 1:3,num.genes = 10,use.full = TRUE,nCol = 3)

#Identify EVL cells (see Manuscript for further detail)
plot(zf.all@pca.rot[,1],zf.all@pca.rot[,2],pch=16,xlab="PC1",ylab="PC2")

x=seq(-0.2,0.2,.01); lines(x,-x*0.5-0.04,lwd=2,lty=2,col="red")
evl.quant=zf.all@pca.rot[,1]+2*zf.all@pca.rot[,2]+0.08; names(evl.quant)=colnames(zf.all@data)
not.evl=names(evl.quant[evl.quant>0]); is.evl=names(evl.quant[evl.quant<0])

points(zf.all@pca.rot[is.evl,1],zf.all@pca.rot[is.evl,2],pch=16,col="red")

#subsetData allows you to create a new Seurat object with a subset of the data 
zf=subsetData(zf.all,cells.use = not.evl)

zf

#So we can move through the tutorial in parts if needed
#Note that you can now send/e-mail/etc. the object file, enabling easy collaboration (receiver uses load)
save(zf,file="~/seurat_files/output_part1.Robj")


#Identify genes to use for building gene expression models
#load previous data object, enables you to start tutorial from Part 2.
load("~/seurat_files/output_part1.Robj")

#recalculate a set of variable genes, now that EVL are removed
zf <- mean.var.plot(zf, y.cutoff = 2, do.plot=FALSE, x.low.cutoff=1, x.high.cutoff=7, fxn.x = expMean, fxn.y=logVarDivMean, set.var.genes = TRUE)
markers.remove=batch.gene(zf,idents.use = c("zf1","zf2","zf3"),genes.use=zf@var.genes)
zf@var.genes=zf@var.genes[!(zf@var.genes%in%markers.remove)]

#redo the PCA on the variable genes.
zf <- pca(zf, do.print = FALSE)

#Run a 'random' PCA 1,000 times - scrambling a random 2.5% of the data each time
#This enables us to identify statistically significant PCs (in this case, 1:3), and genes with significant PC scores
zf <- jackStraw(zf, num.replicate=1000, prop.freq=0.025)
jackStrawPlot(zf)


#Seurat - Guided Zebrafish Tutorial - Part 2


