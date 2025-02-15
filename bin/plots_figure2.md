# Figure 2: Distributions of cell-type proportions across sex. 

### Panel A & C 
```{r}
in_seurat_file = "cell_type.RDS"
in_indiv_phen_file = "phen.RDS"
in_meta_file = "metadata.Rdata"
# Read in seurat object and individual metadata 

obj   <- readRDS(in_seurat_file)
phens <- readRDS(in_indiv_phen_file)
load( in_meta_file) 
obj@meta.data = metadata 

load("umap.Rdata") 
load("palette2.Rdata")
load("umap_featuresXY.Rdata") 

obj@reductions$umap@cell.embeddings = umap1 
 


m = match(levels(factor(obj$predicted.celltype.l2)), colpals[,2] )


png("umap_l2.png", height=2000, width=2000)
DimPlot(obj, 
    reduction = "umap", 
    pt.size = 2, 
    raster=FALSE, 
    group.by = "predicted.celltype.l2" , 
    cols= colpals[m,3]) + NoLegend() 
dev.off()
 


m = match(levels(factor(obj$predicted.celltype.l3)), colpals5[,11] )


png("umap_l3.png", height=2000, width=2000)
DimPlot(obj, 
    reduction = "umap", 
    pt.size = 2, 
    raster=FALSE, 
    group.by = "predicted.celltype.l3" , 
    cols= colpals5[m,12]) + NoLegend() 
dev.off()
 

 

png("umap_sex.png", height=2000, width=2000)
DimPlot(obj, 
    reduction = "umap", 
    pt.size = 2, 
    raster=FALSE, 
    group.by = "sex" , 
    cols= colsex) + NoLegend() 
dev.off()
```




### Panel B total counts + individuals 
```{r}
load("freq_sex_cell_l2.Rdata") 
load("palette2.Rdata")

m = match(colpals[,2], freq_sex_cell[,1])
freq_sex_cell = freq_sex_cell[m,]
#freq_sex_cell = freq_sex_cell[1:27,]
freq_sex_cell

pdf("freq_sex_cell.pdf", height=10) 
barplot( log10( t(t(freq_sex_cell[,c(4:2)])) ) , col=colpals[-32,3], beside=T, hor=T ) 
dev.off() 
```

### Panel D proportions  
```{r}
load(file="props_tests_cell_l2_merged_skip.Rdata") 
load(file="palette2.Rdata")

nn = 75
n =  dim(tempfrac0)[1]  
ni = nn/n
nni = 5
oi = 1:n 

m = match(colpals2[,9], rownames(tempfrac0) )
f.c = !is.na(m)
f.t = m[f.c] 
  
colpals2a = colpals2[f.c,]

 
tempfrac0   =   tempfrac0[f.t,]
tempfrac0lf =   tempfrac0lf[f.t,]
tempfrac0lm =   tempfrac0lm[f.t,]

tempfrac0lf.d   =   tempfrac0lf.d[f.t]
tempfrac0lf.d2  =   tempfrac0lf.d2[f.t]

tempfrac0lm.d   =   tempfrac0lm.d[f.t]
tempfrac0lm.d2  =   tempfrac0lm.d2[f.t]

tempfrac0l.d    =   tempfrac0l.d[f.t]
tempfrac0l.d2   =   tempfrac0l.d2[f.t]

 props0 = props0[f.t,]

xxx = 0*(props0$FDR)
xxx[(props0$FDR  < 0.05)] = "*"
xxx[(props0$FDR  < 0.01)] = "**"
xxx[(props0$FDR  < 0.001)] = "***"
xxx[xxx==0] = ""



 
plot(-10,-10, xlim=c(0,1), ylim=c(0,nn+ni), xlab="Fraction", ylab="Density" ,axes=F)
axis(1)
tt = lapply(oi, function(i) abline( h = ((i-1)*ni)  ,  col=makeTransparent(1), lty=2 ) ) 
tt = lapply(oi, function(i) polygon(tempfrac0lf.d2[[i]][,1],((i-1)*ni)+(tempfrac0lf.d2[[i]][,2]*nni/n) ,  border=NA, col= (colpals2a[i,10] ) ) ) 
tt = lapply(oi, function(i) polygon(tempfrac0lm.d2[[i]][,1],((i-1)*ni)-(tempfrac0lm.d2[[i]][,2]*nni/n) ,  border=NA, col= (colpals2a[i,10] ) ) ) 
tt = lapply(oi, function(i) segments( mean(tempfrac0lf[i,]),  ((i-1)*ni) , 
                                       mean(tempfrac0lf[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= colsex[2]  ) )

tt = lapply(oi, function(i) segments( mean(tempfrac0lm[i,]),  ((i-1)*ni) , 
                                       mean(tempfrac0lm[i,]), ((i-1)*ni) - ni/2 , lwd=2, col= colsex[1]) )


text(0.8, ((oi -1)*ni), colpals2a[,9])
 text(1, ((oi -1)*ni), xxx)


 
pdf("props_tests_cell_l2_merged_skip.pdf")
xmax =  (max(tempfrac0))
plot(-10,-10, xlim=c(0,xmax), ylim=c(0,nn+ni), xlab="Fraction", ylab="Density" ,axes=F)
axis(1)
tt = lapply(oi, function(i) abline( h = ((i-1)*ni)  ,  col=makeTransparent(1), lty=2 ) ) 
tt = lapply(oi, function(i) polygon(tempfrac0lf.d2[[i]][,1],((i-1)*ni)+(tempfrac0lf.d2[[i]][,2]*nni/n) ,  border=NA, col= (colpals2a[i,10] ) ) ) 
tt = lapply(oi, function(i) polygon(tempfrac0lm.d2[[i]][,1],((i-1)*ni)-(tempfrac0lm.d2[[i]][,2]*nni/n) ,  border=NA, col= (colpals2a[i,10] ) ) ) 
tt = lapply(oi, function(i) segments( mean(tempfrac0lf[i,]),  ((i-1)*ni) , 
                                       mean(tempfrac0lf[i,]), ((i-1)*ni) + ni/2 , lwd=2, col= colsex[2]  ) )

tt = lapply(oi, function(i) segments( mean(tempfrac0lm[i,]),  ((i-1)*ni) , 
                                       mean(tempfrac0lm[i,]), ((i-1)*ni) - ni/2 , lwd=2, col= colsex[1]) )


text(xmax, ((oi -1)*ni), xxx)
dev.off() 
```

### Panel E FDRs vesus ratio   
```{r}
m = match( props0$BaselineProp.clusters, colpals2[,9] ) 
sigs = props0$FDR < 0.05

pdf("props_fdrs_ratios.pdf")
plot(props0$PropRatio, -log10(props0$FDR) , pch=19, xlab="Male:Female", ylab="-log10(FDR)", col=colpals2[m,10], cex= (1+sigs)  )
abline(v=1, col="grey", lwd=3, lty=2) 
abline(h=-log10(0.05 ), col=3  )
text(props0$PropRatio[sigs], -log10(props0$FDR)[sigs] + 0.2, props0$BaselineProp.clusters[sigs])
dev.off() 
```

#### Supp figg correlations with age 
```{r}
ff = phens$sex==2
fm = phens$sex==1
 phens$age =  as.numeric(phens$age)
nC = c(5,6)
n = dim(tempfrac0)[1] 

cors = sapply(1:n, function(i) cor.test( phens$age, tempfrac0[i,]  , m= "s")  ) 
cors_e = unlist(cors[4,])
pvals = unlist(cors[3,] ) 
cors_all = cbind(rownames(tempfrac0) , pvals, p.adjust(pvals), cors_e)
collabels = cbind(rownames(tempfrac0) , colpals2a[,10]) 


corsf = sapply(1:n, function(i) cor.test( phens$age[ff], tempfrac0[i,ff]  , m= "s")  ) 
corsf_e = unlist(corsf[4,])
pvalsf = unlist(corsf[3,] ) 

corsm = sapply(1:n, function(i) cor.test( phens$age[fm], tempfrac0[i,fm]  , m= "s")  ) 
corsm_e = unlist(corsm[4,])
pvalsm = unlist(corsm[3,] ) 

cors_all = cbind( cors_all, pvalsf, p.adjust(pvalsf), corsf_e, pvalsm, p.adjust(pvalsm), corsm_e)


temp = apply(cors_all[,c(4,7,10)],2, as.numeric ) 
temp2 = apply(cors_all[,c(4,7,10)-1],2, as.numeric ) 
rownames(temp2) = rownames(tempfrac0) 
rownames(temp) = rownames(tempfrac0) 
rownames(cors_all) = rownames(tempfrac0)  

save(cors_all, file="cors_counts_age.Rdata") 

a =  (c(temp[,1], temp[,2], temp[,3]) )
b =  (c(rep("All",n), rep("Female",n), rep("Male",n))) 
d = rep(colpals2a[,10],3) 
 

aa =  -log10(c(temp2[,1], temp2[,2], temp2[,3]) )  
bb = 2 *(aa >= -log10(0.05) ) + 1 
pdf("correlations_age_l2_merged_skip.pdf")
beeswarm(a~b, pch=19, pwcol=d , ylab="Correlations (rho)", xlab="", bty="n",pwcex=bb) 
dev.off()
```
 


 


