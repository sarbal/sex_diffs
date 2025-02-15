# Figure 2: Distributions of cell-type proportions across sex. 

### Files needed 
```{r} 
in_seurat_file = "cell_type.RDS"
in_indiv_phen_file = "phen.RDS"
in_meta_file = "metadata.Rdata"
umap = "umap.Rdata"
umapXY = "umap_featuresXY.Rdata"
source("helper.r")
load("palette2.Rdata")


```
### Read in seurat object and individual metadata 
```{r}
obj   <- readRDS(in_seurat_file)
phens <- readRDS(in_indiv_phen_file)
load( in_meta_file) 
obj@meta.data = metadata 
load(umap) 
load(umapXY) 
```

### Cell counts 
```{r}
freq1 = plyr:: count( obj$predicted.celltype.l2[ obj$sex=="1"] )
freq2 = plyr:: count( obj$predicted.celltype.l2[ obj$sex=="2"] )
freq = plyr:: count( obj$predicted.celltype.l2 )
freq_sex_cell = cbind(freq, freq1[,2], freq2[,2]  )
colnames(freq_sex_cell) = c("celltype", "total", "male", "female" )
save(freq_sex_cell, file="freq_sex_cell_l2.Rdata")
```

### Cell proportions analysis 
```{r}
library(speckle)
load("phens.Rdata")
```

#### Props_tests_cell_l2_merged
```{r}
filt = metadata$predicted.celltype.l3!="Doublet" & !is.na(metadata$predicted.celltype.l3) 
freq1 = plyr::count( cbind( metadata$individual ,as.character(metadata$predicted.celltype.l3))[filt,] ) 
freq_sexp0 = spread(freq1, key = "x.2", value = "freq")
tenmp = t(freq_sexp0[ ,-1])
tenmp[is.na(tenmp)] = 0
tenmp = as.matrix(tenmp)
tempfrac0 = (sapply(1:dim(tenmp)[2], function(i) tenmp[,i]/sum(tenmp[,i], na.rm=T) ) )
colnames(tempfrac0) = freq_sexp0[,1]
m = match(colnames(tempfrac0), phens[,2] )
phens = phens[m,]

props0 =  propeller(clusters=metadata$predicted.celltype.l3[filt], sample=metadata$individual[filt], group=metadata$sex[filt])

m = match( rownames(tempfrac0),props0[,1]  )
props0 = props0[m,]
 
tempfrac0lf = tempfrac0[,phens$sex==2]
tempfrac0lm = tempfrac0[,phens$sex==1] 

bw=0.05
tempfrac0lf.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,phens$sex==2], bw=bw, from=0, to=1) )
tempfrac0lm.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,phens$sex==1], bw=bw, from=0, to=1) )
tempfrac0l.d  = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,], bw=bw, from=0, to=1) )

tempfrac0l.d2  = lapply(1:length(tempfrac0l.d),  function(i) cbind( c(0,tempfrac0l.d[[i]]$x,1),  c(0,tempfrac0l.d[[i]]$y,0) ) )
tempfrac0lf.d2 = lapply(1:length(tempfrac0lf.d), function(i) cbind( c(0,tempfrac0lf.d[[i]]$x,1), c(0,tempfrac0lf.d[[i]]$y,0) ) )
tempfrac0lm.d2 = lapply(1:length(tempfrac0lm.d), function(i) cbind( c(0,tempfrac0lm.d[[i]]$x,1), c(0,tempfrac0lm.d[[i]]$y,0) ) )
 

save(tempfrac0,freq_sexp0, props0, tempfrac0lf, tempfrac0lm, phens , tempfrac0lf.d, tempfrac0lf.d2,  tempfrac0lm.d, tempfrac0lm.d2, tempfrac0l.d, tempfrac0l.d2, file="props_tests_cell_l2_merged.Rdata")
```



#### Props_tests_cell_l2_merged_skip
```{r}
skip = c("Doublet",  "HSPC", "Platelet", "Eryth")
filt = is.na(match( metadata$predicted.celltype.l3 , skip)  ) & !is.na(metadata$predicted.celltype.l3)  
freq1 = plyr::count( cbind( metadata$individual ,as.character(metadata$predicted.celltype.l3))[filt,] ) 
freq_sexp0 = spread(freq1, key = "x.2", value = "freq")
tenmp = t(freq_sexp0[ ,-1])
tenmp[is.na(tenmp)] = 0
tenmp = as.matrix(tenmp)
tempfrac0 = (sapply(1:dim(tenmp)[2], function(i) tenmp[,i]/sum(tenmp[,i], na.rm=T) ) )
colnames(tempfrac0) = freq_sexp0[,1]
m = match(colnames(tempfrac0), phens[,2] )
phens = phens[m,]

props0 =  propeller(clusters=metadata$predicted.celltype.l3[filt], sample=metadata$individual[filt], group=metadata$sex[filt])

m = match( rownames(tempfrac0),props0[,1]  )
props0 = props0[m,]
 
tempfrac0lf = tempfrac0[,phens$sex==2]
tempfrac0lm = tempfrac0[,phens$sex==1] 

bw=0.05
tempfrac0lf.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,phens$sex==2], bw=bw, from=0, to=1) )
tempfrac0lm.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,phens$sex==1], bw=bw, from=0, to=1) )
tempfrac0l.d  = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,], bw=bw, from=0, to=1) )

tempfrac0l.d2  = lapply(1:length(tempfrac0l.d),  function(i) cbind( c(0,tempfrac0l.d[[i]]$x,1),  c(0,tempfrac0l.d[[i]]$y,0) ) )
tempfrac0lf.d2 = lapply(1:length(tempfrac0lf.d), function(i) cbind( c(0,tempfrac0lf.d[[i]]$x,1), c(0,tempfrac0lf.d[[i]]$y,0) ) )
tempfrac0lm.d2 = lapply(1:length(tempfrac0lm.d), function(i) cbind( c(0,tempfrac0lm.d[[i]]$x,1), c(0,tempfrac0lm.d[[i]]$y,0) ) )
 

save(tempfrac0,freq_sexp0, props0, tempfrac0lf, tempfrac0lm, phens,  tempfrac0lf.d, tempfrac0lf.d2,  tempfrac0lm.d, tempfrac0lm.d2, tempfrac0l.d, tempfrac0l.d2, file="props_tests_cell_l2_merged_skip.Rdata")
```


#### Props_tests_cell_l2_merged_skip2
```{r}
skip = c("Doublet",  "HSPC", "Platelet", "Eryth", "ILC", "CD4 Proliferating", "CD8 Proliferating")
 filt = is.na(match( metadata$predicted.celltype.l3 , skip)  ) & !is.na(metadata$predicted.celltype.l3)  
freq1 = plyr::count( cbind( metadata$individual ,as.character(metadata$predicted.celltype.l3))[filt,] ) 
freq_sexp0 = spread(freq1, key = "x.2", value = "freq")
tenmp = t(freq_sexp0[ ,-1])
tenmp[is.na(tenmp)] = 0
tenmp = as.matrix(tenmp)
tempfrac0 = (sapply(1:dim(tenmp)[2], function(i) tenmp[,i]/sum(tenmp[,i], na.rm=T) ) )
colnames(tempfrac0) = freq_sexp0[,1]
m = match(colnames(tempfrac0), phens[,2] )
phens = phens[m,]

props0 =  propeller(clusters=metadata$predicted.celltype.l3[filt], sample=metadata$individual[filt], group=metadata$sex[filt])

m = match( rownames(tempfrac0),props0[,1]  )
props0 = props0[m,]
 
tempfrac0lf = tempfrac0[,phens$sex==2]
tempfrac0lm = tempfrac0[,phens$sex==1] 

bw=0.05
tempfrac0lf.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,phens$sex==2], bw=bw, from=0, to=1) )
tempfrac0lm.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,phens$sex==1], bw=bw, from=0, to=1) )
tempfrac0l.d  = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,], bw=bw, from=0, to=1) )

tempfrac0l.d2  = lapply(1:length(tempfrac0l.d),  function(i) cbind( c(0,tempfrac0l.d[[i]]$x,1),  c(0,tempfrac0l.d[[i]]$y,0) ) )
tempfrac0lf.d2 = lapply(1:length(tempfrac0lf.d), function(i) cbind( c(0,tempfrac0lf.d[[i]]$x,1), c(0,tempfrac0lf.d[[i]]$y,0) ) )
tempfrac0lm.d2 = lapply(1:length(tempfrac0lm.d), function(i) cbind( c(0,tempfrac0lm.d[[i]]$x,1), c(0,tempfrac0lm.d[[i]]$y,0) ) )
 

save(tempfrac0,freq_sexp0, props0, tempfrac0lf, tempfrac0lm, phens,  tempfrac0lf.d, tempfrac0lf.d2,  tempfrac0lm.d, tempfrac0lm.d2, tempfrac0l.d, tempfrac0l.d2, file="props_tests_cell_l2_merged_skip2.Rdata")
```

#### Props_tests_cell_l2_original_skip
```{r}
skip = c("Doublet",  "HSPC", "Platelet", "Eryth")
filt = is.na(match( metadata$predicted.celltype.l2 , skip)  ) & !is.na(metadata$predicted.celltype.l2)  
freq1 = plyr::count( cbind( metadata$individual ,as.character(metadata$predicted.celltype.l2))[filt,] ) 
freq_sexp0 = spread(freq1, key = "x.2", value = "freq")
tenmp = t(freq_sexp0[ ,-1])
tenmp[is.na(tenmp)] = 0
tenmp = as.matrix(tenmp)
tempfrac0 = (sapply(1:dim(tenmp)[2], function(i) tenmp[,i]/sum(tenmp[,i], na.rm=T) ) )
colnames(tempfrac0) = freq_sexp0[,1]
m = match(colnames(tempfrac0), phens[,2] )
phens = phens[m,]

props0 =  propeller(clusters=metadata$predicted.celltype.l2[filt], sample=metadata$individual[filt], group=metadata$sex[filt])

m = match( rownames(tempfrac0),props0[,1]  )
props0 = props0[m,]
 
tempfrac0lf = tempfrac0[,phens$sex==2]
tempfrac0lm = tempfrac0[,phens$sex==1] 

bw=0.05
tempfrac0lf.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,phens$sex==2], bw=bw, from=0, to=1) )
tempfrac0lm.d = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,phens$sex==1], bw=bw, from=0, to=1) )
tempfrac0l.d  = lapply(1:dim(tempfrac0)[1], function(i) density(tempfrac0[i,], bw=bw, from=0, to=1) )

tempfrac0l.d2  = lapply(1:length(tempfrac0l.d),  function(i) cbind( c(0,tempfrac0l.d[[i]]$x,1),  c(0,tempfrac0l.d[[i]]$y,0) ) )
tempfrac0lf.d2 = lapply(1:length(tempfrac0lf.d), function(i) cbind( c(0,tempfrac0lf.d[[i]]$x,1), c(0,tempfrac0lf.d[[i]]$y,0) ) )
tempfrac0lm.d2 = lapply(1:length(tempfrac0lm.d), function(i) cbind( c(0,tempfrac0lm.d[[i]]$x,1), c(0,tempfrac0lm.d[[i]]$y,0) ) )
 

save(tempfrac0,freq_sexp0, props0, tempfrac0lf, tempfrac0lm, phens,  tempfrac0lf.d, tempfrac0lf.d2,  tempfrac0lm.d, tempfrac0lm.d2, tempfrac0l.d, tempfrac0l.d2, file="props_tests_cell_l2_original_skip.Rdata")
```


 
