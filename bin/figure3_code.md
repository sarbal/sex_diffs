# Figure 3: Sex-differential expression

### Load in data
```{r}
source("helper.r")
# library(tidyverse)
exprs_dir = "cellranger_outs/"
out_dir = "markerDE/"
load("metadata.Rdata") 
setwd(exprs_dir)


load("annotations/annots_x.Rdata")
extra = read.table("sex_biased_genes_immune_schemidel", header=T, sep="\t")
data = extra[-c(617,781),1:3]
data_mat = as.matrix(pivot_wider(data, names_from="Cell.type", values_from = "Log2.fold.change.F.vs..M") )
genes_mat = data_mat[,1]
data_mat = data_mat[,-1]
data_mat = apply(data_mat, 2, as.numeric)
data_x = 1*(data_mat > 0)
data_x2 = -1*(data_mat < 0)
data_x = data_x2 + data_x
data_x[is.na(data_x) ] = 0 
rownames(data_x) = genes_mat
save(data_x, data_mat, genes_mat, file="annotations/annots_extra_sch.Rdata")


```
### Calculate DE by celltype*
Note, "l3" is not the current l3 from Azimuth but a modified (ie merged cells) version of l2. 

```{r}
freq = plyr::count( metadata$predicted.celltype.l3 ) 
celltypes = freq[,1]
celltypes = celltypes[-c(17,29)]

 for(celltype in celltypes){ 
    filt = !is.na(obj@meta.data$predicted.celltype.l3) &  obj@meta.data$predicted.celltype.l3 == celltype & !is.na(obj@meta.data$sex) 

    temp = obj[,filt]
    test = FindMarkers(temp, ident.1="1", ident.2="2", logfc.threshold=0, test.use="wilcox", group.by="sex")
    save(test, file = paste0("wilcox_DE_", celltype, ".Rdata"))
 
   test = FindMarkers(temp, ident.1="1", ident.2="2", logfc.threshold=0, test.use="bimod", group.by="sex")
   save(test, file = paste0("bimod_DE_", celltype, ".Rdata"))

   test = FindMarkers(temp, ident.1="1", ident.2="2", logfc.threshold=0, test.use="roc", group.by="sex")
   save(test, file = paste0("roc_DE_", celltype, ".Rdata"))

  test = FindMarkers(temp, ident.1="1", ident.2="2", logfc.threshold=0, test.use="MAST", group.by="sex")
  save(test, file = paste0("mast_DE_", celltype, ".Rdata"))

} 
```

### Parse DEGs 
```{r}
sex_genes_mat_wilcox = sex_genes_mat[,grep("wilcox", colnames(sex_genes_mat)  ) ]
sex_genes_mat_wilcox_up = sex_genes_mat_wilcox[,grep("up", colnames(sex_genes_mat_wilcox)  ) ]
sex_genes_mat_wilcox_down = sex_genes_mat_wilcox[,grep("down", colnames(sex_genes_mat_wilcox)  ) ]
```

#### DEG enrichment analysis  
```{r}
enrich_degs_x = lapply(1:dim(sex_genes_mat_wilcox)[2], function(i) gene_set_enrichment( names(which(sex_genes_mat_wilcox[,i] > 0) ) , annots_x, voc) )
pvals = sapply(1:length(enrich_degs_x), function(i) enrich_degs_x[[i]][,5])

enrich_degs_x_exp = lapply(1:dim(sex_genes_mat_wilcox)[2], function(i) gene_set_enrichment( names(which(sex_genes_mat_wilcox[,i] > 0) ) , annots_x_exp, voc_exp) )
pvals_exp = sapply(1:length(enrich_degs_x_exp), function(i) enrich_degs_x_exp[[i]][,5])

pvals_exp = apply(pvals_exp, 2, as.numeric)
rownames(pvals_exp) =  enrich_degs_x_exp[[1]][,1]
colnames(pvals_exp) =  colnames(sex_genes_mat_wilcox) 

padj = apply(pvals_exp, 2, p.adjust)
rownames(padj) =  rownames(pvals_exp) 
colnames(padj) =  colnames(pvals_exp) 

pp = sapply(1:length(enrich_degs_x_exp), function(i) enrich_degs_x_exp[[i]][,3])
pp = apply(pp, 2, as.numeric)
rownames(pp) =  rownames(pvals_exp) 
colnames(pp) =  colnames(pvals_exp) 

pp = pp + 1  


skip = c("Eryth", "HSPC", "Platelet")
filt =  is.na(match(  celltypes, skip) ) &  grepl("up|down", colnames(padj)  ) 

filt2 =  (rowSums(pp)  > 0 ) 
tp = rowSums(pp)
tp[ grep("Xi", rownames(pp))  ]  = 0 
tp[14] = 1 # escape genes
tp[34] = 1 # whlbld specific escape
tp[51:80] = 0 
filt3  = tp > 0


row_cols = (rbind(gene_set_cols, rep("white", length(gene_set_cols))))
col_cols = (rbind(colcols, colsex2))
col_filt = which(filt)
row_filt = which(filt3)

pdf("sex_degs_enrich_dotplot.pdf", width=15, height=10)
dotplot_extra_dendro( padj, pvals_exp, pp , rownames(padj), celltypes,  filt_row = row_filt , filt_col =col_filt   , colcol= col_cols , rowcol= row_cols  , rowden = T , colden = T ) 
dev.off() 


save(pvals_exp, padj, row_filt, col_filt, col_cols, row_cols, pp, celltypes, enrich_degs_x_exp, file="data/enrich_degs_x_exp.Rdata")
```



### Identify marker genes
```{r]
ni = 1 
annotf.list = list() 
annotm.list = list() 
annot.list = list() 

for(ni in c(1:77)){ 
  print(ni)
   load(files[ni])
    barcodes = colnames(data)
    genes = rownames(data)
    pool = gsub("Sample", "pool_", sample )
    umis = substr( rownames(metadata)[which(metadata$pool == pool)], 0, 16)
    sex  = metadata$sex[which(metadata$pool == pool)]    
    indv =  metadata$individual[which(metadata$pool == pool)]
    cell_type =   metadata$predicted.celltype.l3[which(metadata$pool == pool)]
    m = match(umis, barcodes)
    fu = !is.na(m)
    fd = m[fu]
    data = data[,fd]
    cell_type = cell_type[fu]
		indv = indv[fu]
    sex = sex[fu]
    umis = umis[fu]
    meta = data.frame( sex, cell_type, indv)
    rownames(meta) = umis

    acells <- CreateSeuratObject(counts = data, project = "onek1k_coexp",  min.cells = 3, min.features = 200, meta.data = meta)
    acells[["percent.mt"]] <- PercentageFeatureSet(acells, pattern = "^MT-")
    acells <- NormalizeData(acells, normalization.method = "LogNormalize", scale.factor = 10000)
    all.genes <- rownames(acells)
    acells <- ScaleData(acells, features = all.genes)
    
    acells@active.ident = factor(acells$cell_type) 
    filtf = acells$sex == 2 
    filtm = acells$sex == 1 
    
    all.markers_female <- FindAllMarkers(acells[,filtf], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
    all.markers_male <- FindAllMarkers(acells[,filtm], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )
    all.markers <- FindAllMarkers(acells[, ], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, )


    save(all.markers,all.markers_female,all.markers_male,  file=paste0(out_dir, pool, ".sex_markers.Rdata"))
     rm(acells)
   }

   annotf =  make_annotations( all.markers_female[,c(7,6)], unique(all.markers_female$gene),  unique(all.markers_female$cluster) )
   annotf.list[[ni]] = annotf
	 annotm =  make_annotations( all.markers_male[,c(7,6)], unique(all.markers_male$gene),  unique(all.markers_male$cluster) )
   annotm.list[[ni]] = annotm
   annot =  make_annotations( all.markers[,c(7,6)], unique(all.markers$gene),  unique(all.markers$cluster) )
   annot.list[[ni]] = annot
    
}

files = grep("sex", dir() , val=T)
annotf = list() 
annotm = list() 
annot = list() 

for(file in files){
	  temp =  gsub(".sex_markers.Rdata", "", file)
		load(file)
		annotf[[file]] = cbind(all.markers_female[all.markers_female$p_val_adj < 0.05,c(7,6)], temp)
		annotm[[file]] = cbind(all.markers_male[all.markers_male$p_val_adj < 0.05,c(7,6)], temp)
		annot[[file]] = cbind(all.markers[all.markers$p_val_adj < 0.05,c(7,6)], temp)

} 

f_mat = do.call(rbind,annotf)
m_mat = do.call(rbind,annotm)
a_mat = do.call(rbind,annot)


 test_a = tapply(a_mat[,1], a_mat[,2], plyr::count )
 test_f = tapply(f_mat[,1], f_mat[,2], plyr::count )
 test_m = tapply(m_mat[,1], m_mat[,2], plyr::count )

test_f2 = do.call(rbind, test_f)
test_m2 = do.call(rbind, test_m)
test_a2 = do.call(rbind, test_a)

 f2 = do.call(rbind, strsplit(rownames(test_f2), "\\." ) )[,1]
 m2 = do.call(rbind, strsplit(rownames(test_m2), "\\." ) )[,1]
 a2 = do.call(rbind, strsplit(rownames(test_a2), "\\." ) )[,1]

test_f3 = cbind(f2,test_f2) 
test_m3 = cbind(m2,test_m2) 
test_a3 = cbind(a2,test_a2) 

nnn = 60
freq_f = test_f3[test_f3[,3] >= nnn ,]
freq_m = test_m3[test_m3[,3] >= nnn ,]
freq_a = test_a3[test_a3[,3] >= nnn,]


annot_f = spread( freq_f, key=1, value=3, fill=0)
annot_m = spread( freq_m, key=1, value=3, fill=0)
annot_a = spread( freq_a, key=1, value=3, fill=0)
 
rownames(annot_f) = annot_f[,1]
rownames(annot_m) = annot_m[,1]
rownames(annot_a) = annot_a[,1]
 
annot_a2 = (annot_a[,-1]>0) * 1
annot_f2 = (annot_f[,-1]>0) * 1
annot_m2 = (annot_m[,-1]>0) * 1
 
nnn = 75
freq_f = test_f3[test_f3[,3] >= nnn ,]
freq_m = test_m3[test_m3[,3] >= nnn ,]
freq_a = test_a3[test_a3[,3] >= nnn,]
```



#### Marker enrichment analysis  
```{r}
load("annotations/annots_az.Rdata")
maz.voc =   cbind(colnames(annots_az),  colSums(annots_az), c(rep("findmarkers", 23), rep("azim", 31)) ) 

enrich.marker.c3 = lapply(1:dim(annots_az)[2], function(i) gene_set_enrichment(rownames(annots_az)[annots_az[,i]>0],  annotsc3,  c3.voc   ))
enrich.marker.c2 = lapply(1:dim(annots_az)[2], function(i) gene_set_enrichment(rownames(annots_az)[annots_az[,i]>0],  annotsc2,  c2.voc   ))
enrich.marker.c7 = lapply(1:dim(annots_az)[2], function(i) gene_set_enrichment(rownames(annots_az)[annots_az[,i]>0],  annotsc7,  c7.voc   ))
enrich.marker.h = lapply(1:dim(annots_az)[2], function(i) gene_set_enrichment(rownames(annots_az)[annots_az[,i]>0],  annotsh,  h.voc   ))

load(file="immport.Rdata")
enrich.marker.immp = lapply(1:dim(annots_az)[2], function(i) gene_set_enrichment(rownames(annots_az)[annots_az[,i]>0],  annots,  voc   ))

load("annotations/annots_hgnc.Rdata")
enrich.marker.hgnc = lapply(1:dim(annots_az)[2], function(i) gene_set_enrichment(rownames(annots_az)[annots_az[,i]>0],  annots,  voc   ))

save(enrich.marker.immp,enrich.marker.hgnc, enrich.marker.c3,  enrich.marker.c2, enrich.marker.c7 ,  enrich.marker.h,  maz.voc,  file="enrich.markers.Rdata")
```

#### MASH analysis 
```{r}
library(mashr)
library(ashr)
# Read in data results 
all_genes = list()
res_all = list() 
for(celltype in celltypes){ 
  load(paste0("sexMAST/mast_DE_v5A_", celltype ,".Rdata")) 
  mast = test 
  load(paste0("sexMAST/mast_DE_v5B_", celltype ,".Rdata")) 
  mast_adj = test 
  load(paste0("sexDE/wilcox_DE_", celltype, ".Rdata")) 
  orig = test
  all_genes = append(all_genes, rownames(orig))
  res_all[[celltype]] = list(mast, mast_adj, orig) 
}

all_genes = unique(unlist(all_genes))
orig_all = matrix(NA,nrow=length(all_genes), ncol = length(celltypes)) 
rownames(orig_all) = all_genes
colnames(orig_all) = celltypes
effsize_all = orig_all * 0 
for( celltype in celltypes){ 
  orig = res_all[[celltype]][[3]]
  m = match(all_genes, rownames(orig) )
  f1 = !is.na(m)
  f2 = m[f1]
  orig_all[f1,celltype] = orig$p_val[f2]  
  effsize_all[f1,celltype] = orig$avg_log2FC[f2]  
}
 effsize_all[is.na(effsize_all)]  = 0  # set missing to 0  
 orig_all[is.na(orig_all)] = 1 # set missing to 1 

# remove some cell types 
celltypes = colnames(zscore_data)
skip = c("Eryth", "HSPC", "Platelet")
filt =  is.na(match(  celltypes, skip) )
 
tests = lapply(1:length(orig_all[,1]), function(i) rank(orig_all[i,filt])* orig_all[i,]/rank(orig_all[i,filt]) ) 
simes_pvals = sapply(1:length(tests), function(i) min(tests[[i]]))
simes_padj = p.adjust(simes_pvals, method = "bonf") # can use the adjusted p-values in simes instead, or adjust post 
save(orig_all,simes_padj,effsize_all, file="data_for_mash.Rdata")

zscore_data = apply( orig_all , 2, scale) 
data = mash_set_data(zscore_data, (zscore_data*0)+1)  


# Canonical covariance 
U.c = cov_canonical(data) ## Set up the covariance matrices  
m.c = mash(data, U.c) ## Fit the model 
print(get_loglik(m),digits = 10)

# Data driven 
m.1by1 = mash_1by1(data) #  select strong signals 
strong = get_significant_results(m.1by1,0.05)
U.pca = cov_pca(data,5,subset=strong) # Obtain initial data-driven covariance matrices
U.ed = cov_ed(data, U.pca, subset=strong) # Apply Extreme Deconvolution
m.ed = mash(data, U.ed) ## Fit the model 
print(get_loglik(m.ed),digits = 10)

# Data driven and canonical analysis 
U.c = cov_canonical(data)  
m   = mash(data, c(U.c,U.ed))
print(get_loglik(m),digits = 10)

## Accounting for correlations among measurements - simple 
V.simple = estimate_null_correlation_simple(data) # estimate correlations:
data.Vsimple = mash_update_data(data, V=V.simple)
U.c = cov_canonical(data.Vsimple) 
m.Vsimple = mash(data.Vsimple, U.c) # fits with correlations because data.V includes correlation information 
print(get_loglik(m.Vsimple),digits=10) # log-likelihood of the fit with correlations set to V

m.orig = mash(data, U.c) # fits without correlations because data object was set up without correlations
print(get_loglik(m.orig),digits=10)


## Accounting for correlations among measurements - Estimate the residual correlations using EM method
V.em = mash_estimate_corr_em(data, U.c, details = TRUE)
m.Vem = V.em$mash.model
print(get_loglik(m.Vem),digits=10) # log-likelihood of the fit


loglik = c(get_loglik(m.orig), get_loglik(m.Vsimple), get_loglik(m.Vem), 
           get_loglik(m), get_loglik(m.ed), get_loglik(m.c) ) 

significant = c( length(get_significant_results(m.orig)), length(get_significant_results(m.Vsimple)),
                  length(get_significant_results(m.Vem)), length(get_significant_results(m)),
                  length(get_significant_results(m.ed)), length(get_significant_results(m.c)) )


tb = rbind(loglik, significant)
colnames(tb) = c('without cor', 'V simple', 'V EM', 'joint', 'data driven', 'canonical' )
row.names(tb) = c('log likelihood', '# significance' )
 
x = get_pairwise_sharing(m.Vem, factor=0.5)
corrplot(x, method='color', col.lim=c(0,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45, title = 'Pairwise Sharing by Magnitude', mar = c(4,0,4,0), number.cex = 0.5)


## filter/compare with original results 
mm = match( names(get_significant_results(m.Vem)),rownames(sex_genes_mat_wilcox) )
f.v= !is.na(mm)
f.oo = mm[f.v]
recur_down = rowSums(wilcox_gene_mat_down[f.oo,filt])
recur_up = rowSums(wilcox_gene_mat_up[f.oo,filt])

 

```

 


