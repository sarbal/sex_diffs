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


