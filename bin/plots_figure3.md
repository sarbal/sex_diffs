# Figure 3: Sex-differential expression.

### Panel A - umap features XY 
```{r}
load("umap_featuresXY.Rdata") 
obj@reductions$umap  = umap_v2 
 

filt =  !is.na(obj$predicted.celltype.l3) & !is.na(obj$sex)
temp = obj[,filt]



png("umap_xy.png", height=2000, width=2000)
DimPlot(temp, 
    reduction = "umap", 
    pt.size = 2, 
    raster=FALSE, 
    group.by = "sex" , 
    cols= colsex) + NoLegend()
dev.off()
 

 
m = match(levels(factor(obj$predicted.celltype.l3)), colpals5[,11] )
```
 
### Panel B - umap features XY 
```{r}
png("umap_l3xy.png", height=2000, width=2000)
DimPlot(obj, 
    reduction = "umap", 
    pt.size = 2, 
    raster=FALSE, 
    group.by = "predicted.celltype.l3" , 
    cols= colpals5[m,12]) + NoLegend() 
dev.off()
``` 
### Panel C - DEGs 
```{r} 
load("seurat_sex_degs.Rdata")
sex_genes_mat_mast = sex_genes_mat[,grep("mast-", colnames(sex_genes_mat)  ) ]
sex_genes_mat_mast_up = sex_genes_mat_mast[,grep("up", colnames(sex_genes_mat_mast)  ) ]
sex_genes_mat_mast_down = sex_genes_mat_mast[,grep("down", colnames(sex_genes_mat_mast)  ) ]

nsums = cbind( -colSums(sex_genes_mat_mast_up), colSums(sex_genes_mat_mast_down) )
colnames(nsums) = c("up", "down")
rownames(nsums) = gsub("mast-up ", "", rownames(nsums)  ) 
m = match(  colpals2[,9], rownames(nsums))
m = m[!is.na(m) ] 
nsums = nsums[m,] 
skip = c("Eryth", "HSPC", "Platelet")
filt =  is.na(match( rownames(nsums), skip) )
celltypes = rownames(nsums) 
pdf("sex_degs_celltype.pdf")
bp = barplot( t( nsums[filt,]) ,beside=T, space=c(-1,0.5), col=colsex, border=NA, xlab="Cell-type", ylab="Number of DE genes", hori=T)
text(20, bp[1,], celltypes[filt],  adj=0) 
dev.off() 




#
sex_genes_mat_wilcox = sex_genes_mat[,grep("wilcox", colnames(sex_genes_mat)  ) ]
sex_genes_mat_wilcox_up = sex_genes_mat_wilcox[,grep("up", colnames(sex_genes_mat_wilcox)  ) ]
sex_genes_mat_wilcox_down = sex_genes_mat_wilcox[,grep("down", colnames(sex_genes_mat_wilcox)  ) ]

nsums = cbind( -colSums(sex_genes_mat_wilcox_up), colSums(sex_genes_mat_wilcox_down) )
colnames(nsums) = c("up", "down")
rownames(nsums) = gsub("wilcox-up ", "", rownames(nsums)  ) 
m = match(  colpals2[,9], rownames(nsums))
m = m[!is.na(m) ] 
nsums = nsums[m,] 
skip = c("Eryth", "HSPC", "Platelet")
filt =  is.na(match( rownames(nsums), skip) )
celltypes = rownames(nsums) 
pdf("../sex_degs_celltype_wilcox.pdf")
bp = barplot( t( nsums[filt,]) ,beside=T, space=c(-1,0.5), col=colsex[-3], border=NA, xlab="Cell-type", ylab="Number of DE genes", hori=T)
text(20, bp[1,], celltypes[filt],  adj=0) 
dev.off() 

#
load("mast_sex_degs.Rdata")
sex_genes_mat_mast_up = sex_deg_mast_sc[,grep("down", colnames(sex_deg_mast_sc)  ) ] # down is male biased here, flipping to keep consistent with other analyses 
sex_genes_mat_mast_down = sex_deg_mast_sc[,grep("up", colnames(sex_deg_mast_sc)  ) ]

nsums = cbind( -colSums(sex_genes_mat_mast_up), colSums(sex_genes_mat_mast_down) ) 
colnames(nsums) = c("up", "down")
rownames(nsums) = gsub("down - ", "", rownames(nsums)  ) 
m = match(  colpals2[,9], rownames(nsums))
m = m[!is.na(m) ] 
nsums = nsums[m,] 
skip = c("Eryth", "HSPC", "Platelet")
filt =  is.na(match( rownames(nsums), skip) )
celltypes = rownames(nsums) 
pdf("../../figures_2022_05/sex_degs_mast_celltype.pdf")
bp = barplot( t( nsums[filt,]) ,beside=T, space=c(-1,0.5), col=colsex[-3], border=NA, xlab="Cell-type", ylab="Number of DE genes", hori=T)
text(5, bp[1,], celltypes[filt],  adj=0) 
dev.off() 

sex_genes_mat_mast_up = sex_genes_mat_mast_up[,m]
sex_genes_mat_mast_down  = sex_genes_mat_mast_down[,m]
```



### Panel D - recurrence 
```{r}
recur_up = rowSums(sex_genes_mat_mast_up[,filt]) 
recur_down = rowSums(sex_genes_mat_mast_down[,filt]) 

recur_up = recur_up[recur_up>0] 
recur_down = recur_down[recur_down>0] 


hist(recur_down, breaks=c(0:24)) 
hist(recur_down, breaks=c(0:24)) 

fu = plyr::count(recur_up ) 
fd = plyr::count(recur_down ) 
x = 1:24
recur_temp = cbind(x*0, x*0)  
recur_temp[match(fu[,1], x) ,1] = -fu[,2] 
recur_temp[match(fd[,1], x) ,2] = fd[,2] 
rownames(recur_temp) = x 



pdf("sex_degs_recurrence.pdf")
bp = barplot( t(recur_temp), beside=T, space=c(-1,0), col=colsex[-3], border=NA, ylab="Number of genes", xlab="Recurrence") 
dev.off() 


## wilcox version 
recur_up = rowSums(sex_genes_mat_wilcox_up[,filt]) 
recur_down = rowSums(sex_genes_mat_wilcox_down[,filt]) 

recur_up = recur_up[recur_up>0] 
recur_down = recur_down[recur_down>0] 


hist(recur_down, breaks=c(0:24)) 
hist(recur_down, breaks=c(0:24)) 

fu = plyr::count(recur_up ) 
fd = plyr::count(recur_down ) 
x = 1:24
recur_temp = cbind(x*0, x*0)  
recur_temp[match(fu[,1], x) ,1] = -fu[,2] 
recur_temp[match(fd[,1], x) ,2] = fd[,2] 
rownames(recur_temp) = x 



pdf("sex_degs_recurrence_wilcox.pdf")
bp = barplot( t(recur_temp), beside=T, space=c(-1,0), col=colsex[-3], border=NA, ylab="Number of genes", xlab="Recurrence") 
dev.off() 

## 

recur_up = rowSums(sex_genes_mat_mast_up[,filt]) 
recur_down = rowSums(sex_genes_mat_mast_down[,filt]) 

recur_up = recur_up[recur_up>0] 
recur_down = recur_down[recur_down>0] 


hist(recur_down, breaks=c(0:24)) 
hist(recur_down, breaks=c(0:24)) 

fu = plyr::count(recur_up ) 
fd = plyr::count(recur_down ) 
x = 1:24
recur_temp = cbind(x*0, x*0)  
recur_temp[match(fu[,1], x) ,1] = -fu[,2] 
recur_temp[match(fd[,1], x) ,2] = fd[,2] 
rownames(recur_temp) = x 



pdf("sex_degs_mast_recurrence.pdf")
bp = barplot( t(recur_temp), beside=T, space=c(-1,0), col=colsex[-3], border=NA, ylab="Number of genes", xlab="Recurrence") 
dev.off() 

mast_down = lapply(1:length(temp), function(i) names(which(sex_genes_mat_mast_down[,i] > 0 ) )   ) 
mast_up = lapply(1:length(temp), function(i) names(which(sex_genes_mat_mast_up[,i] > 0 ) )   ) 
wilcox_up = lapply(1:length(temp2), function(i) names(which(sex_genes_mat_wilcox_up[,i] > 0 ) )   ) 
wilcox_down = lapply(1:length(temp2), function(i) names(which(sex_genes_mat_wilcox_down[,i] > 0 ) )   ) 

temp = gsub("up - ", "", colnames(sex_genes_mat_mast_down)  ) 
temp2 = gsub("wilcox-up ", "", colnames(sex_genes_mat_wilcox_up) ) 
temp3 = temp2[is.na(match( temp2, skip) ) ]

names(wilcox_down)  = temp2
names(wilcox_up)  = temp2
names(mast_up)  = temp
names(mast_down)  = temp

overlap_down = lapply(1:length(temp3), function(i) intersect( unlist(wilcox_down[temp3[i]]), unlist(mast_down[temp3[i]]) )  ) 
overlap_up = lapply(1:length(temp3), function(i) intersect( unlist(wilcox_up[temp3[i]]), unlist(mast_up[temp3[i]]) )  ) 

unique_down_wilcox = lapply(1:length(temp3), function(i) setdiff( unlist(wilcox_down[temp3[i]]), unlist(mast_down[temp3[i]]) )  ) 
unique_up_wilcox = lapply(1:length(temp3), function(i) setdiff( unlist(wilcox_up[temp3[i]]), unlist(mast_up[temp3[i]]) )  ) 

unique_down_mast = lapply(1:length(temp3), function(i) setdiff( unlist(mast_down[temp3[i]]), unlist(wilcox_down[temp3[i]]) )  ) 
unique_up_mast = lapply(1:length(temp3), function(i) setdiff( unlist(mast_up[temp3[i]]), unlist(wilcox_up[temp3[i]]) )  ) 


names(overlap_down)  = temp3
names(overlap_up)  = temp3
names(unique_down_wilcox)  = temp3
names(unique_up_wilcox)  = temp3
names(unique_down_mast)  = temp3
names(unique_up_mast)  = temp3
```




### Panel E - enrichment heatmap
```{r}
genes = rownames(sex_genes_mat) 
sex_genes_mat_mast = sex_genes_mat[,grep("mast", colnames(sex_genes_mat) ) ] 


load("annotations/annots_x.Rdata")
enrich_degs = lapply(1:dim(sex_genes_mat_mast)[2], function(i) gene_set_enrichment( genes[sex_genes_mat_mast[,i]>0] , annots_x, voc) )
enrich_degs = lapply(1:dim(sex_genes_mat_wilcox)[2], function(i) gene_set_enrichment( genes[sex_genes_mat_wilcox[,i]>0] , annots_x, voc) )



sex_genes_mat_mast = sex_genes_mat[,grep("mast", colnames(sex_genes_mat) ) ] 
sex_degs_mast = lapply(  grep("up|down", colnames(sex_genes_mat_mast), inver=T) , function(i) rownames(sex_genes_mat_mast)[sex_genes_mat_mast[,i] > 0]  )  
names(sex_degs_mast) =  gsub("mast ", "", colnames(sex_genes_mat_mast)) [grep("up|down", colnames(sex_genes_mat_mast), inver=T)] 
celltypes = names(sex_degs_mast)

 

pvals = sapply(1:length(enrich_degs), function(i) enrich_degs[[i]][,5])
pvals = apply(pvals, 2, as.numeric)
rownames(pvals) =  enrich_degs[[1]][,1]
#colnames(pvals) =  colnames(sex_genes_mat_mast) 
colnames(pvals) =  colnames(sex_genes_mat_wilcox) 


padj = apply(pvals, 2, p.adjust)
rownames(padj) =  rownames(pvals) 
colnames(padj) =  colnames(pvals) 

pp = sapply(1:length(enrich_degs), function(i) enrich_degs[[i]][,3])
pp = apply(pp, 2, as.numeric)
rownames(pp) =  rownames(pvals) 
colnames(pp) =  colnames(pvals) 

 
celltypes =  gsub("mast-up |mast-down |mast ", "", colnames(padj)) 
skip = c("Eryth", "HSPC", "Platelet")
filt =  is.na(match( celltypes , skip) )   

m = match( celltypes   , colpals2[,9]  ) 
f.p = !is.na(m)
f.c = m[f.p]
colcols = colpals2[f.c,10]

pdf("sex_degs_enrich_heatmap.pdf")
heatmap.3(-log10(padj[,filt] ), col=cols6, cexRow=0.5, ColSideCol=colcols[filt]  )
dev.off() 
 

celltypes = gsub("wilcox-up |wilcox-down |wilcox ", "", colnames(padj) )  
skip = c("Eryth", "HSPC", "Platelet")
filt =  is.na(match( celltypes, skip) )   

m = match( celltypes  , colpals2[,9]  ) 
f.p = !is.na(m)
f.c = m[f.p]
colcols = colpals2[f.c,10]

pdf("sex_degs_enrich_heatmap_wilcox.pdf")
heatmap.3(-log10(padj[,filt] ), col=cols6, cexRow=0.5, ColSideCol=colcols[filt]  )
dev.off() 

 
sex_genes_mat_mast = sex_genes_mat_wilcox
sex_genes_mat_mast_up = sex_genes_mat_wilcox_up
sex_genes_mat_mast_down = sex_genes_mat_wilcox_down


load(file="annotations/markers_azimuth.Rdata")
 
enrich_degs_az = lapply(1:dim(sex_genes_mat_mast)[2], function(i) gene_set_enrichment( genes[sex_genes_mat_mast[,i]>0] , markers_mat, voc2) )

load(file="data/markers_finddegs.Rdata")
enrich_degs2 = lapply(1:dim(sex_genes_mat_mast)[2], function(i) gene_set_enrichment( genes[sex_genes_mat_mast[,i]>0] , marker_mat_2, voc2) )
enrich_degs3 = lapply(1:dim(sex_genes_mat_mast)[2], function(i) gene_set_enrichment( genes[sex_genes_mat_mast[,i]>0] , marker_mat_3, voc3) )
 


pvalsaz = sapply(1:length(enrich_degs_az), function(i) enrich_degs_az[[i]][,5])
paz = sapply(1:length(enrich_degs_az), function(i) enrich_degs_az[[i]][,3])
pvalsaz = apply(pvalsaz, 2, as.numeric)
rownames(pvalsaz) =  enrich_degs_az[[1]][,1]
colnames(pvalsaz) =  colnames(sex_genes_mat_mast) 

padjaz = apply(pvalsaz, 2, p.adjust)
paz = apply(paz, 2, as.numeric)

rownames(padjaz) =  rownames(pvalsaz) 
colnames(padjaz) =  colnames(pvalsaz) 
rownames(paz) =  rownames(pvalsaz) 
colnames(paz) =  colnames(pvalsaz) 


pvals2 = sapply(1:length(enrich_degs2), function(i) enrich_degs2[[i]][,5])
p2 = sapply(1:length(enrich_degs2), function(i) enrich_degs2[[i]][,3])
pvals2 = apply(pvals2, 2, as.numeric)
p2 = apply(p2, 2, as.numeric)
rownames(pvals2) =  enrich_degs2[[1]][,1]
colnames(pvals2) =  colnames(sex_genes_mat_mast) 

padj2 = apply(pvals2, 2, p.adjust)
rownames(padj2) =  rownames(pvals2) 
colnames(padj2) =  colnames(pvals2) 
rownames(p2) =  rownames(pvals2) 
colnames(p2) =  colnames(pvals2) 

pvals3 = sapply(1:length(enrich_degs3), function(i) enrich_degs3[[i]][,5])
p3 = sapply(1:length(enrich_degs3), function(i) enrich_degs3[[i]][,3])
pvals3 = apply(pvals3, 2, as.numeric)
p3 = apply(p3, 2, as.numeric)

rownames(pvals3) =  enrich_degs3[[1]][,1]
colnames(pvals3) =  colnames(sex_genes_mat_mast) 

padj3 = apply(pvals3, 2, p.adjust)
rownames(padj3) =  rownames(pvals3) 
colnames(padj3) =  colnames(pvals3) 
rownames(p3) =  rownames(pvals3) 
colnames(p3) =  colnames(pvals3) 


m = match( colpals2[,9], rownames(padj2) ) 
f.c = !is.na(m)
f.p = m[f.c]

padj2 = padj2[f.p,]
pvals2 = pvals2[f.p,]
p2 = p2[f.p,]

m = match( colpals2[,9], rownames(padj3) ) 
f.c = !is.na(m)
f.p = m[f.c]

padj3 = padj3[f.p,]
pvals3 = pvals3[f.p,]
p3 = p3[f.p,]


m = match( colpals2[,9], rownames(padjaz) ) 
f.c = !is.na(m)
f.p = m[f.c]

padjaz = padjaz[f.p,]
pvalsaz = pvalsaz[f.p,]
paz = paz[f.p,]

paz = paz + 1 
p2 = p2 + 1 
p3 = p3 + 1 
pp = pp + 1 


skip = c("Eryth", "HSPC", "Platelet")
#filt =  is.na(match( celltypes , skip) )  & grepl("mast ", colnames(padj2) )   
filt =  is.na(match( celltypes , skip) )   & grepl("wilcox ", colnames(padj2) )   

m = match( colpals2[,9], celltypes[filt]    ) 
f.c = !is.na(m)
f.p = m[f.c]
colcols = colpals2[f.c,10]

padj_all = rbind( padj3[,filt][,f.p], NA , padjaz[,filt][,f.p])
p_all = rbind( p3[,filt][,f.p], NA , paz[,filt][,f.p])

filt2 =  is.na( match( rownames(padj_all)   , skip) )  
colnames(padj_all)  = gsub("wilcox-up ", "", colnames(padj_all) ) 
m = match(rownames(padj_all)[filt2], colpals2[,9])


heatmap.3(-log10(padj_all[filt2,] ), col=cols6, Rowv=F, Colv=F, cexRow=0.5, ColSideCol=colcols, RowSideCol = colpals2[m,10] )
heatmap.3( (p_all[filt2,] ), col=cols6, Rowv=F, Colv=F, cexRow=0.5, ColSideCol=colcols, RowSideCol = colpals2[m,10] )


pdf("sex_degs_enrich_markers_heatmap.pdf")
heatmap.3(-log10(padj_all[filt2,] ), col=cols6, Rowv=F, Colv=F, cexRow=0.5, ColSideCol=colcols, RowSideCol = colpals2[m,10] )
dev.off() 



 
  
dotplot <- function(padj, pval , pp, rowids, colids, filt_row, filt_col, rowcol="", colcol="") {
  
  rowids = rowids[filt_row]
  colids = colids[filt_col] 
   

  pval = pval[filt_row,filt_col]
  pp = pp[filt_row,filt_col]
  padj = padj[filt_row,filt_col]
 
  nc = length(colids)   
  nr = length(rowids)    
 
  dotsize = range(log(as.numeric(pp)+1), na.rm=T  )
  min_dot = dotsize[1]
  max_dot = dotsize[2]
  dots = seq( min_dot, max_dot, length.out=4 )[-1] 
  
   
  par(oma = c(1, 3, 2, 1))
  par(mar = c(2, 2, 2, 1))
  zones <- matrix(c(1, 3, 2, 4), ncol = 2, byrow = TRUE)
  layout(zones, widths = c(1,5,1,5), heights = c(1, 5, 1,1))
 

  nmax = range( ceiling(-log10(padj)    ) ) 
  temp = cbind( seq( nmax[1], nmax[2], length.out=101) , nmax[1] ,nmax[1] )
  colssig <- gplots::colorpanel(length(temp[,1]), "white", "red", "darkmagenta")

  # Plot legend 
  image(temp, axes=F, col=colssig,  xlab="", ylab="")
  text(0.5,0.5, "-log10 P-adj", font=2 ,  xpd=T )
  axis(1, at=(0:5/5), labels=  round(temp[100 * 0:5/5  + 1 ,1] , 1)  )

  # Plot rows   
  plot(rep(0,nr), 1:nr, axes=F, col=0, xlab="", ylab=" ",  ylim=c(0,nr) )
  text(x=1, y=1:nr,  labels =  rowids  ,  adj=1,  xpd=T)
  
  # Plot columns      
  plot(0, axes=F, col=0, main=" ", xlab="", ylab="", xlim=c(0,nc))
  text( x=1:nc , y=0, cex=1, labels = colids , srt=90, xpd=T, adj=1)

  # Plot legend
  points( rep(-nc*0.05,3), -(1:3)/2, cex=dots, pch=21, xpd=T )
  text( x=-nc*0.01 ,y=-(1:3)/2, cex=1, labels = round(exp(dots)) ,  xpd=T )


  # Plot box 
  plot(0,0, ylim=c(0,nr), xlim=c(0,nc), col=0, axes=F, xlab=" ", ylab="")
  segments ( y0 = 0:nr + 0.5,  x0=0.5 , x1 = nc+0.5, col="lightgrey")
  segments ( x0 = 0:nc + 0.5,  y0=0.5 , y1 = nr+0.5, col="lightgrey")
 
  
  temp = cbind( seq( nmax[1], nmax[2] +1  ) , nmax[1]/100 ,nmax[1]/100 )
  colssig <- gplots::colorpanel(length(temp[,1]), "white", "red", "darkmagenta")

  # Plot dots 
  for(i in 1:nr){ points( 1:nc, rep(i,nc),  cex= log(as.numeric(pp[i,])+1),  bg=  colssig[ match(round( -log10(padj[i,] ) + 1 )  , round(temp[,1]) ) ] ,    pch=21  ) }
 

  if( length(colcol) > 1  ){
     colcol = colcol[filt_col]
     points(  x=1:nc, y=rep(nr+1.5,nc), col = colcol,  pch=15, cex=3, xpd=T)    
  }
    
  if( length(rowcol) > 1  ){
      rowcol = rowcol[filt_row] 
      points(  y=1:nr, x=rep(-0.5,nr), col = rowcol,  pch=15, cex=3, xpd=T)
    }
  } 


 
dotplot_extra <- function(padj, pval , pp, rowids, colids, filt_row, filt_col, rowcol="", colcol="") {
  
  rowids = rowids[filt_row]
  colids = colids[filt_col] 
   

  pval = pval[filt_row,filt_col]
  pp = pp[filt_row,filt_col]
  padj = padj[filt_row,filt_col]
 
  nc = length(colids)   
  nr = length(rowids)    
 
  dotsize = range(log(as.numeric(pp)+1), na.rm=T  )
  min_dot = dotsize[1]
  max_dot = dotsize[2]
  dots = seq( min_dot, max_dot, length.out=4 )[-1] 
  
   
  par(oma = c(1, 3, 2, 1))
  par(mar = c(2, 2, 2, 1))
  zones <- matrix(c(1, 3, 2, 4), ncol = 2, byrow = TRUE)
  layout(zones, widths = c(1,5,1,5), heights = c(1, 5, 1,1))
 

  nmax = range( ceiling(-log10(padj)    ) ) 
  temp = cbind( seq( nmax[1], nmax[2], length.out=101) , nmax[1] ,nmax[1] )
  colssig <- gplots::colorpanel(length(temp[,1]), "white", "red", "darkmagenta")

  # Plot legend 
  image(temp, axes=F, col=colssig,  xlab="", ylab="")
  text(0.5,0.5, "-log10 P-adj", font=2 ,  xpd=T )
  axis(1, at=(0:5/5), labels=  round(temp[100 * 0:5/5  + 1 ,1] , 1)  )

  # Plot rows   
  plot(rep(0,nr), 1:nr, axes=F, col=0, xlab="", ylab=" ",  ylim=c(0,nr) )
  text(x=1, y=1:nr,  labels =  rowids  ,  adj=1,  xpd=T)
  
  # Plot columns      
  plot(0, axes=F, col=0, main=" ", xlab="", ylab="", xlim=c(0,nc))
  text( x=1:nc , y=0, cex=1, labels = colids , srt=90, xpd=T, adj=1)

  # Plot legend
  points( rep(-nc*0.05,3), -(1:3)/2, cex=dots, pch=21, xpd=T )
  text( x=-nc*0.01 ,y=-(1:3)/2, cex=1, labels = round(exp(dots)) ,  xpd=T )


  # Plot box 
  plot(0,0, ylim=c(0,nr), xlim=c(0,nc), col=0, axes=F, xlab=" ", ylab="")
  segments ( y0 = 0:nr + 0.5,  x0=0.5 , x1 = nc+0.5, col="lightgrey")
  segments ( x0 = 0:nc + 0.5,  y0=0.5 , y1 = nr+0.5, col="lightgrey")
 
  
  temp = cbind( seq( nmax[1], nmax[2] +1  ) , nmax[1]/100 ,nmax[1]/100 )
  colssig <- gplots::colorpanel(length(temp[,1]), "white", "red", "darkmagenta")

  # Plot dots 
  for(i in 1:nr){ points( 1:nc, rep(i,nc),  cex= log(as.numeric(pp[i,])+1),  bg=  colssig[ match(round( -log10(padj[i,]) + 1  )  , round(temp[,1]) ) ] ,    pch=21  ) }
 

  if( !is.null(dim(colcol))   ){
     colcol = colcol[,filt_col]
     nnr = dim(colcol)[1]
      for(i in 1:nnr){ 
          points(  x=1:nc, y=rep(nr+1.5+i,nc), col = colcol[i,],  pch=15, cex=3, xpd=T)    
      } 
  }
    
  if( !is.null(dim(rowcol))   ){
      rowcol = rowcol[,filt_row] 
      nnr = dim(rowcol)[1]
      for(i in 1:nnr){ 
         points(  y=1:nr, x=rep(-0.5-i,nr), col = rowcol[i,],  pch=15, cex=3, xpd=T)
       } 
    }
  } 





 
dotplot_extra_dendro <- function(padj, pval , pp, rowids, colids, filt_row, filt_col, rowcol="", colcol="", rowden=F, colden=F) {
  
  rowids = rowids[filt_row]
  colids = colids[filt_col] 
   

  pval = pval[filt_row,filt_col]
  pp = pp[filt_row,filt_col]
  padj = padj[filt_row,filt_col]
 
  nc = length(colids)   
  nr = length(rowids)    
 
  dotsize = range(log(as.numeric(pp)+1), na.rm=T  )
  min_dot = dotsize[1]
  max_dot = dotsize[2]
  dots = seq( min_dot, max_dot, length.out=4 )[-1] 
  
  if(rowden==T)  { 
    hc_row =  hclust( dist( padj) ) 
    padj = padj[hc_row$order,]
    pval = pval[hc_row$order,]
    pp   = pp[hc_row$order,]
    rowids= rowids[hc_row$order]
    filt_row = filt_row[hc_row$order]
  } 

  if(colden==T)  { 
    hc_col = hclust( dist( t(padj) )) 
    padj = padj[,hc_col$order]
    pval = pval[,hc_col$order]
    pp   = pp[,hc_col$order]
    colids= colids[hc_col$order]
    filt_col = filt_col[hc_col$order]
  } 


  par(oma = c(1, 3, 2, 1))
  par(mar = c(2, 2, 2, 1))
  zones <- matrix(c(1, 3, 2, 4), ncol = 2, byrow = TRUE)
  layout(zones, widths = c(1,5,1,5), heights = c(1, 5, 1,1))
 

  nmax = range( ceiling(-log10(padj)    ) ) 
  temp = cbind( seq( nmax[1], nmax[2], length.out=101) , nmax[1] ,nmax[1] )
  colssig <- gplots::colorpanel(length(temp[,1]), "white", "red", "darkmagenta")

  # Plot legend 
  image(temp, axes=F, col=colssig,  xlab="", ylab="")
  text(0.5,0.5, "-log10 P-adj", font=2 ,  xpd=T )
  axis(1, at=(0:5/5), labels=  round(temp[100 * 0:5/5  + 1 ,1] , 1)  )

  # Plot rows   
   if( rowden == T  ){
      plot( as.dendrogram(hc_row) , axes=F,  main=" ", xlab="", ylab="", ylim=c(0,nr), hori=T)
   } else { 
      plot(rep(0,nr), 1:nr, axes=F, col=0, xlab="", ylab=" ",  ylim=c(0,nr) )
      text(x=1, y=1:nr,  labels =  rowids  ,  adj=1,  xpd=T)
   }
  
  # Plot columns      
   if( colden == T  ){
     plot(as.dendrogram(hc_col) , axes=F,  main=" ", xlab="", ylab="", xlim=c(0,nc))

  } else { 
      plot(0, axes=F, col=0, main=" ", xlab="", ylab="", xlim=c(0,nc))
      text( x=1:nc , y=0, cex=1, labels = colids , srt=90, xpd=T, adj=1)
  }
    
  # Plot legend
  points( rep(-nc*0.05,3), -(1:3)/2, cex=dots, pch=21, xpd=T )
  text( x=-nc*0.01 ,y=-(1:3)/2, cex=1, labels = round(exp(dots)) ,  xpd=T )


  # Plot box 
  plot(0,0, ylim=c(0,nr), xlim=c(0,nc), col=0, axes=F, xlab=" ", ylab="")
  segments ( y0 = 0:nr + 0.5,  x0=0.5 , x1 = nc+0.5, col="lightgrey")
  segments ( x0 = 0:nc + 0.5,  y0=0.5 , y1 = nr+0.5, col="lightgrey")
 
  
  temp = cbind( seq( nmax[1], nmax[2] +1  ) , nmax[1]/100 ,nmax[1]/100 )
  colssig <- gplots::colorpanel(length(temp[,1]), "white", "red", "darkmagenta")

  # Plot dots 
  for(i in 1:nr){ points( 1:nc, rep(i,nc),  cex= log(as.numeric(pp[i,])+1),  bg=  colssig[ match(round( -log10(padj[i,]) + 1 )  , round(temp[,1]) ) ] ,  col= 1*(padj[i,]<0.05) ,   pch=21  ) }
 

  if( !is.null(dim(colcol))   ){
     colcol = colcol[,filt_col]
     nnr = dim(colcol)[1]
      for(i in 1:nnr){ 
          points(  x=1:nc, y=rep(nr+1.5+i,nc), col = colcol[i,],  pch=15, cex=3, xpd=T)    
      } 
  }
    
  if( !is.null(dim(rowcol))   ){
      rowcol = rowcol[,filt_row] 
      nnr = dim(rowcol)[1]
      for(i in 1:nnr){ 
         points(  y=1:nr, x=rep(-0.5-i,nr), col = rowcol[i,],  pch=15, cex=3, xpd=T)
       } 
    }
  } 


m =  match(rownames(padj_all), colpals2[,9])
rowcols = colpals2[m,10]
rowcols[is.na(rowcols)] = "grey"

pdf("sex_degs_enrich_markers_dotplot.pdf")
dotplot_extra_dendro( padj_all, padj_all, p_all , rownames(padj_all), colnames(padj_all),  filt_row = which(filt2 & !is.na(padj_all[,1])) , filt_col = 1:dim(padj_all)[2] , rowcol= rowcols, colcol= colcols ) 
dev.off() 

celltypes =  gsub("wilcox-down |wilcox-up |wilcox ", "", colnames(padj))
celltypes =  gsub("mast-down |mast-up |mast ", "", colnames(padj))
m =  match(celltypes, colpals2[,9])
colcols2 = colpals2[m,10]
colcols2[is.na(colcols2)] = "grey"
dotplot_extra_dendro( padjaz, pvalsaz, paz*2 , rownames(padjaz), colnames(padjaz),  filt_row =   1:dim(padjaz)[1]  , filt_col = 1:dim(padjaz)[2] , rowcol= "", colcol= colcols2 ) 
dotplot_extra_dendro( padj3, pvals3, p3 , rownames(padj3), colnames(padj3),  filt_row =   1:dim(padj3)[1]  , filt_col = 1:dim(padj3)[2] , rowcol= "", colcol= colcols2 ) 

dotplot_extra_dendro( padjaz, pvalsaz, paz*2 , rownames(padjaz), colnames(padjaz),  filt_row =   1:dim(padjaz)[1]  , filt_col = filt , rowcol= "", colcol= colcols2 ) 

  
skip = c("Eryth", "HSPC", "Platelet")
filt =  is.na(match(  celltypes, skip) )   
filt2 =  (rowSums(pp)  > 0 ) 
tp = rowSums(pp)
tp[ grep("Xi", rownames(pp))  ]  = 0 
tp[14] = 1 # escape genes
tp[34] = 1 # whlbld specific escape
filt3  = tp > 0



gene_set_cols  = rep("#692C89", dim(padj)[1] ) 
gene_set_cols[grep("chr", rownames(padj)) ] = "#FDD700"
gene_set_cols[grep("age", rownames(padj)) ] = "grey"
gene_set_cols[grep("pop", rownames(padj)) ] = "#6ABD45"
gene_set_cols[grep("sex|SEX", rownames(padj)) ] = "#8A181A"
 
 hc = heatmap.3(-log10(padj[ filt2 ,filt] ), col=cols6, cexRow=0.5 , RowSideCol= gene_set_cols[filt2], ColSideCol=colcols2[filt] )
 hc = heatmap.3(-log10(padj[ filt2 ,filt] ), col=cols6, cexRow=0.5 , RowSideCol= gene_set_cols[filt2], ColSideCol=colsex2[filt] )

 
dotplot( padj, pvals, pp , rownames(padj), colnames(padj),  filt_row =  which(filt2)[hc$rowInd ] , filt_col =  which(filt)[hc$colInd] , 
rowcol=gene_set_cols, colcol= colcols2 ) 
 
colsex2 = rep( viridis(11)[5], length(colcols2)) 
colsex2[grep("up", colnames(padj))] = colsex[1]
colsex2[grep("down", colnames(padj))] = colsex[2]

pdf("sex_degs_enrich_dotplot.pdf", width=15, height=10)
dotplot_extra( padj, pvals, pp , rownames(padj), celltypes,  filt_row =  which(filt2)[hc$rowInd ] , filt_col =  which(filt)[hc$colInd] , rowcol= rbind(gene_set_cols, rep("white", length(gene_set_cols))), colcol= rbind(colcols2, colsex2) ) 
 dev.off() 


pdf("sex_degs_enrich_dotplot_dendro.pdf", width=15, height=10)
dotplot_extra_dendro( padj, pvals, pp , rownames(padj), celltypes,  filt_row = which(filt2), filt_col =  which(filt), 
rowcol= rbind(gene_set_cols, rep("white", length(gene_set_cols))), colcol= rbind(colcols2, colsex2) , rowden = T , colden = T ) 

 dev.off() 

pdf("sex_degs_enrich_dotplot_dendro_wilcox.pdf", width=15, height=10)
dotplot_extra_dendro( padj, pvals, pp , rownames(padj), celltypes,  filt_row = which(filt2), filt_col =  which(filt), rowcol= rbind(gene_set_cols, rep("white", length(gene_set_cols))), colcol= rbind(colcols2, colsex2) , rowden = T , colden = T )
 dev.off()  

pdf("sex_degs_enrich_dotplot_dendro_wilcox2.pdf", width=10, height=6)
dotplot_extra_dendro( padj, pvals, pp , rownames(padj), celltypes,  filt_row = which(filt3), filt_col =  which(filt), rowcol= rbind(gene_set_cols, rep("white", length(gene_set_cols))), colcol= rbind(colcols2, colsex2) , rowden = T , colden = T ) 
 dev.off()  


 dotplot_extra_dendro( padj, pvals, pp , rownames(padj), celltypes,  filt_row = 1:length(gene_set_cols), filt_col =  1:dim(padj)[2], rowcol= rbind(gene_set_cols, rep("white", length(gene_set_cols))), colcol= rbind(colcols2, colsex2) , rowden = F , colden = F ) 



dotplot_extra_dendro( padj, pvals, pp , rownames(padj), celltypes,  filt_row = which(filt3), filt_col =  which(filt), rowcol= rbind(gene_set_cols, rep("white", length(gene_set_cols))), colcol= rbind(colcols2, colsex2) , rowden = T , colden = T ) 
 


filt4 = is.na( match(rownames(paz), skip)   ) 
dotplot_extra_dendro( padjaz, pvalsaz, paz , rownames(padjaz), colnames(padjaz),  filt_row = filt4  , filt_col = filt , rowcol= "", colcol= colcols2 ) 
 
```


### Panel F - density plots of genes 
```{r}
library(Nebulosa) 
filt = !is.na(obj$predicted.celltype.l3) & obj$predicted.celltype.l3 != "Doublet"
temp = obj[,filt] 


png("neb_umapXY_RPS4X.png", height=2000, width=2000)
plot_density(temp,  "RPS4X"   )
dev.off() 

png("neb_umapXY_XIST.png", height=2000, width=2000)
plot_density(temp,  "XIST"   )
dev.off() 

png("neb_umapXY_RPS4Y1.png", height=2000, width=2000)
plot_density(temp,  "RPS4Y1"   )
dev.off() 

png("neb_umapXY_IL2RG.png", height=2000, width=2000)
plot_density(temp,  "IL2RG"   )
dev.off() 

png("neb_umapXY_CD99.png", height=2000, width=2000)
plot_density(temp,  "CD99"   )
dev.off() 


png("neb_umapXY_CYBB.png", height=2000, width=2000)
plot_density(temp,  "CYBB"   )
dev.off() 

png("neb_umapXY_CCL5.png", height=2000, width=2000)
plot_density(temp,  "CCL5"   )
dev.off() 


genes_of_interest  = c("XIST" , "RPS4X", "RPS4Y1" , "CD99",  "IL2RG" ,   "CYBB", "CCL5",  "CXCR4", "CD79A", "GZMK", "NKG7")

# Panel G - violin plots of genes 
test <- factor( temp$predicted.celltype.l3, levels = colpals2[,9] )
temp$predicted.celltype.l3 <- test

pdf("sex_degs_vln_plots.pdf", width=15)
VlnPlot(temp  , features= rev(genes_of_interest),stack=T,  group.by = "predicted.celltype.l3", split.by = "sex" , col=colsex,pt.size=0)
dev.off() 
```









