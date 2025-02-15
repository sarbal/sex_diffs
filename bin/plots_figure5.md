# Fig. 5: Sex chromosome eQTL analysis. 

### Panel A - illustrations, using ideograms and karyotype 
```{r}
BiocManager::install("karyoploteR")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(karyoploteR)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
all.genes <- genes(txdb)

kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chrX", "chrY"))
kp <- kpPlotDensity(kp, all.genes)
kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chrX", "chrY")) 
kp <- kpPlotDensity(kp, all.genes, window.size = 1e4 )


escape = read.table("escape_genes_tuki2017")
all.genes <- genes(txdb, columns = "GENEID")
escape_entrez = EGAD::attr.human$entrezID[match( escape[,1], EGAD::attr.human$name )]
x_genes = subset(all.genes, GENEID %in% escape_entrez)

kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chrX", "chrY"))
kp <- kpPlotDensity(kp, x_genes, window.size = 1e5)

kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chrX", "chrY"))
kp <- kpPlotDensity(kp, x_genes, window.size = 1e4, col="darkblue")
kp <- kpPlotDensity(kp, all.genes, window.size = 1e4, col="darkblue")
   

kpDataBackground(kp, data.panel = 2)
kpDataBackground(kp, data.panel = 1)
kpPlotDensity(kp, all.genes, window.size = 1e4, col="darkblue")
kpPlotDensity(kp, x_genes, window.size = 1e4, col="red", data.panel=2) 
kpAxis(kp,cex=0.5, data.panel=2, side=1)
kpAxis(kp,cex=0.5, data.panel=1, side=1)
kpPlotMarkers(kp, x_genes, labels= attr.human$name[match( labels(x_genes), attr.human$entrezID ) ] ) 


kp <- plotKaryotype(genome="hg19", plot.type=4, chromosomes=c("chrX", "chrY"))
kpDataBackground(kp, data.panel = 2)
kpDataBackground(kp, data.panel = 1)
kpPlotDensity(kp, all.genes, window.size = 1e4, col="darkblue")
kpPlotDensity(kp, x_genes, window.size = 1e4, col="red", data.panel=2) 
 



kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chrX", "chrY"))
kpDataBackground(kp, data.panel = 2)
kpDataBackground(kp, data.panel = 1)
kpPlotCoverage(kp, all.genes,  col="darkblue")
kpPlotCoverage(kp, x_genes,  col="red", data.panel=2) 
kpAxis(kp,cex=0.5, data.panel=2, side=1)
kpAxis(kp,cex=0.5, data.panel=1, side=1)

degs = c("RPS4Y1",
"DDX3Y",
"EIF1AY",
"RPL39",
"TSC22D3",
"CD99",
"KDM5D",
"RPL36A",
"TTTY15",
"RPS4X",
"XIST",
"EIF1AX",
"EIF2S3",
"SEPT6",
"FLNA",
"IL2RG",
"PLP2","SAT1","SMC1A","SYAP1",
"ZRSR2")


#degs = c("CD99","RPS4X","XIST","EIF1AX","EIF2S3","SMC1A","SYAP1","ZRSR2")
degs2 = attr.human$entrezID[match( degs, attr.human$name ) ]

x_genes2 = subset(x_genes, gene_id %in% degs2)
y_genes2 = subset(y_genes, gene_id %in% degs2)

kp <- plotKaryotype(genome="hg19", plot.type=2, chromosomes=c("chrX", "chrY"))
kpPlotDensity(kp, x_genes, window.size = 1e4, col="darkblue")
kpPlotDensity(kp, y_genes, window.size = 1e4, col="darkblue")

kpPlotMarkers(kp, x_genes2, labels= attr.human$name[match( labels(x_genes2), attr.human$entrezID ) ], data.panel=2 ) 
kpPlotMarkers(kp, y_genes2, labels= attr.human$name[match( labels(y_genes2), attr.human$entrezID ) ], data.panel=2 ) 
```




### Panel B - illustrations 

### Panel C - PAR results 
```{r}
chrxeqtls =  read.table(file="chrX_cis_eQTLs_20220314_v2.tsv", header=T, sep="\t")
chrxeqtls0 =  read.table("chrX_cis_eQTLs_20220314.tsv", header=T, sep="\t")

chrxeqtls0 = chrxeqtls0[(chrxeqtls0$Chromosome !=""),] 
chrxeqtls = chrxeqtls[(chrxeqtls$Chromosome !=""),] 

chrmat_snp[chrmat_snp[,1] == "X",1] = 23 
chrmat_gene[chrmat_gene[,1] == "X",1] = 23

os = order(as.numeric(chrmat_snp[,1] )) 
og = order(as.numeric(chrmat_gene[,1] )) 
 
par(mfrow=c(1,2))
barplot(t(chrmat_gene[og,-1][,f.t]), names=chrmat_gene[og,1], hor=T, col = colpals2[f.c,10], xlab="Number of egenes", border = NA)
barplot(t(chrmat_snp[os,-1][,f.t]), names=chrmat_snp[os,1], hor=T, col = colpals2[f.c,10], xlab="Number of eSNPs", border = NA)

 
par(mfrow=c(1,2))
barplot(t(chrmat_gene[og,-1][23,f.t]),names = colpals2[f.c,9], hor=T, col = colpals2[f.c,10], xlab="Number of egenes", border = NA, beside=T)
barplot(t(chrmat_snp[os,-1][23,f.t]),  names = colpals2[f.c,9], hor=T, col = colpals2[f.c,10], xlab="Number of eSNPs", border = NA, beside=T)
```

 

### Panel D
```{r}
pdf("sex_biased_eqtls_vln_plots.pdf", width=15)
VlnPlot(temp, features=sbgenes, stack=T, group.by = "predicted.celltype.l3", split.by = "sex" , col=colsex, pt.size=0)     
dev.off() 

pdf("vln_chrxrec_unique.pdf", width=15)
VlnPlot(temp, features= chrxgenes, stack=T, group.by = "predicted.celltype.l3", split.by = "sex" , col=colsex, pt.size=0)     
dev.off() 


colsex = c("#F7A51D","#671D5B","#2C788F")
 
files = dir() 

library(vioplot)

for(file in files){ 
  print(file)
  res = read.table(file, header=T)
  pdffile = gsub(".tsv", ".pdf", file)
  pdf(pdffile, width=20, height=16)
  celltypes = sort(unique(res$cellType) )
  ff = res$sex == "female" 
  fm = res$sex == "male" 
    par(mfrow=c(7,9) ) 
  for(celltypei in celltypes)  { 
      #par(mfrow=c(1,3) ) 
      fc = res$cellType == celltypei   

      

      zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
      zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
      zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
      yrange = range(res$expression[fc]) 

      vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=celltypei, xlab="Genotype", ylab="Residual", border= (1) , col=0, colMed = c(1), rectCol= (1) , lineCol= (1)) 
      abline(zlm, lwd=2)
      vioplot(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= (colsex[2]) , lineCol= (colsex[2]))  
      abline(zlmf, col=colsex[2], lwd=2)
      vioplot(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= (colsex[1]) , lineCol= (colsex[1]))   
      abline(zlmm, col=colsex[1], lwd=2)
  }
    dev.off() 

}



colsex = c("#F7A51D","#671D5B","#2C788F")
 
sbs = read.table("sb_cis_eQTLs_at_25perc_fdr_20220928.tsv", header=T ) 
 
pdffile = "sb_cis_eQTLs_at_25perc_fdr_20220928.pdf"


pdf(pdffile, width=8, height=5)
  
for(j in 1:length(sbs[,1]) ){ 
  #if( j == 28) { next; } 
  if( j == 44) { next; }
  if( j == 51) { next; }

  file = grep( sbs[j,2], files, val=T) 

  res = read.table(file, header=T)
  celltypes = sort(unique(res$cellType) )
  ff = res$sex == "female" 
  fm = res$sex == "male" 
  #  par(mfrow=c(7,9) ) 
  celltypei = sbs[j,1] 
  rsid = strsplit(file, "_")[[1]][2]
  geneid =  strsplit(file, "_")[[1]][1]
   par(mfrow=c(1,3) ) 
   fc = res$cellType == celltypei   
   zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
   zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
   zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
   yrange = range(res$expression[fc]) 

   vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=paste( celltypei ,  geneid, rsid), xlab="Genotype", ylab="Residual", border= (1) , col=0, colMed = c(1), rectCol= make_transparent(1) , lineCol= (1)) 
   abline(zlm, lwd=2)
   vioplot(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= make_transparent(colsex[2]) , lineCol= (colsex[2]))  
   abline(zlmf, col=colsex[2], lwd=2)
   vioplot(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= make_transparent(colsex[1]) , lineCol= (colsex[1]))   
   abline(zlmm, col=colsex[1], lwd=2)


}
 
    dev.off() 


pdffile = "sb_cis_eQTLs_at_25perc_fdr_20220928_missing.pdf"
pdf(pdffile, width=8, height=5)

j = 28 
file = grep( sbs[j,2], files, val=T) 


j = 44 
file = grep( sbs[j,2], files, val=T) [1] 
file = grep( sbs[j,2], files, val=T) [2] 

j = 51 
file = grep( sbs[j,2], files, val=T) [1] 
file = grep( sbs[j,2], files, val=T) [2] 


  res = read.table(file, header=T)
  celltypes = sort(unique(res$cellType) )
  ff = res$sex == "female" 
  fm = res$sex == "male" 
  celltypei = sbs[j,1] 
  rsid = strsplit(file, "_")[[1]][2]
  geneid =  strsplit(file, "_")[[1]][1]

   par(mfrow=c(1,3) ) 
   fc = res$cellType == celltypei   
   zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
   zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
   zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
   yrange = range(res$expression[fc]) 

   vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=paste( celltypei ,  geneid, rsid), xlab="Genotype", ylab="Residual", border= (1) , col=0, colMed = c(1), rectCol= make_transparent(1) , lineCol= (1)) 
    abline(zlm, lwd=2)
   vioplot(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= make_transparent(colsex[2]) , lineCol= (colsex[2]))  
   abline(zlmf, col=colsex[2], lwd=2)
   vioplot(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= make_transparent(colsex[1]) , lineCol= (colsex[1]))   
   abline(zlmm, col=colsex[1], lwd=2)


dev.off() 




  file = grep( sbs[j,2], files, val=T) 

  res = read.table(file, header=T)
  celltypes = sort(unique(res$cellType) )
  ff = res$sex == "female" 
  fm = res$sex == "male" 
  celltypei = sbs[j,1] 

   par(mfrow=c(1,3) ) 
   fc = res$cellType == celltypei   
   zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
   zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
   zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
   yrange = range(res$expression[fc]) 

   boxplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=paste( celltypei , sbs[j,2]), xlab="Genotype", ylab="Residual", border= (1) , col=0, colMed = c(1), rectCol= make_transparent(1) , lineCol= (1)) 
   abline(zlm, lwd=2)
   boxplot(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= make_transparent(colsex[2]) , lineCol= (colsex[2]))  
   abline(zlmf, col=colsex[2], lwd=2)
   boxplot(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= make_transparent(colsex[1]) , lineCol= (colsex[1]))   
   abline(zlmm, col=colsex[1], lwd=2)


   par(mfrow=c(1,3) ) 
   fc = res$cellType == celltypei   
   zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
   zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
   zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
   yrange = range(res$expression[fc]) 

   beeswarm(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=paste( celltypei , sbs[j,2]), xlab="Genotype", ylab="Residual", border= (1) , col=colsex[3], colMed = c(1), rectCol= make_transparent(1) , lineCol= (1)) 
   abline(zlm, lwd=2)
   beeswarm(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=colsex[2], colMed = c(colsex[2]), rectCol= make_transparent(colsex[2]) , lineCol= (colsex[2]))  
   abline(zlmf, col=colsex[2], lwd=2)
   beeswarm(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=colsex[1], colMed = c(colsex[1]), rectCol= make_transparent(colsex[1]) , lineCol= (colsex[1]))   
   abline(zlmm, col=colsex[1], lwd=2)





  temp = do.call(rbind, strsplit(file, "_")  ) 

  genename = temp[ ,1] 
  rsid = temp[ ,2] 
  #celltypei = sbs[j,1]
  
  par(mfrow=c(1,3))

  ggplot(res[fc,], aes(x = genotype, y = expression, color=cellType)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ cellType, scales = "free", nrow=2) + 
  theme_classic() +
  #scale_color_manual(values=colsex[1], drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  geom_smooth(aes(group=cellType), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))


  ggplot(res[fc&fm,], aes(x = genotype, y = expression, color=cellType)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ cellType, scales = "free", nrow=2) + 
  theme_classic() +
  #scale_color_manual(values=colsex[1], drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  # strip.text.x = element_blank(),
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  geom_smooth(aes(group=cellType), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))



  ggplot(res[fc&ff,], aes(x = genotype, y = expression, color=cellType)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ cellType, scales = "free", nrow=2) + 
  theme_classic() +
  #scale_color_manual(values=colsex[1], drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  # strip.text.x = element_blank(),
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  geom_smooth(aes(group=cellType), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))




  ggplot(res, aes(x = genotype, y = expression, color=cellType)) + geom_violin(trim=FALSE) + #geom_bar(stat = "identity") + 
  geom_boxplot(width=0.1, outlier.size = 0.5, fill="white") +
  facet_wrap(~ cellType, scales = "free", nrow=2) + 
  theme_classic() +
  #scale_color_manual(values=colsex[1], drop=FALSE) +
  theme(legend.position='none',
  strip.background = element_blank(),
  plot.title = element_text(hjust = 0.5)) +
  # strip.text.x = element_blank(),
  #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  geom_smooth(aes(group=cellType), method="lm", fullrange=T, color="red",size=0.5, se=FALSE) +
  labs(title=paste0(genename,"(",rsid,")"))





  file = grep( sbs[j,2], files, val=T) 

  res = read.table(file, header=T)
  celltypes = sort(unique(res$cellType) )
  ff = res$sex == "female" 
  fm = res$sex == "male" 
  celltypei = sbs[j,1] 
  rsid = strsplit(file, "_")[[1]][2]
  geneid =  strsplit(file, "_")[[1]][1]
   par(mfrow=c(1,3) ) 
   fc = res$cellType == celltypei   
   zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
   zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
   zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
   yrange = range(res$expression[fc]) 

 temp3 =  cbind( cor(res$expression[fc] ,  as.numeric(as.factor(res$genotype[fc] ) )) ,
   cor(res$expression[fc & ff],  as.numeric(as.factor(res$genotype[fc & ff]))) ,
   cor(res$expression[fc & fm], as.numeric(as.factor(res$genotype[fc & fm]) ))) 

colnames(temp3) = c("joint", "female", "male")
 rownames(temp3) = file 






temp3 = list() 
for(j in 1:length(sbs[,1]) ){ 
  #if( j == 28) { next; } 
  if( j == 44) { next; }
  if( j == 51) { next; }

  file = grep( sbs[j,2], files, val=T) 

  res = read.table(file, header=T)
  celltypes = sort(unique(res$cellType) )
  ff = res$sex == "female" 
  fm = res$sex == "male"
  fj = res$sex == "joint"
  #  par(mfrow=c(7,9) ) 
  celltypei = sbs[j,1] 
   fc = res$cellType == celltypei   
 rsid = strsplit(file, "_")[[1]][2]
  geneid =  strsplit(file, "_")[[1]][1]
  temp3[[j]] = cbind( cor(res$expression[fc & fj] ,  as.numeric(as.factor(res$genotype[fc &fj] ) ), m="s") ,
   cor(res$expression[fc & ff],  as.numeric(as.factor(res$genotype[fc & ff])), m="s") ,
   cor(res$expression[fc & fm], as.numeric(as.factor(res$genotype[fc & fm]) ), m="s"), 
cor(res$expression[fc & fj] ,  as.numeric(as.factor(res$genotype[fc &fj] ) ), m="p") ,
   cor(res$expression[fc & ff],  as.numeric(as.factor(res$genotype[fc & ff])), m="p") ,
   cor(res$expression[fc & fm], as.numeric(as.factor(res$genotype[fc & fm]) ), m="p"), 
   sum(fc&fj), sum(fc&ff), sum(fc & fm)   ) 

}
 
temp5 = list() 
for(j in 1:length(sbs[,1]) ){ 
  #if( j == 28) { next; } 
  if( j == 44) { next; }
  if( j == 51) { next; }

  file = grep( sbs[j,2], files, val=T) 

  res = read.table(file, header=T)
  celltypes = sort(unique(res$cellType) )
  ff = res$sex == "female" 
  fm = res$sex == "male"
  fj = res$sex == "joint"
  celltypei = sbs[j,1] 
  fc = res$cellType == celltypei   
  rsid = strsplit(file, "_")[[1]][2]
  geneid =  strsplit(file, "_")[[1]][1]

 
  test1 = aov(res$expression[fc & fj] ~ as.numeric(as.factor(res$genotype[fc &fj]))) 
  test2 = aov(res$expression[fc & ff] ~ as.numeric(as.factor(res$genotype[fc & ff]))) 
  test3 = aov(res$expression[fc & fm] ~ as.numeric(as.factor(res$genotype[fc & fm])))  
  temp5[[j]] = cbind( test1$coefficients[2] , test2$coefficients[2], test3$coefficients[2])  

}  
 
temp6 = lapply( 1:51, function(i) spread ( do.call(rbind, lapply(1:3, function(j) cbind(j, temp4[[i]][[j]]) ) ), key=1, value=3 ) )   
 
checks = c(7,10,15,39,40)

 


j = 51 
file = "FAM118A_rs104664_plot_data.tsv"  
celltypei="Plasmablast"   



j = 44 
file = "FAM118A_rs170505_plot_data.tsv"  
celltypei="gdT"
 

 





file = files[99]
pdffile = gsub(".tsv", "_2.pdf", file)
res = read.table(file, header=T)
   celltypes = sort(unique(res$cellType) )
  
 pdf(pdffile, width=8, height=5)
 
 for(celltypei in celltypes)  { 
    print (celltypei)
    ff = res$sex == "female" 
    fm = res$sex == "male" 
    rsid = strsplit(file, "_")[[1]][2]
    geneid =  strsplit(file, "_")[[1]][1]
    par(mfrow=c(1,3) ) 
     fc = res$cellType == celltypei   
     zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
     zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
     zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
     yrange = range(res$expression[fc]) 

     vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=paste( celltypei ,  geneid, rsid), xlab="Genotype", ylab="Residual", border= (1) , col=0, colMed = c(1), rectCol= make_transparent(1) , lineCol= (1)) 
     abline(zlm, lwd=2)
     vioplot(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= make_transparent(colsex[2]) , lineCol= (colsex[2]))  
     abline(zlmf, col=colsex[2], lwd=2)
     vioplot(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= make_transparent(colsex[1]) , lineCol= (colsex[1]))   
     abline(zlmm, col=colsex[1], lwd=2)
}
dev.off() 

 files = grep("tsv", dir("2022_11/plots_2/") , val=T)
 
 file = files[99]
pdffile = gsub(".tsv", "_3.pdf", file)
res = read.table(file, header=T)
   celltypes = sort(unique(res$cellType) )
  
 pdf(pdffile, width=8, height=5)
 
 for(celltypei in celltypes)  { 
    print (celltypei)
    ff = res$sex == "female" 
    fm = res$sex == "male" 
    rsid = strsplit(file, "_")[[1]][2]
    geneid =  strsplit(file, "_")[[1]][1]
    par(mfrow=c(1,3) ) 
     fc = res$cellType == celltypei   
     zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
     zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
     zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
     yrange = range(res$expression[fc]) 

     vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=paste( celltypei ,  geneid, rsid), xlab="Genotype", ylab="Residual", border= (colsex[3]) , col=0, colMed = c(colsex[3]), rectCol= make_transparent(colsex[3]) , lineCol= (colsex[3])) 
     beeswarm(res$expression[fc] ~ res$genotype[fc],   add=T, col=colsex[3] , pch=19, cex=0.5)    
     abline(zlm, lwd=2)
     vioplot(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= make_transparent(colsex[2]) , lineCol= (colsex[2]))  
     beeswarm(res$expression[fc & ff] ~ res$genotype[fc & ff],  col=(colsex[2]) , add=T, pch=19, cex=0.5)     
     abline(zlmf, col=colsex[2], lwd=2)
     vioplot(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= make_transparent(colsex[1]) , lineCol= (colsex[1]))   
     beeswarm(res$expression[fc & fm] ~ res$genotype[fc & fm],col= (colsex[1]) , add=T, pch=19, cex=0.5)      
     abline(zlmm, col=colsex[1], lwd=2)
}
dev.off() 

 



 files = grep("tsv", dir("2022_11/plots_2/") , val=T)
  i = grep("CD52", files )
 

 file = files[i] 
pdffile = gsub(".tsv", "_2.pdf", file)
res = read.table(file, header=T)
   celltypes = sort(unique(res$cellType) )
  
 pdf(pdffile, width=8, height=5)
 
 for(celltypei in celltypes)  { 
    print (celltypei)
    ff = res$sex == "female" 
    fm = res$sex == "male" 
    rsid = strsplit(file, "_")[[1]][2]
    geneid =  strsplit(file, "_")[[1]][1]
    par(mfrow=c(1,3) ) 
     fc = res$cellType == celltypei   
     zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
     zlmf = lm(res$expression[fc & ff] ~  as.numeric(as.factor(res$genotype[fc & ff]) ))
     zlmm = lm(res$expression[fc & fm] ~  as.numeric(as.factor(res$genotype[fc & fm]) ))
     yrange = range(res$expression[fc]) 

     vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=paste( celltypei ,  geneid, rsid), xlab="Genotype", ylab="Residual", border= (1) , col=0, colMed = c(1), rectCol= make_transparent(1) , lineCol= (1)) 
     abline(zlm, lwd=2)
     vioplot(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= make_transparent(colsex[2]) , lineCol= (colsex[2]))  
     abline(zlmf, col=colsex[2], lwd=2)
     vioplot(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= make_transparent(colsex[1]) , lineCol= (colsex[1]))   
     abline(zlmm, col=colsex[1], lwd=2)
}
dev.off() 
```








### Panel G 
```{r}
chrmlength = read.table(chrNameLength.txt")
chrmlength = chrmlength[grep("chr", (chrmlength[,1]) ) ,]
startpos = cumsum( as.numeric(chrmlength[,2]) ) 
startpos = c(0,startpos)
startpos = startpos[-26]
chrmlength = cbind(chrmlength, startpos) 

cellcols = colpals2[match(gsub("NK C", "NK_C", gsub("_", " ", eqtls$cell_type) ) , colpals2[,9] ),10]

alleqtlschrms = as.numeric(alleqtlschrms)
alleqtlschrms[alleqtls$Chromosome=="X"]  = 23 
alleqtlschrms = as.numeric(alleqtlschrms)

male_nk = read.table("male_chrX_correlations_all_regions_mod")
female_nk = read.table("female_chrX_correlations_all_regions_mod")

# manhattan plot, mising non-sig, might ask for them
plot( alleqtls$Position  + chrmlength[alleqtlschrms ,3] , -log10(as.numeric(alleqtls$FDR)), col = c("violet", "darkorchid")[alleqtlschrms%%2  + 1 ]  , pch=19  )
 


  
celltypes = unique(alleqtls$cell_type)
par(mfrow=c(5,1))
for( celltypei in celltypes[1:5]){ 
  filt = alleqtls$cell_type == celltypei
plot( alleqtls$Position[filt]  + chrmlength[alleqtlschrms ,3][filt] , -log10(as.numeric(alleqtls$FDR))[filt], col = c("violet", "darkorchid")[alleqtlschrms[filt]%%2  + 1 ]  , pch=19  )

}

   
eSNP_ranks = unique(alleqtls$eSNP_rank)
par(mfrow=c(5,1))
for( eSNP_rank in eSNP_ranks[1:5]){ 
  filt = alleqtls$eSNP_rank == eSNP_rank
plot( alleqtls$Position[filt]  + chrmlength[alleqtlschrms ,3][filt] , -log10(as.numeric(alleqtls$FDR))[filt], col = c("violet", "darkorchid")[alleqtlschrms[filt]%%2  + 1 ]  , pch=19  )

}



# NK specific list
library(vioplot)

files = dir() 


pdf("NK_X_results.pdf", width=8, height=5)
for(file in files){ 
  print(file)
  res = read.table(file, header=T)
  colnames(res)  = c("expression", "individual", "sex", "genotype")


  celltypei = "NK"
  geneid = gsub("plot_data_|.tsv", "", file)

  ff = res$sex == "female" 
  fm = res$sex == "male" 
  fc = res$sex == "joint"
      

      zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
      zlmf = lm(res$expression[ff] ~  as.numeric(as.factor(res$genotype[ff]) ))
      zlmm = lm(res$expression[fm] ~  as.numeric(as.factor(res$genotype[fm]) ))
      yrange = range(res$expression) 
     par(mfrow=c(1,3) ) 
      vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=geneid, xlab="Genotype", ylab="Residual", border= (colsex[3]) , col=0, colMed = c(colsex[3]), rectCol= (colsex[3]) , lineCol= (colsex[3])) 
      abline(zlm, col=colsex[3],lwd=2)
      vioplot(res$expression[ff] ~ res$genotype[ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= (colsex[2]) , lineCol= (colsex[2]))  
      abline(zlmf, col=colsex[2], lwd=2)
      vioplot(res$expression[fm] ~ res$genotype[fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= (colsex[1]) , lineCol= (colsex[1]))   
      abline(zlmm, col=colsex[1], lwd=2)
      
}
dev.off() 





pdf("sept6_X_results.pdf", width=8, height=5)
for(file in files){ 
  print(file)
  res = read.table(file, header=T)
  colnames(res)  = c("expression", "individual", "sex", "genotype")


  celltypei = "NK"
  geneid = gsub("plot_data_|.tsv", "", file)

  ff = res$sex == "female" 
  fm = res$sex == "male" 
  fc = res$sex == "joint"
      

      zlm = lm(res$expression[fc] ~  as.numeric(as.factor(res$genotype[fc])) )
      zlmf = lm(res$expression[ff] ~  as.numeric(as.factor(res$genotype[ff]) ))
      zlmm = lm(res$expression[fm] ~  as.numeric(as.factor(res$genotype[fm]) ))
      yrange = range(res$expression) 
     par(mfrow=c(1,3) ) 
      vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=geneid, xlab="Genotype", ylab="Residual", border= (colsex[3]) , col=0, colMed = c(colsex[3]), rectCol= (colsex[3]) , lineCol= (colsex[3])) 
      abline(zlm, col=colsex[3],lwd=2)
      vioplot(res$expression[ff] ~ res$genotype[ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= (colsex[2]) , lineCol= (colsex[2]))  
      abline(zlmf, col=colsex[2], lwd=2)
      vioplot(res$expression[fm] ~ res$genotype[fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= (colsex[1]) , lineCol= (colsex[1]))   
      abline(zlmm, col=colsex[1], lwd=2)
      
}
dev.off() 
```





