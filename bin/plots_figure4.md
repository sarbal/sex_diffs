# Figure 4: Sex-linked cis-eQTLs and sex-interacting eQTLs

### Panel B - total counts of eqtls across male, female and joint 
```{r}
load("palette2.Rdata")
source("helper.r")
load("eqtls_so_far.Rdata")
sbgenes = unique(sb_cis[,2])
func <- function(x) { length(unique(x))}

skip = c("Eryth", "Platelet", "HSPC", "ILC", "CD4 Proliferating", "CD8 Proliferating")
filt = is.na(match(eqtls[,1] , skip))
eqtls = eqtls[filt,]

filt = is.na(match(sb_cis[,1] , skip))
sb_cis = sb_cis[filt,]

filt = is.na(match(chrxeqtls[,1] , skip))
chrxeqtls = chrxeqtls[filt,]

filt = is.na(match(male_eqtls[,1] , skip))
male_eqtls = male_eqtls[filt,]

filt = is.na(match(female_eqtls[,1] , skip))
female_eqtls = female_eqtls[filt,]

alleqtls = rbind(eqtls, chrxeqtls) 
alleqtls = alleqtls[(alleqtls$Chromosome !=""),] 

esnps = plyr::count(alleqtls[,1]) 
egenes = tapply( alleqtls[,2], alleqtls[,1], func)
 
esnps_sb = plyr::count(sb_cis[,1]) 
egenes_sb = tapply( sb_cis[,2], sb_cis[,1], func)
 
temp = gsub("_", " ", names(egenes))
temp = gsub("NK CD", "NK_CD", temp)

m = match(colpals2[,9], temp)
f.c = !is.na(m)
f.t = m[f.c]

esnps_female = plyr::count(female_eqtls[,1]) 
egenes_female = tapply( female_eqtls[,3], female_eqtls[,1], func)

esnps_male = plyr::count(male_eqtls[,1]) 
egenes_male = tapply( male_eqtls[,3], male_eqtls[,1], func)

temp = plyr::count( cbind(alleqtls$Chromosome, alleqtls$cell_type) )
chrmat_snp = spread(temp, key="x.2", value="freq") 
chrmat_snp[is.na(chrmat_snp)] = 0 

temp = plyr::count( unique(cbind(alleqtls$Chromosome, alleqtls$cell_type, alleqtls$GeneID))[,1:2] )
chrmat_gene = spread(temp, key="x.2", value="freq") 
chrmat_gene[is.na(chrmat_gene)] = 0 

chrmat_snp[chrmat_snp[,1] == "X",1] = 23 
chrmat_gene[chrmat_gene[,1] == "X",1] = 23

# Joint analysis 
pdf("joint_analysis_eqtls.pdf")
par(mfrow=c(1,2))
barplot(esnps[f.t,2], hor=T, col = colpals2[f.c,10], xlab="Number of eSNPs", border = NA)
barplot(egenes[f.t], hor=T, col = colpals2[f.c,10], xlab="Number of egenes", border = NA)
dev.off() 

# Stratified analysis 
pdf("female_analysis_eqtls.pdf")
par(mfrow=c(1,2))
barplot(esnps_female[f.t,2], hor=T, col = colpals2[f.c,10], xlab="Number of eSNPs", border = NA)
barplot(egenes_female[f.t], hor=T, col = colpals2[f.c,10], xlab="Number of egenes", border = NA)
dev.off() 

pdf("male_analysis_eqtls.pdf")
par(mfrow=c(1,2))
barplot(esnps_male[f.t,2], hor=T, col = colpals2[f.c,10], xlab="Number of eSNPs", border = NA)
barplot(egenes_male[f.t], hor=T, col = colpals2[f.c,10], xlab="Number of egenes", border = NA)
dev.off() 

# Sex-biased eqtls 
pdf("sb_eqtls.pdf")
par(mfrow=c(1,2))
barplot(esnps_sb[f.t2,2], hor=T, col = colpals2[f.c2,10], xlab="Number of eSNPs", border = NA)
barplot(egenes_sb[f.t2], hor=T, col = colpals2[f.c2,10], xlab="Number of egenes", border = NA)
dev.off() 
```

### Panel C - upset plots 

### Panel D - comparing rho estimates 

### Panel E - example plots  

### Panel F - differential expression overlap 
#### plotting violin plots of genes FCGR3A and ITGB2 
```{r}
# FCGR3A 
file = files[7]
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

     vioplot(res$expression[fc] ~ res$genotype[fc], ylim=yrange, main=paste( celltypei ,  geneid, rsid), xlab="Genotype", ylab="Residual", border= (colsex[3]) , col=0, colMed = c(colsex[3]), rectCol= make_transparent(colsex[3]) , lineCol= (colsex[3])) 
     abline(zlm, col=colsex[3],  lwd=2)
     vioplot(res$expression[fc & ff] ~ res$genotype[fc & ff], main="female", ylim=yrange,  xlab="Genotype", ylab="Residual" ,border= (colsex[2]) , col=0, colMed = c(colsex[2]), rectCol= make_transparent(colsex[2]) , lineCol= (colsex[2]))  
     abline(zlmf, col=colsex[2], lwd=2)
     vioplot(res$expression[fc & fm] ~ res$genotype[fc & fm], main="male", ylim=yrange,   xlab="Genotype", ylab="Residual", border= (colsex[1]) , col=0, colMed = c(colsex[1]), rectCol= make_transparent(colsex[1]) , lineCol= (colsex[1]))   
     abline(zlmm, col=colsex[1], lwd=2)
}
dev.off() 
```


### plotting volcano plots for NK and CD14 Monocytes  
```{r}
```
