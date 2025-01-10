quiet <- suppressPackageStartupMessages 

if( ! quiet(require("fossil")))
  install.packages("fossil")

if( ! quiet(require("gplots")))
  install.packages("gplots")

if( !quiet(require("plyr")))
  install.packages("plyr")

if( !quiet(require("RColorBrewer")))
  install.packages("RColorBrewer")

if( !quiet(require("viridis")))
  install.packages("viridis")

if( !quiet(require("beanplot")))
  install.packages("beanplot")

if( !quiet(require("beeswarm")))
  install.packages("beeswarm")

if( !quiet(require("corrgram")))
  install.packages("corrgram")

if( !quiet(require("vioplot")))
  install.packages("vioplot")

if( !quiet(require("pheatmap")))
  install.packages("pheatmap")

if( !quiet(require("ape")))
  install.packages("ape")

if( !quiet(require("venn")))
  install.packages("venn")

if( !quiet(require("tidyverse")))
  install.packages("tidyverse")

if( !quiet(require("Seurat")))
  install.packages("Seurat")

if( !quiet(require("igraph")))
  install.packages("igraph")

#if( !quiet(require("devtools")))
#  install.packages("devtools")
 
if (!quiet(require("BiocManager")))
  install.packages("BiocManager")

if( !quiet(require("EGAD")))
 BiocManager::install("EGAD")

if( !quiet(require("WGCNA")))
    BiocManager::install("WGCNA")

if( !quiet(require("DESeq2")))
  BiocManager::install("DESeq2")

if( !quiet(require("edgeR")))
  BiocManager::install("edgeR")

if( !quiet(require("rhdf5")))
  BiocManager::install("rhdf5")

if( !quiet(require("RSkittleBrewer")))
  install_github('alyssafrazee/RSkittleBrewer')

 if( !quiet(require("wesanderson")))
  install.packages("wesanderson")

if( !quiet(require("Nebulosa")))
  BiocManager::install("Nebulosa")

if( !quiet(require("cluster"))) 
  install.packages("cluster")
 
if( !quiet(require("data.table"))) 
  install.packages("data.table")
  

# Colors
cols = colorpanel(16, "red", "blue")
cols2 = brewer.pal(8, "Spectral")
cols3= rainbow(30)
cols4 = colorpanel(63, "lightgrey", "blue", "darkblue")
cols5 = colorpanel(300, "lightgrey", "red", "darkred")
cols6 = colorpanel(100, "lightgrey", "red", "darkmagenta")
cols7 = c("seagreen", "black", "darkmagenta")
cols8 = viridis(10)
cols9 = colorpanel(100, "white", "red", "darkmagenta")
cols10 = colorpanel(100, "white", "blue", "darkcyan")
cols11 = colorpanel(100, "white", "orange", "deeppink4")
cols12 = magma(100)
cols13 = viridis(100)
col_map= cols13
cols14 = c("grey", colorpanel(100, "black", "darkmagenta", "magenta") ) 

#pal <- wes_palette("Zissou1", 21, type = "continuous") 
original = RSkittleBrewer('original')
tropical = RSkittleBrewer('tropical')
wildberry = RSkittleBrewer('wildberry')
mm = RSkittleBrewer('M&M')

# Random functions
rank_std <- function(x)  { r = rank( x, na.last="keep"); r/max(r, na.rm=T)  }
colSD <- function( data){ return(apply( data, 2, sd, na.rm=T))}
rowSD <- function( data){ return(apply( data, 1, sd, na.rm=T))}
colSE <- function( data){ return( apply( data, 2, sd, na.rm=T)/sqrt(dim(data)[2]))}
rowSE <- function( data){ return( apply( data, 1, sd, na.rm=T)/sqrt(dim(data)[1]))}
se    <- function(x){ sd(x,na.rm=T)/sqrt(length(!is.na(x))) }
rmse  <- function(error){sqrt(mean(error^2, na.rm=T) )}
mae   <- function(error){ mean(abs(error), na.rm=T)}

geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

geo_sd <- function(data) {
  log_data <- log(data)
  gs <- exp(sd(log_data[is.finite(log_data)]))
  return(gs)
}

geo_se <- function(data) {
  gs <- geo_sd(data)
  log_data <- log(data)
  gse <- gs/sqrt(sum(is.finite(log_data)))
  return(gse)
}


lm.studentized  <- function(x,y){
  z = lm(y ~ x )
  z = rstudent(z)
  return( rank(abs(z)) )
}

lm.function  <- function(x,y){
  z = lm(y ~ x )
  return( rank(abs(z$residuals)) )
}


residuals <- function(x,y,A,B,C){ (A*x + B*y + C) }
residuals2 <- function(x,y,A,B,C) { (A*x + B*y + C)/sqrt(A^2+B^2) }
residuals3 <- function(x,y,A,B,C) { abs(A*x + B*y + C)/sqrt(A^2+B^2) }

z_scores <- function(x) {
  mean_x = mean(x, na.rm=T)
  sd_x = sd(x, na.rm=T)
  z =  (x - mean_x) / (sd_x)
  return(z)
}

z_scores_mod <- function(x) {
  med_x = median(x, na.rm=T)
  mad_x = median(abs(x-med_x), na.rm=T)
  z =  0.6745 * (x - med_x) / (mad_x)
  return(z)
}

calc_cpm <-function(X){
  K  = colSums(X)
  X.cpm = sapply(1:length(K), function(k) 10^6*X[,k]/K[k] )
  return(X.cpm)
}


heatmap.3 <- function(mat, ...){
  heatmap.2( mat, ..., density="none", trace="none")
}

vioplot.2 <- function(x,colin=cols3,plotpoints=FALSE,...){
  if( class(x) == "list") {
    temp = unlist(x)
    yrange = range(temp, na.rm=T)
    xmax = length(x) + 1
    plot(0,0, xlim=c(0,xmax), ylim=yrange, axes=F, col=0, ..., bty="n")
    axis(1, labels = names(x), at = 1:(xmax-1)  )
    for( i in 1:(xmax-1)) {
      vioplot(x[[i]][ !is.na(x[[i]])], col=colin[i], horizontal=FALSE, at=i, add=TRUE, rectCol="gray", axes=F)
      if(plotpoints==TRUE){  points( jitter( rep(i, length(x[[i]])),amount=0.1 ), x[[i]],   pch=19 , ... ) } 
    }
    axis(2)
  }
}



vioplot.3 <- function(x,colin=cols3,plotpoints=FALSE,xlab="", ylab="",...){
  if( class(x) == "list") {
    temp = unlist(x)
    yrange = range(temp, na.rm=T)
    xmax = length(x) + 1
    plot(0,0, ylim=c(0,xmax), xlim=yrange, axes=F, ylab=ylab, xlab=xlab, col=0, ..., bty="n")
    axis(2, labels = names(x), at = 1:(xmax-1) , las=2 , cex.axis=0.2)
    
    for( i in 1:(xmax-1)) {
      vioplot(x[[i]][ !is.na(x[[i]])], col=colin[i], horizontal=TRUE, at=i, add=TRUE, rectCol="gray", axes=F)
      if(plotpoints==TRUE){  points( x[[i]], jitter( rep(i, length(x[[i]])),amount=0.1 ),   pch=19 , ... ) } 
    }
    axis(1)
  }
}





# Transparent colors
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


# Given x and two points
get_value <- function( x1, x2, y1,y2, x) {
  m = (y2 - y1) / (x2 - x1 )
  y = y1 + m *( x - x1)
  return(y)
}

# Given y and two points
get_value_x <- function( x1, x2, y1,y2, y) {
  m = (y2 - y1) / (x2 - x1 )
  x = x1 + (y - y1)/m
  return(x)
}


## Formats the density distribution from the histogram function
get_density <- function(hist)
{
  x = sort(rep(hist$breaks,2))
  y = matrix(rbind(hist$density, hist$density))
  y = c(0,y,0)

  return(cbind(x,y))
}


## Formats the counts distribution from the histogram function
get_counts <- function(hist)
{
  x = sort(rep(hist$breaks,2))
  y = matrix(rbind(hist$counts, hist$counts))
  y = c(0,y,0)

  return(cbind(x,y))
}

# Tic toc functions
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self"))
{
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function()
{
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}

convolve_nets <- function(netA,netB,f){
  n <- order(netA)
  temp_netA <- netA[n]

  temp_netB = convolve( netB[n], rep(1,f),type="filter")
  convolved = cbind(temp_netA[(f/2):(length(temp_netA)-f/2)],temp_netB/f)
}

gene_set_enrichment <- function(genes, genes.labels, voc){

  genes.names = rownames(genes.labels)
  labels.names = colnames(genes.labels)
  genes.counts = rowSums(genes.labels)
  labels.counts = colSums(genes.labels)              			# p

  m = match ( genes, genes.names )
  filt.genes  = !is.na(m)
  filt.labels = m[filt.genes]

  labels.counts.set = rep( sum(filt.genes), length(labels.counts) )	# g

  m = match (labels.names, voc[,1])
  v.f = !is.na(m)
  v.g = m[v.f]

  universe = rep ( dim(genes.labels)[1], dim(genes.labels)[2])
  if(  length(filt.labels) == 1 ) { 
      genes.counts.set = genes.labels[filt.labels,]  
      test =  cbind(   (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set) # check again 
  } else { 
    genes.counts.set = colSums(genes.labels[filt.labels,]) 
    test =  cbind( (genes.counts.set -1) , labels.counts, universe-labels.counts, labels.counts.set)
  }             ## does weird things with 0 sets

  
  pvals = phyper(test[,1], test[,2], test[,3], test[,4], lower.tail=F)
  sigs = pvals < ( 0.05/length(pvals) )
  pvals.adj = p.adjust( pvals, method="BH")

  results = cbind(voc[v.g,1:2], test[v.f,c(1,2)], pvals[v.f], pvals.adj[v.f], sigs[v.f])
  colnames(results) = c("term", "descrp","p", "q", "pvals", "padj", "sig" )

  return (results)

}






dotplot_extra_dendro <- function(padj, pval , pp, rowids, colids, filt_row, filt_col, rowcol="", colcol="", rowden=F, colden=F) {
  rowids = rowids[filt_row]
  colids = colids[filt_col] 
  pval = pval[filt_row,filt_col]
  pp = pp[filt_row,filt_col]
  padj = padj[filt_row,filt_col]
  nc = length(colids)   
  nr = length(rowids) 
  
  if( !is.null(dim(colcol))){  colcol = colcol[filt_col,] } 
  if( !is.null(dim(rowcol))){  rowcol = rowcol[filt_row,] } 



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
    if( !is.null(dim(rowcol))){  rowcol = rowcol[hc_row$order,]} 

  } 
  if(colden==T)  { 
    hc_col = hclust( dist( t(padj) )) 
    padj = padj[,hc_col$order]
    pval = pval[,hc_col$order]
    pp   = pp[,hc_col$order]
    colids= colids[hc_col$order]
    filt_col = filt_col[hc_col$order]
    if( !is.null(dim(colcol))){  colcol = colcol[hc_col$order,]} 

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
     points(  x=1:nc, y=rep(nr*1.05 ,nc), col = colcol[,2],  pch=15, cex=3, xpd=T)
  }
  if( !is.null(dim(rowcol))   ){
     points(  y=1:nr, x=rep(-0.5,nr), col = rowcol[,2],  pch=15, cex=3, xpd=T)
      
    }
  }


