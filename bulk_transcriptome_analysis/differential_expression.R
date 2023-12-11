#########################################################################
############## Differential expression analysis script ##################
#########################################################################


### Déclarer le répertoire de travail { A REMPLIR }
setwd('')

#-----------------------------------------------------------------------#
#-----------------------Libraries installation--------------------------#
#-----------------------------------------------------------------------#
{
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
  BiocManager::install(c("limma,edgeR,Glimma"))
  install.packages(c("ggplot2","RColorBrewer","plotly","stringr","NMF",
                     "gplots","stringr","gprofiler2"))
}
}

#-----------------------------------------------------------------------#
#-------------------------Libararies import-----------------------------#
#-----------------------------------------------------------------------#
{
  library(limma,quietly=T)
  library(edgeR,quietly=T)
  library(Glimma,quietly=T)
  library(ggplot2,quietly=T)
  library(gplots,quietly=T)
  library(RColorBrewer,quietly=T)
  library(NMF,quietly=T)
  library(plotly,quietly=T)
  library(gprofiler2,quietly=T)
  library(stringr,quietly=T)
}

#-----------------------------------------------------------------------#
#-----------------------------Functions---------------------------------#
#-----------------------------------------------------------------------#
{
diff_expr <- function(countdata,metadata,conditions){
  metadata=subset(metadata,status==conditions[1] | status==conditions[2])
  countdata=countdata[,rownames(metadata)]
  
  metadata$bayes='cond1'
  metadata$bayes[metadata$status==conditions[2]]='cond2'
  y=DGEList(countdata)
  group=metadata$bayes
  group=factor(group)
  myCPM=cpm(countdata)
  thresh=myCPM>1
  keep=rowSums(thresh)>=2
  y=y[keep,keep.lib.sizes=FALSE]
  logcounts=cpm(y,log=TRUE)
  y <- calcNormFactors(y)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  v <- voom(y,design,plot = FALSE)
  fit <- lmFit(v)
  cont.matrix <- makeContrasts(CtrlVsCond=cond2 - cond1,levels=design)
  fit.cont <- contrasts.fit(fit, cont.matrix)
  fit.cont <- eBayes(fit.cont)
  summa.fit <- decideTests(fit.cont,adjust.method='bonferroni')
  
  #glXYPlot(x=fit.cont$coefficients[,1], y=fit.cont$lods[,1],
  #         xlab="logFoldChange", ylab="-log10(adjusted_pvalue)", main="P25 vs Starved",
  #         counts=v$E, groups=group, status=summa.fit[,1],
  #         anno=fit.cont$genes, side.main="gene_name", folder="volcano")
  
  count_norm=as.data.frame(cpm(y))
  dataCon1=count_norm[,c(which(metadata$bayes=='cond1'))]
  dataCon2=count_norm[,c(which(metadata$bayes=='cond2'))]
  foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
  temp=as.matrix(as.data.frame(summa.fit))
  pvals=fit.cont$p.value
  final=cbind(temp,foldChanges)
  final=cbind(final,pvals)
  colnames(final)[1]='DE_status'
  colnames(final)[3]='pvalue'
  return(as.data.frame(final))
}

expression_levels <- function(countdata,metadata){
  y=DGEList(countdata)
  group=metadata$status
  group=factor(group)
  myCPM=cpm(countdata)
  thresh=myCPM>1
  keep=rowSums(thresh)>=2
  y=y[keep,keep.lib.sizes=FALSE]
  logcounts=cpm(y,log=TRUE)
  y <- calcNormFactors(y)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  v <- voom(y,design,plot = FALSE)
}


gprofil <- function(genes,organism){
  ids=genes
  gostres = gost(query=ids,organism=organism)
  return(gostres)
}


plot_gprofil <- function(gostres){
  gostplot(gostres,capped=TRUE,interactive=TRUE)
}


visu_mds = function(countdata,metadata,title){
  group=metadata$status
  group=factor(group)
  col_colors=c('magenta','purple','orange','red','blue','darkgreen','green','black','cyan')
  assigns=col_colors[1:length(unique(metadata$status))]
  names(assigns)=unique(metadata$status)
  colors_used=c()
  for (i in 1:length(metadata$status)){
    new_color=assigns[[metadata$status[i]]]
    colors_used=c(colors_used,new_color)
  }
  y=DGEList(countdata, group=group)
  myCPM=cpm(countdata)
  thresh=myCPM>1
  keep=rowSums(thresh) >=2
  y <- y[keep, keep.lib.sizes=FALSE]
  logcounts <- cpm(y,log=TRUE)
  labels=metadata$status
  plotMDS(y,labels=labels,cex=0.8,main=title,col=colors_used)
}

impact_GO = function(genes_list,GOtype,organism){
  df=gprofil(genes_list,organism)$result
  df=df[df$source==GOtype,]
  ggplot(df, aes(y=reorder(str_wrap(term_name,width=90),intersection_size),x=intersection_size,fill=p_value))+
    geom_bar(stat='identity') + xlab('Nombre de gènes') +
    #theme(axis.text.y = element_text(size = 5)) +
    ylab(GOtype)
}
}

#-----------------------------------------------------------------------#
#-----------------------Déclarer les données----------------------------#
#-----------------------------------------------------------------------#

### Expression matrix file name
path_data=''

### Metadata file name
path_metadata=''

### Load data
countdata=read.table(path_data,header=T,row.names=1,sep=',',stringsAsFactors=T)
metadata=read.table(path_metadata,header=T,row.names=1,sep='\t',stringsAsFactors=T)


#-----------------------------------------------------------------------#
#----------------------------Data analysis------------------------------#
#-----------------------------------------------------------------------#

##### MDS plot #####
visu_mds(countdata,metadata,'')


#### Differential expression between two conditions #####
# Conditions to compare
cond1=''
cond2=''
conditions=c(cond1,cond2)

# Compare conditions
all_genes=diff_expr(countdata,metadata,conditions)
de_genes=all_genes[all_genes$DE_status!=0,]

# Export DE genes
new_file_name=''
write.csv(de_genes,new_file_name)


# Separate over- and under-expressed genes
up_genes=de_genes[de_genes$logFC>0,]
down_genes=de_genes[de_genes$logFC<0,]

# Gene names
up_genes_names=rownames(all_genes[all_genes$logFC>0,])
down_genes_names=rownames(all_genes[all_genes$logFC<0,])



##### Enrichment analysis #####

#Barplots
GO_term='' # GO:BP or GO:CC or GO:MF
model='' #hsapiens or celegans
impact_GO(up_genes_names,GO_term,model)
impact_GO(down_genes_names,Go_term,model)


#Manhattan plots
plot_gprofil(gprofil(up_genes_names))
plot_gprofil(gprofil(down_genes_names))



##### Plot expression of specific gene #####
v=expression_levels(countdata,metadata)
gene=''
plot_name=''

stripchart(v$E[gene,]~group,vertical=TRUE,las=1,cex.axis=0.8,pch=16,cex=0.8,
           ylab="Normalized log2 expression",main=plot_name)
boxplot(v$E[gene,]~group,vertical=TRUE,las=1,cex.axis=0.8,pch=16,cex=0.8,
        ylab="Normalized log2 expression", main=plot_name)
