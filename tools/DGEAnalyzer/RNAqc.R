#############################################################################################################################
# 
# RNASeqQC v1.0
# Toolbox for loading RNAseq data and performing QC analyses
# (c) Maximilian Zwiebel, 2018
# 
#############################################################################################################################

library(LSD)
library(pheatmap)
library(limma)
library(edgeR)
library(DESeq2)
library(hexbin)
library(plot3D)

# load data of an RNAseq experiment and store it to an internal expression data object
# 
# dge_file: path to the DGE Matrix File ("DGE_Matrix.txt" by default)
# targets_file: path to targets file containing information about the experiment (e.g. group descriptors; "Targets.txt" by default)
# rowMin: minimal number of counts within one row (gene); the default value is 2
# 
# value: a list with three elements: "Targets" contains experimental information, rowMin the minimal number of counts within one row,
# "Matrix" the raw count matrix
loadData <- function(dge_file="DGE_Matrix.txt",targets_file="Targets.txt",rowMin=2){
	dge = read.csv(dge_file,header=T,sep='\t')
	# replace any '+' in the targets file by '.' ('+' in a name cannot be handled in R)
	targets = as.data.frame(apply(read.csv(targets_file,header=T,sep='\t'),2,function(x){x=gsub("\\+","\\.",x)}),stringsAsFactors=F)
	# if there are any fields with label 'Control', change these to 'control' (the default value for reference samples)
	for(i in 1:nrow(targets)){
	  for(j in 1:ncol(targets)){
	    if(targets[i,j] == "Control") targets[i,j] = "control"
	  }
	}
	colnames(targets)[1] = "Barcode"
	rownames(targets) = targets$Barcode
	
	rownames(dge) = dge$GENE
	dge = dge[,colnames(dge)!='GENE']
	dge = dge[,targets$Barcode]
	dge = dge[apply(dge,1,sum) >= rowMin,colnames(dge) %in% targets$Barcode]
	
	return(list(Targets=targets,rowMin=rowMin,Matrix=dge))
}


# transforms the raw count matrix of an expression data object by different methods (log, rlog, VST and CPM)
# 
# data: internal expression data object
# design: name of the column in the targets file that contains group information
# 
# value: internal expression data object with additional elements for each transformation and the design column
transformDGE <- function(data,design){
	data$log = log2(data$Matrix+1)
	
	dds <- DESeqDataSetFromMatrix(countData=data$Matrix,colData=data$Targets,as.formula(paste0("design=~",design)))
	data$rlog = assay(rlog(dds,blind=F))
	data$vst = assay(vst(dds,blind=F))
	data$cpm = cpm(data$Matrix,log=T)
	data$design = design
	
	return(data)
}


# show experimental design by printing the first few lines of the targets file
# 
# data: internal expression data object
# fields_only: logical indicating if only column names of the targets file should be shown; the default is FALSE
showDesign = function(data,fields_only=F){
	if(fields_only){
		return(colnames(data$Targets))
	} else {
		return(rbind(head(data$Targets,3),'...'))
	}
}


# rename samples by using informative columns from the targets information; if the columns are not sufficient to provide unique
# names, an enumeration for identical names is added; this can be done before or after transformation
# 
# data: internal expression data object
# fields: vector of column names to be considered for sample renaming
# 
# value: internal expression data object with changed sample names
createNames = function(data,fields){
	ids = c()
	for(i in rownames(data$Target)){
		ids = c(ids,paste(data$Targets[i,fields],collapse='_'))
	}
	if(anyDuplicated(ids)){
	  warning('Given fields are not sufficient to provide unique sample names. Enumeration is being added.')
	  for(i in ids[!duplicated(ids)]){
  	  nr=1
  	  for(j in which(ids==i)){
  	    ids[j] = paste0(ids[j],"_",nr)
  	    nr = nr+1
  	  }
  	  nr=1
  	}
	}
	names(ids) = data$Targets$Barcode
	
	rownames(data$Targets) = ids
	cn=ids[colnames(data$Matrix)]
	for(i in names(data)[!(names(data)%in%c('Targets','rowMin','design'))]){
		colnames(data[[i]]) = cn
	}
	
	return(data)
}


# plot dispersion of the data and its transformations
# 
# data: internal expression data object
dispersionPlot = function(data){
  x11()
  par(mfrow=c(2,2))
  heatscatter(log2(apply(data$Matrix,1,mean)),log2(apply(data$Matrix,1,var)),pch=20,main="original",xlab="log(µ)",ylab="log(s²)")
  abline(0,1,lty=2,col="red")
  heatscatter(apply(dat$log,1,mean),apply(dat$log,1,var),pch=20,main="log-transform",xlab="µ",ylab="s²")
  heatscatter(apply(dat$rlog,1,mean),apply(dat$rlog,1,var),pch=20,main="rlog-transform",xlab="µ",ylab="s²")
  heatscatter(apply(dat$vst,1,mean),apply(dat$vst,1,var),pch=20,main="VST",xlab="µ",ylab="s²")
}


# MA Plot
# 
# data: internal expression data object
# x: either specific sample name or group
# y: if x is a specific sample, y gives the name of a second sample to compare x to
# transformation: transformation to use; one of either raw matrix ("Matrix", the default), log, rlog, VST or CPM
MAplotter = function(data,x,y=NULL,transformation=c("log","rlog","vst","cpm")){
	transformation = match.arg(transformation)
	mat = data[[transformation]]
	x11()
	
	if(!is.null(y)){
	  if(!(x%in%colnames(mat))|!(y%in%colnames(mat))) stop('One of the specified samples does not exist.')
	  x = mat[,x]
	  y = mat[,y]
		M = x-y
		A = 0.5*(x+y)
		
		heatscatter(A,M,pch=20,xlab='mean(log2(expression))',ylab='log2 foldchange')
		abline(h=0,lty=2)
		
		df = data.frame(A,M)
		df = df[order(A),]
		lines(predict(loess(M~A,data=df,span=1,family='gaussian')),lwd=2,col='red')
	} else{
	  if(!(x%in%data$Targets[[data$design]])) stop('The specified group does not exist.')
	  targets = data$Targets[data$Targets[[data$design]]==x,]
	  size = dim(targets)[1]
	  
	  par(mfrow=rep(size,2),mar=rep(2,4))
	  mat = mat[,colnames(mat)%in%rownames(targets)]
	  for(i in 1:size){
	    for(j in 1:size){
	      if(j<i){
	        plot(x=0,y=0,xaxt='n',yaxt='n',xlab='',ylab='',type='n',frame.plot=F)
	      } else if (i==j){
	        plot(x=0,y=0,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',type='n',frame.plot=F)
	        text(0.5,0.5,colnames(mat)[j],cex=1.5)
        } else{
	        M = mat[,j]-mat[,i]
	        A = 0.5*(mat[,j]+mat[,i])
	        
	        heatscatter(A,M,pch=20,xlab='mean(log2(expression))',ylab='log2 foldchange',main='')
	        abline(h=0,lty=2)
	        #plot(hexbin(A,M),xlab='mean(log2(expression))',ylab='log2 foldchange',colramp=colorRampPalette(rev(c('blue','white','red','yellow'))))
	      }
	    }
	  }
	}
}


# QQ Plot
# 
# data: internal expression data object
# x: either specific sample name or group
# y: if x is a specific sample, y gives the name of a second sample to compare x to
# transformation: transformation to use; one of either raw matrix ("Matrix", the default), log, rlog, VST or CPM
QQplot = function(data,x,y=NULL,transformation=c("log","rlog","vst","cpm")){
  transformation = match.arg(transformation)
  mat = data[[transformation]]
  x11()
  
  if(!is.null(y)){
    if(!(x%in%colnames(mat))|!(y%in%colnames(mat))) stop('One of the specified samples does not exist.')
    xvalue = mat[,x]
    xvalue = (xvalue - mean(xvalue))/sd(xvalue)
    yvalue = mat[,y]
    yvalue = (yvalue - mean(yvalue))/sd(yvalue)
    
    qqplot(xvalue,yvalue,pch=20,xlab=x,ylab=y)
    abline(0,1,lty=2)
  } else{
    if(!(x%in%data$Targets[[data$design]])) stop('The specified group does not exist.')
    targets = data$Targets[data$Targets[[data$design]]==x,]
    size = dim(targets)[1]
    
    par(mfrow=rep(size,2),mar=rep(2,4))
    mat = mat[,colnames(mat)%in%rownames(targets)]
    for(i in 1:size){
      for(j in 1:size){
        if(j<i){
          plot(x=0,y=0,xaxt='n',yaxt='n',xlab='',ylab='',type='n',frame.plot=F)
        } else if (i==j){
          plot(x=0,y=0,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',type='n',frame.plot=F)
          text(0.5,0.5,colnames(mat)[j],cex=1.5)
        } else{
          xvalue = mat[,i]
          xvalue = (xvalue - mean(xvalue))/sd(xvalue)
          yvalue = mat[,j]
          yvalue = (yvalue - mean(yvalue))/sd(yvalue)
          
          qqplot(xvalue,yvalue,pch=20)
          abline(0,1,lty=2)
        }
      }
    }
  }
}


# Boxplot
# 
# data: internal expression data object
# x: either specific sample name or group identifier
# transformation: transformation to use; one of either log (the default), rlog, VST or CPM
# col: color to draw the density boxplots
boxplotDGE = function(data,x,transformation=c("log","rlog","vst","cpm"),col="white"){
	transformation = match.arg(transformation)
  mat = data[[transformation]]
	
  if(all(x%in%colnames(mat))){
		x11()
		boxplot(mat[,x],col=col,main='Readcount Distribution',ylab='log(expression)')
	} else if(x%in%data$Targets[[data$design]]){
	  targets = data$Targets[data$Targets[[data$design]]==x,]
	  mat = mat[,colnames(mat)%in%rownames(targets)]
	  
	  x11()
		boxplot(mat,col=col,main='Readcount Distribution',ylab='log(expression)')
	} else if(x=="all"){
	  col.code = as.numeric(as.factor(dat$Targets[[data$design]]))[order(colnames(mat))]
	  mat = mat[order(colnames(mat))]
	  color=rainbow(max(col.code))
	  
	  x11()
	  boxplot(mat,col=color[col.code],main='Readcount Distribution',ylab='log(expression)',las=2)
	} else stop("The specified sample / group does not exist.")
}


# Density Plot
# 
# data: internal expression data object
# x: either specific sample name or group identifier
# add: logical, indicating if new density plot should be added to the existing one (FALSE by default)
# transformation: transformation to use; one of either raw matrix ("Matrix", the default), log, rlog, VST or CPM
# col: color(group) to draw the density plot
densityPlot = function(data,x,add=F,transformation=c("log","rlog","vst","cpm"),col="black"){
	transformation = match.arg(transformation)
  mat = data[[transformation]]
	
  if(x%in%colnames(mat)){
    if(!add){
      x11()
      plot(x=NULL,y=NULL,xlim=c(0,max(mat[,x])),ylim=c(0,max(density(mat[,x])$y)+0.05),main='Readcount Densities',xlab='log(readcounts)',ylab='density',type='n')
    }
    lines(density(mat[,x]),col=col)
  } else if(x%in%data$Targets[[data$design]]){
    targets = data$Targets[data$Targets[[data$design]]==x,]
    size = dim(targets)[1]
    color=colorRampPalette(rev(c(col,"lightgray")))(size)
    mat = mat[,colnames(mat)%in%rownames(targets)]
    
    if(!add){
      x11()
      plot(x=NULL,y=NULL,xlim=c(0,max(mat)),ylim=c(0,max(apply(mat,2,function(x){return(density(x)$y)})+0.05)),main='Readcount Densities',xlab='log(readcounts)',ylab='density',type='n')
    }
    for(i in colnames(mat)){
      lines(density(mat[,i]),col=color[which(colnames(mat)==i)])
    }
  } else stop("The specified sample / group does not exist.")
}


# Cumulative Distribution
# 
# data: internal expression data object
# x: either specific sample name or group identifier
# add: logical, indicating if new density plot should be added to the existing one (FALSE by default)
# transformation: transformation to use; one of either raw matrix ("Matrix", the default), log, rlog, VST or CPM
# col: color(group) to draw the cumulative distribution
cumulDist = function(data,x,add=F,transformation=c("log","rlog","vst","cpm"),col="black"){
  transformation = match.arg(transformation)
  mat = data[[transformation]]
  
  if(x%in%colnames(mat)){
    if(!add){
      x11()
      plot(x=NULL,y=NULL,xlim=c(0,max(mat[,x])),ylim=c(0,1),main='Readcount Densities',xlab='log(readcounts)',ylab='density',type='n')
    }
    xvalue = mat[,x]
    xvalue = xvalue[order(xvalue)]
    cumul <- cumsum(xvalue)/sum(xvalue)
    
    lines(xvalue,cumul,col=col)
  } else if(x%in%data$Targets[[data$design]]){
    targets = data$Targets[data$Targets[[data$design]]==x,]
    size = dim(targets)[1]
    color=colorRampPalette(rev(c(col,"lightgray")))(size)
    mat = mat[,colnames(mat)%in%rownames(targets)]
    
    if(!add){
      x11()
      plot(x=NULL,y=NULL,xlim=c(0,max(mat)),ylim=c(0,1),main='Readcount Densities',xlab='log(readcounts)',ylab='density',type='n')
    }
    for(i in colnames(mat)){
      xvalue = mat[,i]
      xvalue = xvalue[order(xvalue)]
      cumul <- cumsum(xvalue)/sum(xvalue)
      
      lines(xvalue,cumul,col=color[which(colnames(mat)==i)])
    }
  } else stop("The specified sample / group does not exist.")
}


# Rank Scatter Plot
# evtl. umstellen auf Hexbin Plot
# 
# data: internal expression data object
# x: either specific sample name or group
# y: if x is a specific sample, y gives the name of a second sample to compare x to
# transformation: transformation to use; one of either raw matrix ("Matrix", the default), log, rlog, VST or CPM
rankScatter = function(data,x,y=NULL,transformation=c("log","rlog","vst","cpm")){
  transformation = match.arg(transformation)
  mat = data[[transformation]]
  x11()
  
  if(!is.null(y)){
		heatscatter(rank(mat[,x]),rank(mat[,y]),pch=20,xlab=x,ylab=y,main='')
    text(0.2*max(rank(mat[,x])),max(rank(mat[,y])),labels=paste0('r = ',signif(cor(mat[,x],mat[,y],method='spearman'),2)))
		abline(0,1,lty=2)
	} else{
	  if(!(x%in%data$Targets[[data$design]])) stop('The specified group does not exist.')
	  targets = data$Targets[data$Targets[[data$design]]==x,]
	  size = dim(targets)[1]
	  
	  par(mfrow=rep(size,2),mar=rep(2,4))
	  mat = mat[,colnames(mat)%in%rownames(targets)]
	  for(i in 1:size){
	    for(j in 1:size){
	      if(j<i){
	        plot(x=0,y=0,xaxt='n',yaxt='n',xlab='',ylab='',type='n',frame.plot=F)
	      } else if (i==j){
	        plot(x=0,y=0,xlim=c(0,1),ylim=c(0,1),xaxt='n',yaxt='n',xlab='',ylab='',type='n',frame.plot=F)
	        text(0.5,0.5,colnames(mat)[j],cex=1.5)
	      } else{
	        heatscatter(rank(mat[,i]),rank(mat[,j]),pch=20,main="",xlab="",ylab="")
	        text(0.3*max(rank(mat[,i])),0.9*max(rank(mat[,j])),labels=paste0('p = ',signif(cor(mat[,i],mat[,j],method='spearman'),2)))
	      }
	    }
		}
	}
}


# Correlation Heatmap
# 
# data: internal expression data object
# transformation: transformation to use; one of either raw matrix ("Matrix", the default), log, rlog, VST or CPM
# method: either Spearman rank correlation (default) or Pearson correlation
# groupBy: vector of colum names from the targets information file to group data by
corMap = function(data,transformation=c("Matrix","log","rlog","vst","CPM"),method=c('spearman','pearson'),groupBy=NULL){
	transformation = match.arg(transformation)
	method = match.arg(method)
	mat = data[[transformation]]
	
	x11()
	if(!is.null(groupBy)){
		df=data.frame(data$Targets[,groupBy])
		rownames(df) = rownames(data$Targets); colnames(df) = groupBy
		pheatmap(cor(mat,method=method),clustering_distance_rows='correlation',clustering_distance_cols='correlation',clustering_method='average',annotation_col=df,color=colorRampPalette(rev(c('red','yellow')))(100))
	} else 	pheatmap(cor(mat,method=method),clustering_distance_rows='correlation',clustering_distance_cols='correlation',clustering_method='average',color=colorRampPalette(rev(c('red','yellow')))(100))
}



# Principle Component Analysis
# 
# data: internal expression data object
# groupBy: name of the column of the targets file to be used for coloring the data points; by default, the design variable
# of the expression object (containing the group information) is used
# mode: either "2D" (first two principle components, default) or "3D" (first three principle components)
# transformation: transformation to use; one of either raw matrix ("Matrix", the default), log, rlog, VST or CPM
expressionPCA = function(data,groupBy="design",mode=c("2D","3D"),transformation=c("log","rlog","vst","cpm")){
  if(groupBy=="design") groupBy = data$design
  transformation = match.arg(transformation)
  mode = match.arg(mode)
  mat = data[[transformation]]
  col.code = as.numeric(as.factor(data$Targets[[groupBy]]))
  color = rainbow(max(col.code))
  
  pca = prcomp(mat,center=T)
  percentage = signif(pca$sdev^2/sum(pca$sdev^2)*100,2)
  x11()
  
  if(mode=="2D"){
    plot(pca$rotation[,1:2],pch=20,col=color[col.code],xlab=paste0("PC1 (",percentage[1]," % variance)"),ylab=paste0("PC2 (",percentage[2]," % variance)"),cex=2)
  } else{
    scatter3D(pca$rotation[,1],pca$rotation[,2],pca$rotation[,3],colvar=col.code,col=color,pch=20,colkey=F,xlab=paste0("PC1 (",percentage[1]," % variance)"),ylab=paste0("PC2 (",percentage[2]," % variance)"),zlab=paste0("PC3 (",percentage[3]," % variance)"),cex=2)#,theta=0,phi=90)
    legend("bottomleft",legend=as.factor(unique(data$Targets[[groupBy]])),fill=color)
  }
}

