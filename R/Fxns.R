library(RColorBrewer)
library(statmod)
library(sva)
library(devtools)

"%||%" <- devtools:::`%||%`

brewer16 = c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2"))
brewer16[6] = "khaki2"
brewer16[8] = "lightskyblue2"

cubehelix1.16 = c('#000000', '#1B0F00', '#411704', '#681B20', 
                  '#85214B', '#932D7E', '#9042AF', '#8160D2', '#6F83E3', 
                  '#63A6E2', '#65C5D3', '#78DBC2', '#99E9B9', '#C1F0BF', '#E6F5D8', '#FFFFFF')

### Compute the group-wise mean of a dataset.
group.means <- function(counts, groups, fn=mean, use.data.table=F)
{
  counts <- aggregate(t(counts), by=list(groups), FUN=fn)
  rownames(counts) = counts$Group.1
  counts$Group.1 = NULL
  r = t(counts)
  return(r)
}

# Logging utility function
info <- function(text, ...)
{
  cat(sprintf(paste(Sys.time(),"INFO:", text,"\n")))
}

# Logging utility function
warn <- function(text, ...)
{
  cat(sprintf(paste(Sys.time(),"WARN:", text,"\n")))
}

### Compute TPM expression values from raw UMI counts
tpm <- function(counts, mult=10000)
{
  info("Running TPM normalisation")
  total.counts = colSums(counts)
  scaled.counts = t(t(counts) / total.counts) 
  scaled.counts * mult
}
###
fread.zip <- function(zipfile,sep,header, quote) { 
  # Create a name for the dir where we'll unzip 
  zipdir <- tempfile() 
  # Create the dir using that name 
  dir.create(zipdir) 
  # Unzip the file into the dir 
  unzip(zipfile, exdir=zipdir) 
  # Get the files into the dir 
  files <- list.files(zipdir) 
  # Throw an error if there's more than one 
  if(length(files)>1) stop("More than one data file inside zip") 
  # Get the full name of the file 
  file <- paste(zipdir, files[1], sep="/") 
  # Read the file 
  fread(file, sep = sep,header = header, quote=quote,check.names=FALSE)
} 

###
download_FTP <- function(url,taskdir){
  require("RCurl")
  res <- c(status=0,path="")
  url <- gsub("https://","",url)
  url <- gsub("http://","",url)
  
  #1 判断是否存在
  public <- paste0(taskdir,"public/")
  
  dirlists <- getURL(url,verbose=F,ftp.use.epsv=F,dirlistonly = T)
  getfiles <- unlist(strsplit(dirlists,"\n",fixed = T))
  
  t1 <- sum(grepl("matrix.mtx.gz",getfiles))
  t2 <- sum(grepl("barcodes.tsv.gz",getfiles))
  t3 <- sum(grepl("genes.tsv.gz",getfiles))
  t4 <- sum(grepl("features.tsv.gz",getfiles))
  
  if(t1==1 && t2==1 && (t3==1 || t4 ==1)){
    mtx <- getfiles[grepl("matrix.mtx.gz",getfiles)]
    barcodes <- getfiles[grepl("barcodes.tsv.gz",getfiles)]
    if(t3 == 1){
      genes <- getfiles[grepl("genes.tsv.gz",getfiles)]
    }else if(t4 == 1){
      genes <- getfiles[grepl("features.tsv.gz",getfiles)]
    }
    GSM <- unlist(strsplit(mtx,"-"))[1]
    
    filepath <- paste0(public,GSM)
    
    if(!file.exists(filepath)){
      #not exist-download
      dir.create(file.path(filepath),recursive = TRUE)
      download.file(paste0(url,mtx),paste0(filepath,"/matrix.mtx.gz"))
      download.file(paste0(url,barcodes),paste0(filepath,"/barcodes.tsv.gz"))
      download.file(paste0(url,genes),paste0(filepath,"/genes.tsv.gz"))
    }
    
    if(file.exists(paste0(filepath,"/matrix.mtx.gz")) && 
       file.exists(paste0(filepath,"/barcodes.tsv.gz")) &&
       file.exists(paste0(filepath,"/genes.tsv.gz"))){
      res["status"] <- 3
      res["path"] <- filepath
      return(res)
    }else{
      res["status"] <- 2
      return(res)
    }
  }else{
    res["status"] <- 1
    return(res)
  }
  
}


### random select the subset
ransubset <- function(count,group = FALSE,maxcol = 10,maxrow = 10,subset_size = 0.1){
  if(group){
    
  }else{
    if(ncol(count) < maxcol) maxcol = ncol(count)
    if(nrow(count) < maxrow) maxrow = nrow(count)
    
    subset_idc <- sample.int(n = ncol(count), size = maxcol, replace=FALSE)
    subset_idr <- sample.int(n = nrow(count), size = maxrow, replace=FALSE)
    tmp <-  count[subset_idr,subset_idc]
  }
  return(tmp)
  
}

### as data frame
As.DF <- function(data){
  data <- as.data.frame(data)
  #data[is.na(data)] <- 0
  rownames(data) <- data[,1]
  data <- data[,-1]
  data <- as(as.matrix(data),"dgCMatrix")
  return(data)
}
### plot hvg
hvgPlotly <- function(object,cols = c("black", "red"),pt.size = 1){
  hvf.info <- HVFInfo(object = object)
  var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1]
  hvgs <- VariableFeatures(object)
  
  hvf.info <- hvf.info[, c(1, 3)]
  
  axis.labels <- switch(EXPR = colnames(x = hvf.info)[2], 
                        variance.standardized = c("Average Expression", "Standardized Variance"), 
                        dispersion.scaled = c("Average Expression", "Dispersion"), 
                        residual_variance = c("Geometric Mean of Expression", 
                                              "Residual Variance"))
  nvc <- rownames(hvf.info) %in% hvgs
  labels = paste(c("Non-variable","Variable"), "count:", table(nvc))
  nvc[nvc] <- labels[2]
  nvc[nvc==FALSE] <- labels[1]
  hvf.info$hvgs <- nvc
  
  top10 <- hvf.info[rownames(hvf.info) %in% head(hvgs,10),]
  
  a <- list(
    x = log10(top10$mean),
    y = top10$variance.standardized,
    text = rownames(top10),
    xref = "x",
    yref = "y",
    showarrow = TRUE,
    arrowhead = 4,
    arrowsize = .3,
    ax = 15,
    ay = -20
  )
  
  py <- plot_ly(hvf.info, x = ~mean, y = ~variance.standardized,color = ~hvgs,colors = cols,
                text= ~paste("Genes: ", rownames(hvf.info))) %>% add_markers() %>%
    layout(xaxis = list(type = "log",title = axis.labels[1]),
           yaxis = list(title = axis.labels[2]),annotations=a) %>% config(displaylogo = FALSE)
    
  return(py)
}

DimPlotly <- function(object,reduction,label=TRUE,pt.size = 0.5,label.size = 14){
 
  dt <- as.data.frame(Embeddings(object,reduction=reduction))
  a <- NULL
  if(label){
  data <- as.data.table(dt)
  data$group <- as.character(Idents(object))
  labels <- data[,lapply(.SD, mean, na.rm=TRUE),by=group]
  a <- list(
    x = as.numeric(unlist(labels[,2])),
    y = as.numeric(unlist(labels[,3])),
    text = labels$group,
    xref = "x",
    yref = "y",
    showarrow = F,
    font = list(family = 'sans serif',
                size = label.size,weight='bolder')
  )
  }
  
  p <- plot_ly(dt, x = dt[,1], y = dt[,2],color = Idents(object),colors = "Set1",
               text= ~paste(rownames(dt)),marker = list(size = pt.size)) %>% add_markers() %>%
    layout(xaxis = list(title = colnames(dt)[1]),
           yaxis = list(title = colnames(dt)[2]),annotations=a) %>% config(displaylogo = FALSE)
  
  return(p)
}

### Run ComBat batch correction from the SVA package
batch.normalise.comBat <- function(data, batch.groups=NULL){
  combat_data <- logcounts(data)
  mod_data <- as.data.frame(t(combat_data))
  
  batch.groups <- batch.groups %||% data$orig.ident
  
  correct.data <- ComBat(
    dat = t(mod_data), 
    batch = factor(batch.groups), 
    par.prior = TRUE,
    prior.plots = FALSE
  )
  return(correct.data)
}


batch.normalise.ruvg <- function(data, batch.groups=NULL, max.val=6,k=10){
  require(RUVSeq)
  matrix.data <- as.matrix(counts(data))
  ruvg <- RUVg(matrix.data, erccs, k)
  correct.data <- log2(
    t(t(ruvg$normalizedCounts) / colSums(ruvg$normalizedCounts) * 1e6) + 1
  )
  return(correct.data)
}

glm_fun <- function(g, batch, indi) {
  model <- glm(g ~ batch + indi)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}

batch.normalise.glm <- function(data,batch,indi){
  effects <- apply(
    logcounts(data), 
    1, 
    glm_fun, 
    batch = batch, 
    indi = indi
  )
  correct.data <- logcounts(data) - t(effects[as.numeric(factor(batch)), ])
  return(correct.data)
}


glm_fun1 <- function(g, batch) {
  model <- glm(g ~ batch)
  model$coef[1] <- 0 # replace intercept with 0 to preserve reference batch.
  return(model$coef)
}

do_glm <- function(data.qc) {
  effects <- apply(
    logcounts(data.qc), 
    1, 
    glm_fun1, 
    batch = data.qc$orig.ident
  )
  corrected <- logcounts(data.qc) - t(effects[as.numeric(factor(data.qc$orig.ident)), ])
  return(corrected)
}

batch.normalise.glm_indi <- function(data,batch,indi){
  labels <- names(summary(indi))
  correct.data<- data.frame()
  for (order in 1:length(labels)) {
    
    label <- labels[order]
    do_indi <- do_glm(data[, indi == label])
    correct.data <- cbind(correct.data,do_indi)
  }
  return(correct.data)
}
#generate the filter info dataframe
filter.info <- function(id,dataset.name,species,tech,scRNA.data,scRNA,is.spike){
  Raw.cells <- ncol(scRNA.data)
  Raw.genes <- nrow(scRNA.data)
  Raw.sparsity <- round(sum(scRNA.data == 0)/(Raw.cells*Raw.genes),3)*100
  Final.cells <- ncol(scRNA@assays$RNA@scale.data)
  Final.genes <- nrow(scRNA@assays$RNA@scale.data)
  Final.sparsity <- round(sum(scRNA@assays$RNA@scale.data == 0)/(Final.cells*Final.genes),3)*100
  spike.num <- sum(is.spike)
  result <- data.frame(ID=id,Dataset.name=dataset.name,Species=species,
                       Protocol =tech,Raw.cells=Raw.cells,Raw.genes=Raw.genes,
                       Raw.sparsity=Raw.sparsity,Final.cells=Final.cells,
                       Final.genes=Final.genes,Final.sparsity=Final.sparsity,
                       is.spike=spike.num)
  return(result)
}

plot_umap <- function(reduce_var,groups,title=NULL,colorbar_name=NULL){
  
  n <- length(table(groups))
  if(n>30){
    colors = NULL
    
  }else{
    colors = hue_pal()(n)
    groups <- as.factor(groups)
  }
  cols <- colnames(reduce_var)
  fig <- plot_ly(x = reduce_var[,1], y = reduce_var[,2], color =  groups, colors = colors,text=rownames(reduce_var),marker = list(size = 2.5))
  
  fig <- fig %>% add_markers() %>% colorbar(title = colorbar_name)
  fig <- fig %>% layout(title=title,
                        xaxis = list(title = cols[1],zeroline = F,showgrid=F),
                        yaxis = list(title = cols[2],zeroline = F,showgrid=F)
                        #legend = list(x = 0.5,y=-0.2, orientation="h",xanchor="center")
                        ) %>% config(displaylogo = FALSE)
  
}


plotly_3d <- function(reduce_var,groups,title,colorbar_name=NULL){
  
  n <- length(table(groups))
  if(n>30){
    colors = NULL
    
  }else{
    colors = hue_pal()(n)
    groups <- as.factor(groups)
  }
  cols <- colnames(reduce_var)
  fig <- plot_ly(x = reduce_var[,1], y = reduce_var[,2], z = reduce_var[,3], color =  groups, colors = colors,size=1,text=rownames(reduce_var),marker = list(size = 2.5))
  fig <- fig %>% add_markers() %>%  colorbar(title = colorbar_name)
  
  fig <- fig %>% layout(title=title,
                        scene = list(xaxis = list(title = cols[1],zeroline = T,showgrid=T),
                                     yaxis = list(title = cols[2],zeroline = T,showgrid=T),
                                     zaxis = list(title = cols[3],zeroline = T,showgrid=T))) %>% config(displaylogo = FALSE)
  
}

estimate.kBET <- function(data,batch){
  subset_size <- 0.2 #subsample to 10% of the data
  subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
  batch.estimate <- kBET(data[subset_id,], batch[subset_id],plot = FALSE)
  #batch.estimate <- kBET(data, batch,plot = FALSE)
  
  #batch.estimate$stats$kBET.expected <- 1- batch.estimate$stats$kBET.expected
  #batch.estimate$stats$kBET.observed <- 1- batch.estimate$stats$kBET.observed
  return(batch.estimate)
}

estimate.ASW <- function(data,batch,max_val=1,min_val=-1){
  library(cluster)
  dis <- dist(data)
  #dis <- dist(data,method='euclidean')
  
  x <- as.numeric(as.factor(batch))
  sil <- silhouette(x, dis)
  asw_batch_norm <- ((sil[,3] - min_val)/(max_val - min_val))
  return(asw_batch_norm)
}

estimate.ARI <- function(data,batch){
  library(mclust)
  #library(NbClust)
  
  subset_size <- 0.8 #subsample to 10% of the data
  nbiters = 20
  ari_avg =0
  for(i in 1:nbiters) {
  subset_id <- sample.int(n = nrow(batch), size = floor(subset_size * nrow(batch)), replace=FALSE)
  sub_data <- data[subset_id,]
  sub_batch <- batch[subset_id,]
  clustering_result <- kmeans(x = sub_data, centers=8, iter.max = 30)
  sub_batch$clusterlb <- clustering_result$cluster
  
  ari_batch <- mclust::adjustedRandIndex(sub_batch$project.PIAS,sub_batch$clusterlb)
  ari_avg <- ari_batch +ari_avg
  }
  ari_avg <-(ari_avg/20)
  
  return(ari_avg)
}

estimate.CH <- function(data,batch){
  library(fpc)
  res <- fpc::calinhara(as.data.frame(data),as.numeric(as.factor(batch$project.PIAS)))
  
  res <-  round(res,digits=2)
  
}


## EWM
entropy <- function(d){
  res <- 0
  for(i in 1:length(d))
  {
    if(d[i]!=0)
      res <- res + d[i]*log(d[i])
  }
  return (-res/log(length(d)))
}
negative.scale <-function(x){
  if(max(x)==0) return(x)
  x=(max(x)-x)/(max(x)-min(x))
  return(round(x,2))
}

positive.scale<-function(x){
  if(max(x)==0) return(x)
  x=(x-min(x))/(max(x)-min(x))
  return(round(x,2))
}

Percentage <-function(x){
  if(sum(x)==0) return(x)
  x=x/sum(x)
  return(x)
}
eweight <- function(e){
  d <- 1-e
  d[d==1] <- 0
  w<-round(d/sum(d),2)
}

##
toMonocle <- function(otherCDS, import_all = FALSE) {
  if (class(otherCDS)[1] == "Seurat") {
    requireNamespace("Seurat")
    if("integrated" %in% names(otherCDS@assays)){
      data <- as.data.frame(otherCDS@assays$integrated@data)
      
    }else{
      data <- as.data.frame(otherCDS@assays$RNA@data)
      
    }
    data <- data[rownames(otherCDS@meta.data)]
    
    data <- as(as.matrix(data), "sparseMatrix")
    
    pd <- tryCatch({
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, error = function(e) {
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      message("This Seurat object doesn't provide any meta data")
      pd
    })
    if (length(setdiff(colnames(data), rownames(pd))) > 
        0) {
      data <- data[, rownames(pd)]
    }
    
    fData <- data.frame(gene_short_name = row.names(data), 
                        row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0
    if (all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    }else if (any(data < 0)) {
      expressionFamily <- uninormal()
    }else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    if(!identical(colnames(data),rownames(pd))){
      stop("Error:names differ between data and phenoData")
    }else if(!identical(rownames(data),rownames(fd))){
      stop("Error:names differ between data and featureData")
    }
    monocle_cds <- newCellDataSet(data, phenoData = pd, 
                                  featureData = fd, lowerDetectionLimit = lowerDetectionLimit, 
                                  expressionFamily = expressionFamily)
    if (import_all) {
      if ("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      }else {
        mist_list <- otherCDS
      }
    }else {
      mist_list <- list()
    }
    if ("var.features" %in% slotNames(otherCDS@assays[["integrated"]])) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays[["integrated"]]@var.features)
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
  }else {
    stop("the object type you want to export to is not supported yet")
  }
  return(monocle_cds)
}

### update homologenes for seurat

UpdateHomoSeurat <- function(obj,inTax,outTax){
  library(homologene)
  
  genelist <- rownames(obj)
  
  molos <- homologene(genelist, inTax = inTax, outTax = outTax)
  
  #molos[duplicated(molos[,2]),2] <- paste0(molos[duplicated(molos[,2]),1],"+",molos[duplicated(molos[,2]),2])
  
  selgenes <- genelist[match(molos[,1],genelist)]
  
  obj <- subset(obj,features = unique(selgenes))
  
  genelist <- rownames(obj)
  newgene <- molos[match(genelist,molos[,1]),2]
  
  #genelist[which(duplicated(genelist))] <- paste0(genelist[which(duplicated(genelist))],"_2")
  #dups <- unique(molos[duplicated(molos[match(molos[,1],genelist),1]),1])
  
  rownames(obj@assays$RNA@counts) <- newgene
  rownames(obj@assays$RNA@data) <- newgene
  rownames(obj@assays$RNA@scale.data) <- newgene
  
  vars <- homologene(obj@assays$RNA@var.features, inTax = inTax, outTax = outTax)
  obj@assays$RNA@var.features <- vars[,2]
  return(obj)
}

### Automatically find a single cell type
AutoSingleType <- function(difMarker,species,tissue,auto){
  
  if(auto){
    
    newmark <- data.frame(cluster = difMarker$cluster,gene = difMarker$gene,avg_logFC = difMarker$avg_logFC)
    
    i=0.1
    pinmarker <- data.frame()
    pincluster <- c()
    while (i<1) {
      cur_newmark <- newmark %>% group_by(cluster) %>% top_n(n() * i,wt= avg_logFC)
      if(!is.null(pincluster)){
        cur_newmark <- subset(cur_newmark,!cluster %in% pincluster)
        cur_newmark <- rbind(cur_newmark,pinmarker)
      }
      clu_annotation <- scCATCH(object = cur_newmark[,1:2],species = species,tissue = tissue)
      
      is.confusion <-  is.na(clu_annotation$cell_type) | grepl(",",clu_annotation$cell_type)
      
      pincluster <- clu_annotation[!is.confusion,]$cluster
      
      pinmarker <- cur_newmark[cur_newmark$cluster %in% pincluster,]
      if(sum(!levels(newmark$cluster) %in% pincluster) == 0) break
      i=i+0.1
    }
    
    order_annotation <- clu_annotation[order(clu_annotation$cluster),]
    
    ffs <- data.frame()
    
    for (j in 1:length(order_annotation$cluster)) {
      tmp <- subset(cur_newmark,cluster == order_annotation$cluster[j])
      if(j==1) ffs = tmp
      else ffs <- rbind(ffs,tmp)
    }
    clu_annotation[is.na(clu_annotation)] <- "unknown"
    res <- list(annotation=clu_annotation,subMarkers=ffs)
    
  }else{
    gene.ids <- difMarker$gene
    gene.ids[is.na(gene.ids)] <- names(gene.ids)[is.na(gene.ids)]
    newmark <- data.frame(cluster = difMarker$cluster,gene = gene.ids)
    clu_annotation <- scCATCH(object = newmark,species = species,tissue = tissue)
    clu_annotation[is.na(clu_annotation)] <- "unknown"
    res <- list(annotation=clu_annotation,subMarkers=difMarker)
    
  }
  
  return(res)
}


### Get variable genes. Code adapted from:
### | Brennecke et al, Accounting for technical noise in single-cell RNA-seq experiments
### | Nature Methods 10, 1093???1095 (2013), doi:10.1038/nmeth.2645
### 	See: https://images.nature.com/original/nature-assets/nmeth/journal/v10/n11/extref/nmeth.2645-S2.pdf
### 	and: http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html
get.variable.genes <- function(ed, min.cv2=2, pdf=NULL, width=9, height=8, do.plot=T, p.thresh=0.05)
{
  
  means <- rowMeans(ed)
  vars <- apply(ed,1,var)
  cv2 <- vars/means^2
  minMeanForFit <- unname( quantile( means[ which( cv2 > min.cv2 ) ], .95 ) )
  useForFit <- means >= minMeanForFit # & spikeins
  info(sprintf("Fitting only the %s genes with mean expression > %s", sum(useForFit), minMeanForFit))
  fit <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ), cv2[useForFit] )
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])
  if(do.plot){par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9); smoothScatter(log(means),log(cv2))}
  xg <- exp(seq( min(log(means[means>0])), max(log(means)), length.out=1000 ))
  vfit <- a1/xg + a0
  if(do.plot){lines( log(xg), log(vfit), col="black", lwd=3 )}
  
  df <- ncol(ed) - 1
  # add confidence interval
  if(do.plot){
    lines(log(xg),log(vfit * qchisq(0.975,df)/df),lty=2,col="black")
    lines(log(xg),log(vfit * qchisq(0.025,df)/df),lty=2,col="black")
  }
  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  varorder <- order(varFitRatio, decreasing=T)
  oed <- ed[varorder,]
  pval <- pchisq(varFitRatio*df,df=df,lower.tail=F)
  adj.pval <- p.adjust(pval,"fdr")
  r = data.frame(rownames(ed), varFitRatio, pval, adj.pval)
  colnames(r) = c("Gene", "VarianceFitRatio", "p", "p.adj")
  v = r[!is.na(r$p.adj),]
  n.sig = sum(v$p.adj<p.thresh)
  info(sprintf("Found %s variable genes (p<0.05)", n.sig))
  
  # add top 100 genes
  if(do.plot){
    points(log(means[varorder[1:n.sig]]),log(cv2[varorder[1:n.sig]]),col=2)
  }
  r = r[order(r$VarianceFitRatio, decreasing=T), ]
  r$Rank = 1:nrow(r)
  return(r)
}


# Test for significant PCs adapted from: 
# 
#	' Permutation Parallel Analysis
#	'
#	' Estimate a number of significant principal components from a permutation test
#   B is the number of permutations
#   threshold is p-value for significance
#'
sig.pcs.perm <- function (dat, B = 100, threshold = 0.05,
                          randomized=F,  
                          verbose=TRUE, seed = NULL,
                          max.pc=100, n.cores=1, 
                          center=T, scale=T) {
  ptm <- proc.time()
  if(B %% n.cores != 0){stop("Permutations must be an integer multiple of n.cores")}
  cat(sprintf("Scaling input matrix [center=%s, scale=%s]\n", center, scale))
  dat = t(dat)
  dat = as.matrix(t(scale(t(dat), center=center, scale=scale)))
  if (!is.null(seed)) set.seed(seed)
  n <- min(max.pc, ncol(dat))
  m <- nrow(dat)
  print(paste0("Considering only the top ", n, " PCs. Supply max.pc if you wish to change"))
  cat(sprintf("Running initial PCA\n"))
  if(randomized){
    library(rsvd)
    uu <- rsvd(as.matrix(dat), k=max.pc)
  }else{
    uu <- corpcor::fast.svd(dat, tol = 0)
  }
  ndf <- n - 1
  dstat <- uu$d[1:ndf]^2/sum(uu$d[1:ndf]^2)
  dstat0 <- matrix(0, nrow = B, ncol = ndf)
  if(verbose==TRUE) message("Estimating number of significant principal components. Permutation: ")
  #permutations
  if(n.cores==1){
    for (i in 1:B) {
      if(verbose==TRUE) cat(paste(i," "))
      dat0 <- t(apply(dat, 1, sample, replace = FALSE))
      if(randomized){
        library(rsvd)
        uu0 <- rsvd(as.matrix(dat0), k=max.pc)
      }else{
        uu0 <- corpcor::fast.svd(dat0, tol = 0)
      }
      dstat0[i, ] <- uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
    }
  }else{
    library(parallel)
    library(foreach)
    library(doParallel)
    cl<-makePSOCKcluster(n.cores, outfile="")
    registerDoParallel(cl, n.cores)
    chunksize = B/n.cores
    vals = split(1:B, ceiling(seq_along(1:B)/chunksize))
    dstat0 = foreach(run.id=1:n.cores, .packages="corpcor", .combine=cbind) %dopar% {
      v = vals[[run.id]]
      #cat(sprintf("Core %s will run perms: %s \n", run.id, paste(v, collapse=",")))
      do.call(rbind, lapply(v, function(i) {
        if(verbose==TRUE) cat(paste(i," "))
        dat0 <- t(apply(dat, 1, sample, replace = FALSE))
        
        if(randomized){
          library(rsvd)
          uu0 <- rsvd(as.matrix(dat0), k=max.pc)
        }else{
          uu0 <- corpcor::fast.svd(dat0, tol = 0)
        }
        uu0$d[1:ndf]^2/sum(uu0$d[1:ndf]^2)
      }))
    }
    cat("\nUnregistering parallel backend..")
    stopCluster(cl)
    registerDoSEQ()
    cat(" done\n");
  }
  p <- rep(1, n)
  for (i in 1:ndf) {
    p[i] <- mean(dstat0[, i] >= dstat[i])
  }
  for (i in 2:ndf) {
    p[i] <- max(p[(i - 1)], p[i])
  }
  r <- sum(p <= threshold)
  y = proc.time() - ptm
  cat(sprintf("\n\n PC permutation test completed. \n %s PCS significant (p<%s, %s bootstraps)\n Runtime: %s s\n ", r,  threshold, B,signif(y[["elapsed"]], 3)))
  return(list(r = r, p = p))
}

build_knn_graph <- function(dm, k=200, verbose=F)
{
  if(k==0)
  {
    k = floor(sqrt(nrow(dm))/2)
  }
  if(verbose)
  {
    info(sprintf("Building %s-nearest [%s] neighbor graph..", k, dist.type))
  }
  g <- nng(dx=dm,k=k)
  V(g)$name = rownames(dm)
  if(verbose)
  {
    info(sprintf("%s %s-NN computed. Average degree: %s", dist.type, k, mean(degree(g))))
  }
  return(g)
}


CheckCellNames <- function (object.list, verbose = TRUE, stop = FALSE) 
{
  cell.names <- unlist(x = lapply(X = object.list, FUN = colnames))
  if (any(duplicated(x = cell.names))) {
    if (stop) {
      stop("Duplicate cell names present across objects provided.")
    }
    if (verbose) {
      warning("Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.")
    }
    object.list <- lapply(X = 1:length(x = object.list), 
                          FUN = function(x) {
                            return(RenameCells(object = object.list[[x]], 
                                               new.names = paste0(Cells(x = object.list[[x]]), 
                                                                  "_", x)))
                          })
  }
  return(object.list)
}


# graph.type can be jaccard, invlogweighted or dice, community detect
# can be louvain, infomap or markov. 
cluster_graph <- function(	g, 
                           graph.type="knn", # can be threshold (binarise the distance matrix), jaccard or knn.
                           dm=NULL,
                           community.detect="infomap", 
                           distance.method="euclidean",
                           k=0)
{
  if(identical(toupper(community.detect), toupper("markov")))
  {
    r = igraph::cluster.markov(g)
    clusters = r$Cluster
  }else{
    if(identical(toupper(community.detect), toupper("louvain")))
    {
      r = igraph::multilevel.community(as.undirected(g))
      clusters = r$membership
    }else{
      if(identical(toupper(community.detect), toupper("infomap")))
      {
        r = igraph::infomap.community(g, modularity=TRUE)
        clusters = r$membership
      }else{
        error(sprintf("Unknown community detection method: %s", community.detect))
        return (FALSE)
      }
    }
  }
  n.clusters =length(unique(clusters))
  f = function(i){as.vector(clusters==i)}
  clist= lapply(1:n.clusters, f)
  m = igraph::modularity(g, clusters)
  return (list("result"=r,
               "clustermethod"=paste(graph.type, "-graph clustering [", community.detect,"]", sep=""), 
               "nc"=n.clusters, 
               "modularity"=m, 
               "clusterlist"=clist,		
               "partition"=clusters))
}


merge_clusters <- function(clustering, clusters.to.merge, new.name=NULL)
{
  if(length(clustering) < 2){cat("ERROR: Must provide 2 or more cluster ID's to merge!");return (clustering)}
  i = 1
  if(!is.null(new.name)){
    use.id = new.name
    levels(clustering) = c(levels(clustering), use.id)
    clustering[which(clustering == clusters.to.merge[1])] = use.id
  }else
  {use.id = clusters.to.merge[1]}
  for(id in clusters.to.merge)
  {
    if(i > 1)
    {
      cat(sprintf("Merging cluster %s into %s ..\n", id, use.id))
      clustering[which(clustering == id)] = use.id
    }
    i = i + 1 
  } 
  factor(clustering)
}


checklogin <- function(memid, userId)
{
  con <- dbConnect(RMariaDB::MariaDB(), username = "root", password = "bmilab",dbname = "pscap",host="127.0.0.1")
  
  id <- memid
  kk <- userId
  sql1 <- paste0("SELECT userID,username FROM members WHERE memberID ='", id,"';")
  
  res <- dbSendQuery(con,sql1)
  
  result <- dbFetch(res)
  if(nrow(result) == 1){
    nuserID <- result$userID
    if(kk == md5(paste0(nuserID,"pscap")))
    {
      ruall <- c(nuserID,result$username)
      return(ruall)
    }else{
      return(FALSE)
    }
    
  }else{
    return(FALSE)
  }
  
  dbClearResult(res)
  
  # Disconnect from the database
  dbDisconnect(con)
  
}

FoldNames <- function(ids,obj){
  names(ids) <- levels(obj@active.ident)
  obj <- RenameIdents(object = obj, ids)
  return(as.character(obj@active.ident))
}

LoadOrgDB <- function(species){
  if(species == "Homo sapiens")
  {
    library(org.Hs.eg.db)
    Org.Db <- org.Hs.eg.db
  }else if(species == "Mus musculus"){
    library(org.Mm.eg.db)
    Org.Db <- org.Mm.eg.db
  }else if(species == "Arabidopsis thaliana"){
    library(org.At.tair.db)
    Org.Db <- org.At.tair.db
  }else if(species == "Xenopus (Silurana) tropicalis"){
    library(org.Xl.eg.db)
    Org.Db <- org.Xl.eg.db
  }else if(species == "Rattus norvegicus"){
    library(org.Rn.eg.db)
    Org.Db <- org.Rn.eg.db
  }else if(species == "Saccharomyces cerevisiae"){
    library(org.Sc.sgd.db)
    Org.Db <- org.Sc.sgd.db
  }else if(species == "Anopheles gambiae"){
    library(org.Ag.eg.db)
    Org.Db <- org.Ag.eg.db
  }else if(species == "Gallus gallus"){
    library(org.Gg.eg.db)
    Org.Db <- org.Gg.eg.db
  }else if(species == "Pan troglodytes"){
    library(org.Pt.eg.db)
    Org.Db <- org.Pt.eg.db
  }else if(species == "Bos taurus"){
    library(org.Bt.eg.db)
    Org.Db <- org.Bt.eg.db
  }else if(species == "Danio rerio"){
    library(org.Dr.eg.db)
    Org.Db <- org.Dr.eg.db
  }
  return(Org.Db)
}


ReOrder <- function(object.list,features){
  for (i in 1:length(object.list)) {
    object.list[[i]] <- object.list[[i]][features,]
    object.list[[i]]@assays$RNA@counts <- object.list[[i]]@assays$RNA@counts[features,]
    object.list[[i]]@assays$RNA@data <- object.list[[i]]@assays$RNA@data[features,]
    object.list[[i]]@assays$RNA@scale.data <- object.list[[i]]@assays$RNA@scale.data[features,]
    object.list[[i]]@assays$RNA@meta.features <- object.list[[i]]@assays$RNA@meta.features[features,]
  }
  return(object.list)
}
