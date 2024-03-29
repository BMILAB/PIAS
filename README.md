# PIAS
a web-based interactive platform for integrative analysis of multi-source single-cell RNA-seq datasets
PIAS is a web-based interactive platform for integrative analysis of multi-source single-cell RNA-seq datasets. Different from many other single-cell RNA-seq analysis 
platforms or pipelines that mainly focus on preprocessing or analysis of one single-cell RNA-seq dataset, PIAS has the unique feature of integrating multi-source datasets 
and incorporates various metrics for comprehensively evaluating the result of data integration. Moreover, PIAS provides rich functions for data preprocessing, comprehensive 
analyses and visualization, including gene name transfer, quality control, normalization, highly variable genes identification, batch-effect removal, dimensionality 
reduction, clustering, differentially expressed, cluster annotation, enrichment analysis, and single-cell trajectories construction. Users can freely choose to perform 
desired functions, visualize results, and transfer data through interactive operations with PIAS. As an easy-to-use interactive website, PIAS is an extendable platform for 
preprocessing and integrative analysis of multi-source single-cell RNA-seq datasets.

## Getting started

The schematic of PIAS workflow is illustrated in the figure. The entire workflow is divided into three major parts: data preprocessing, integration and downstream analysis. 
PIAS allows users to upload multiple gene expression matrices and corresponding metadata to perform data integration, and embeds rich functions for comprehensive single-cell 
analysis. PIAS is highly modular and project-oriented; results from each step can be visually displayed and stored in a project.

<img src=workflow.png height="1000">


## Installation

There are three different ways to launch PIAS: 
1)	Option 1: direct access to PIAS online. 
URL: http://bigbio.xmu.edu.cn/PIAS/
2)	Option 2: running via PIAS software package 
Download package from URL: http://bigbio.xmu.edu.cn/PIAS/download/PIAS.zip.  
Enter the R console or Rstudio, and then use runApp()
```{r}
runApp(dirpath/PIAS)
```
3)	Option 3: run app directly through github
```{r}
install.packages("devtools")
devtools::install_github("rstudio/shiny", dependencies=FALSE)
shiny::runUrl("https://github.com/BMILAB/PIAS/archive/master.zip")
```

## Public data
When running the local version of PIAS, you need to move the public data to the www/task/public folder, where you can download the public data set: https://doi.org/10.6084/m9.figshare.13205924. The list after the move is as follows:
<img src=public.jpg>

## File format
__*It is important that the input file should follow the same format as descibed below.*__

1) Normal gene expression matrix and metadata file, such as: [GSE81076_CEL_seq.csv](www/example_data/GSE81076_CEL_seq.zip), [GSE81076_metadata.csv](www/example_data/GSE81076_metadata.zip).
   
2) Sparse matrix and metadata file (similar to [10x Genomics Single Cell Gene Expression Datasets](https://www.10xgenomics.com/resources/datasets/)), such as:  [matrix.mtx(.gz)](http://bigbio.xmu.edu.cn/PIAS/download/matrix.mtx.gz), [barcodes.tsv(.gz)](http://bigbio.xmu.edu.cn/PIAS/download/barcodes.tsv.gz), [features.tsv(.gz)](http://bigbio.xmu.edu.cn/PIAS/download/features.tsv.gz).

3) Cell markers file,such as: [Pancreas_markers.txt](www/example_data/Pancreas_markers.txt).


## Load projects:
We have uploaded the three public projects in PIAS web server to [Figshare](https://doi.org/10.6084/m9.figshare.13205924). Users can find our database and download it by searching title. PIAS web server saves all public projects. PIAS local version comes with the GSE84133 public project. If you want to use these items in the local version of PIAS, you need to download the item, unzip the zip file and move the RData file to www/task/public to load the corresponding project in PIAS local version.

