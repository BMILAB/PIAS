# PIAS
a web-based interactive platform for integrative analysis of multi-source single-cell RNA-seq datasets
PIAS is a web-based interactive platform for integrative analysis of multi-source single-cell RNA-seq datasets. Different from many other single-cell RNA-seq analysis 
platforms or pipe-lines that mainly focus on preprocessing or analysis of one single-cell RNA-seq dataset, PIAS has the unique feature of integrating multi-source datasets 
and incorporates various metrics for comprehensive-ly evaluating the result of data integration. Moreover, PIAS provides rich functions for data prepro-cessing, comprehensive 
analyses and visualization, including gene name transfer, quality control, nor-malization, highly variable genes identification, batch-effect removal, dimensionality 
reduction, cluster-ing, differentially expressed, cluster annotation, enrichment analysis, and single-cell trajectories con-struction. Users can freely choose to perform 
desired functions, visualize results, and transfer data through interactive operations with PIAS. As an easy-to-use interactive website, PIAS is an extendable platform for 
preprocessing and integrative analysis of multi-source single-cell RNA-seq datasets.

## Getting started

The schematic of PIAS workflow is illustrated in the figure. The entire workflow is divided into three major parts: data preprocessing, integration and downstream analysis. 
PIAS allows users to upload multiple gene expression matrices and corresponding metadata to perform data integration, and embeds rich functions for comprehensive single-cell 
analysis. PIAS is highly modular and project-oriented; results from each step can be visually displayed and stored in a project.

<img src=images/iS-CellR_workflow.png height="800">


## Installation

There are three different ways to launch PIAS:
1)	Option 1: direct access to PIAS online
URL: http://bigbio.xmu.edu.cn/PIAS/
2)	Option 2: running via PIAS software package 
Download package from URL: http://bigbio.xmu.edu.cn/PIAS/download/PIAS.zip
Enter the R console or Rstudio, and then use runApp()
```{r}
runApp(dirpath/PIAS)
```
3)	Option 3: run app directly through github
```{r}
install.packages("devtools")
devtools::install_github("rstudio/shiny", dependencies=FALSE)
shiny::runUrl(“https://github.com/gangeszs/PIAS/archive/master.zip”)
```

## File format
__*It is important that the input file should follow the same format as descibed below.*__

1) 


## Demo data files:

There is one matrix file for Malignant dataset from Tirosh et al., 2016 (Main file for the analysis) and genes file (For STEP: Discriminating marker genes). 

### Count matrix: Maligant50.csv
### Genes file: Genes_MITF_AXL.txt
