#!/usr/bin/Rscript

# Check and install if the package is not installed and then load them into the R session.

checkPackage <- function(pkgs,repo){
  pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(pkgs)){
    if(repo == 'CRAN'){
      install.packages(pkgs, dependencies = TRUE, repos='https://cloud.r-project.org/')
    }
    if(repo == "BiocManager"){
      BiocManager::install(pkgs,dependencies = T)
    }
  }
  sapply(pkgs, library, character.only = TRUE)
}

CRANpkgs <- c("BiocManager","Matrix", "Seurat","shiny","shinyjs","DT","dplyr","data.table","ggplot2",
              "plotly","shinyalert","glue","shinyBS","DBI","RMariaDB","openssl","patchwork","scales",
              "stringi","reactable","shinydashboard","factoextra","cluster","devtools","RColorBrewer","statmod",
              "mclust","homologene","htmltools","fpc","venn")
checkPackage(CRANpkgs, "CRAN")

Biopkgs <- c("gage", "clusterProfiler","SingleCellExperiment","scater","cowplot","batchelor",
              "org.Hs.eg.db","org.Mm.eg.db","pathview",'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats','limma',
              'S4Vectors','SummarizedExperiment', 'Matrix.utils',"sva","AUCell","enrichplot")

checkPackage(Biopkgs, "BiocManager")

# github

pkg <- "SeuratWrappers"
if(!pkg %in% installed.packages()[, "Package"]){
  remotes::install_github('satijalab/seurat-wrappers',dependencies=FALSE)
}
require(pkg,character.only = TRUE)

pkg <- "scCATCH"
if(!pkg %in% installed.packages()[, "Package"]){
  devtools::install_github('ZJUFanLab/scCATCH',dependencies=FALSE)
}
require(pkg,character.only = TRUE)

pkg <- "harmony"
if(!pkg %in% installed.packages()[, "Package"]){
  devtools::install_github("immunogenomics/harmony",dependencies=FALSE)

}
require(pkg,character.only = TRUE)

pkg <- "liger"
if(!pkg %in% installed.packages()[, "Package"]){
  devtools::install_github('MacoskoLab/liger',dependencies=FALSE)
  
}
require(pkg,character.only = TRUE)


pkg <- "kBET"
if(!pkg %in% installed.packages()[, "Package"]){
  devtools::install_github('theislab/kBET',dependencies=FALSE)
}
require(pkg,character.only = TRUE)

pkg <- "monocle3"
if(!pkg %in% installed.packages()[, "Package"]){
  devtools::install_github('cole-trapnell-lab/leidenbase')
  devtools::install_github('cole-trapnell-lab/monocle3',dependencies=FALSE)
}
require(pkg,character.only = TRUE)

pkg <- "scRNA.seq.funcs"
if(!pkg %in% installed.packages()[, "Package"]){
  devtools::install_github("hemberg-lab/scRNA.seq.funcs",dependencies=FALSE)
}
require(pkg,character.only = TRUE)

