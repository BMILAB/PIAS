source(R/install_package.R)

library(Matrix)
library(Seurat)
library(SeuratWrappers)
library(shiny)
library(shinyjs)
library(DT)
library(gage)
library(dplyr)
library(data.table)
#library(cytofkit)
#library(RColorBrewer)
library(clusterProfiler)
library(ggplot2)
#library(monocle)
library(plotly)
#library(Cairo)
library(shinyalert)
library(glue)
library(shinyBS)
library(DBI)
#library(RMariaDB)
library(openssl)
#library(enrichplot)
library(SingleCellExperiment)
library(scater)
#library(shinyWidgets)
library(cowplot)
library(patchwork)
library(batchelor) #mnnCorrect
library(scales)
library(stringi)
#devtools::install_github("immunogenomics/lisi")
library(shinydashboard)
library(reactable)
library(RCurl)
library(shinyBS)
library(htmltools)
options(shiny.maxRequestSize=5000*1024^2,shiny.sanitize.errors = FALSE)
#memory.limit(10000000)
#--load function ---
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = '*.R')) {
    message(nm)
    source(file.path(path, nm), ...)
  }
}

sourceDir(path = 'R',encoding = 'UTF-8')

taskdir <- paste0(getwd(),"/www/task/")
uuiddir <- "/var/www/html/PIAS/webuploader/server/upload/"

##--UI---------------------------------------------
ui <- bootstrapPage(
  useShinyjs(),
  htmlTemplate("index.html")
  
)

##--server---------------------------------------------

server <- function(input, output,session) {
  
  ## Global initialization--------------------------------------------------
  
  setClass("PISASDS", 
           slots = c(
             name = "character",
             species = "character",
             format = "character",
             protocols = "character",
             description = "character",
             create.time="character",
             count = "dgCMatrix",
             meta = "data.frame"
           )
  )
  
  setClass("PISAS", 
           slots = c(
             project = "character", 
             seed = "numeric",
             create.time="character",
             species = "character",
             step="character",
             type="character",
             dataset.list = "data.frame",
             array= "list",
             ori.result = "list",
             pro.result = "list",
             integrate = "list"
           )
  )
  
  pipelines <- reactiveValues()
  PISAS_syn <- reactiveValues()
  #ProjectList <- data.frame() #完整的dataset list
  selectList <- data.frame() #选择的dataset list
  daFormat <- "exp" 
  #最终状态
  #example = 0:new;1:example;
  pipelines <- reactiveValues(example=0,username=NULL) 
  
  PISAS_final <- NULL
  ispublic <- FALSE
  
  VaryPro <- reactiveValues(total = 0)
  VaryDS <- reactiveValues(total = 0)
  VaryDT <- reactiveValues(total = 0)
  Info_out <- reactiveValues(out_pro=c(),out_int=c(),out_dim=c(),out_file=c())
  
  ## Status --------------------------------------------------
  observe({
    
    query <- parseQueryString(session$clientData$url_search)
    #query <- parseQueryString("?o=6&x=0721e004f7aad93c29bbfa95c2f26602")
    
    #login
    if (!is.null(query[['x']])) {
      memid <- query[['x']]
      uuid <- query[['u']]
      style <- as.numeric(query[['s']])
      
      memid2 <- as.numeric(memid)
      uuid2 <- as.character(uuid)
      # nres <- checklogin(memid2,pass)
      
      if(!is.na(memid2)){
        shinyjs::show("exbut")
        pipelines$username <<- memid2
        # shinyalert("Thanks!", "You have successfully logged in!", type = "success")
      }else{
        # shinyalert("Note!", "User does not exist, please register first! You will log in with ip!", type = "success")
      }
      
      if(!is.na(uuid2) && !is.na(style)){
        if(style == 1){
          updateTextInput(
            session = session,
            inputId = "expFTPurl",
            value = paste0("UUID:",uuid2)
          )
        }else if(style==2){
          updateTextInput(
            session = session,
            inputId = "FTPurl",
            value = paste0("UUID:",uuid2)
          )
        }
      }
    }else{
      #shinyalert("Note!", "You are logging in as guest mode and will log in as ip!", type = "success")
      shinyjs::hide("exbut")
      
    }
    
    if(is.null(pipelines$username)){
      
      runjs('var user_ips= returnCitySN["cip"];Shiny.onInputChange("user_ips",user_ips);')
      
      #If you are only using it alone, you do not need to obtain the ip. 
      #pipelines$username <<- "myip"
      #If multiple people use and need to distinguish folders by IP, use the following code
      
      pipelines$username <<- ip2long(input$user_ips)
      
      if(is.null(input$user_ips)) pipelines$username <<- "local"
      
    }
  })
  
  output$users <- renderText({
    return(pipelines$username)
  })
  
  
  
  shinyInput = function(FUN, len, id, label = NULL,...) {
    inputs = character(len)
    for (i in seq_len(len)) {
      inputs[i] = as.character(FUN(paste0(id, i), label = label, ...))
    }
    inputs
  }
  
  
  shinyValue = function(id, len) {
    
    ss <- unlist(lapply(seq_len(len), function(i) {
      value = input[[paste0(id, i)]]
      if (is.null(value)) NA else value
    }))
    return(ss)
  }
  
  
  ##>Public Project------------------------------------------
  
  observeEvent(input$view_PBMC, {
    
    withProgress(message = 'Load the public data', value = 0.1, {
      for (i in 1:8) {
        # Each time through the loop, add another row of data. This a stand-in
        # for a long-running computation.
        
        # Increment the progress bar, and update the detail text.
        incProgress(0.1, detail = "This may take a while...")
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      #rm(list = ls())
      RData_path <- c("www/task/public/PBMC.RData")
      load(RData_path)
      PISAS_syn$PISAS_pro <<- PISAS_pro
      ispublic <<- TRUE
      
      if(PISAS_pro@type == "Single"){
        shinyjs::show("home_next")
        shinyjs::hide("btn_goMarge")
        shinyjs::show("drIs")
        shinyjs::hide("intIs")
      }else{
        shinyjs::show("btn_goMarge")
        shinyjs::hide("home_next")
        shinyjs::hide("drIs")
        shinyjs::show("intIs")
      }
      
    })
    shinyalert("OK!", "The project has been loaded successfully", type = "success")
  })
  
  observeEvent(input$view_HPIC, {
    
    withProgress(message = 'Load the public data', value = 0.1, {
      for (i in 1:8) {
        # Each time through the loop, add another row of data. This a stand-in
        # for a long-running computation.
        
        # Increment the progress bar, and update the detail text.
        incProgress(0.1, detail = "This may take a while...")
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      #rm(list = ls())
      RData_path <- c("www/task/public/HPIC.RData")
      load(RData_path)
      PISAS_syn$PISAS_pro <<- PISAS_pro
      
      ispublic <<- TRUE
      
      if(PISAS_pro@type == "Single"){
        shinyjs::show("home_next")
        shinyjs::hide("btn_goMarge")
        shinyjs::show("drIs")
        shinyjs::hide("intIs")
      }else{
        shinyjs::show("btn_goMarge")
        shinyjs::hide("home_next")
        shinyjs::hide("drIs")
        shinyjs::show("intIs")
      }
      
    })
    shinyalert("OK!", "The project has been loaded successfully", type = "success")
  })
  
  observeEvent(input$view_GSE84133, {
    
    withProgress(message = 'Load the public data', value = 0.1, {
      for (i in 1:8) {
        # Each time through the loop, add another row of data. This a stand-in
        # for a long-running computation.
        
        # Increment the progress bar, and update the detail text.
        incProgress(0.1, detail = "This may take a while...")
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      #rm(list = ls())
      RData_path <- c("www/task/public/GSE84133.RData")
      load(RData_path)
      PISAS_syn$PISAS_pro <<- PISAS_pro
      
      ispublic <<- TRUE
      
      if(PISAS_pro@type == "Single"){
        shinyjs::show("home_next")
        shinyjs::hide("btn_goMarge")
        shinyjs::show("drIs")
        shinyjs::hide("intIs")
      }else{
        shinyjs::show("btn_goMarge")
        shinyjs::hide("home_next")
        shinyjs::hide("drIs")
        shinyjs::show("intIs")
      }
      
    })
    shinyalert("OK!", "The project has been loaded successfully", type = "success")
  })
  
  ##>Personal Project List----
  Project_ID <- c()
  
  ProList <- reactive({
    if(VaryPro$total >0){
      message("Project list : Changed")
    }
    list <- c("www/task/Projectlist.Rds")
    if (!(file.exists(list))){
      return(NULL)
    }else{
      prolist <- readRDS(list)
      usernameTmp <- as.character(pipelines$username)
      #message("3",pipelines$username)
      df <- subset(prolist,username == usernameTmp)
      len <- nrow(df)
      
      if(len>0){
        
        ranids <- stri_rand_strings(1,3,'[a-z]')
        
        ProList <- data.frame(
          Project = df$Project,
          Type = df$Type,
          Number = df$Number,
          Seed = df$Seed,
          Size=df$Size,
          Create.time = df$create.time,
          Update.time = df$update.time,
          Delete = paste0('<a class="btn btn-md btn-danger" onclick="Shiny.onInputChange(\'delPro\',\'',df$ID,'\')"><span class="glyphicon glyphicon-trash"></span></a>'),
          Export = paste0('<a class="btn btn-md btn-info" onclick="Shiny.onInputChange(\'exPro\',\'',paste0("E",ranids,"_",df$ID),'\')"><span class="glyphicon glyphicon-download-alt"></span></a>'),
          Load = paste0('<a class="btn btn-md btn-success" onclick="Shiny.onInputChange(\'loadPro\',\'',paste0("L",ranids,"_",df$ID),'\')"><span class="glyphicon glyphicon-play"></span></a>'),
          stringsAsFactors = FALSE
        )
        return(ProList)
      }else{
        return(NULL)
      }
    }
  })
  
  output$proTab <- DT::renderDataTable(ProList(), server = FALSE, escape = FALSE, selection = 'none')
  
  ## >> Load project ----
  
  observeEvent(input$loadPro, {
    load_ids <- unlist(strsplit(input$loadPro, "_"))[2]
    withProgress(message = 'Load the public data', value = 0.1, {
      for (i in 1:8) {
        # Each time through the loop, add another row of data. This a stand-in
        # for a long-running computation.
        
        # Increment the progress bar, and update the detail text.
        incProgress(0.1, detail = "This may take a while...")
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      #rm(list = ls())
      prodir <- paste0("www/task/",pipelines$username,"/project/",load_ids,".RData")
      message(prodir)
      load(prodir)
      PISAS_syn$PISAS_pro <<- PISAS_pro
      if(PISAS_pro@type == "Single"){
        shinyjs::show("home_next")
        shinyjs::hide("btn_goMarge")
        shinyjs::show("drIs")
        shinyjs::hide("intIs")
      }else{
        shinyjs::show("btn_goMarge")
        shinyjs::hide("home_next")
        shinyjs::hide("drIs")
        shinyjs::show("intIs")
      }
    })
    VaryPro$total <- VaryPro$total+1
    Project_ID <<- load_ids
    shinyalert("OK!", "The project has been loaded successfully", type = "success")
  })
  
  ## >> Delete projest ----
  observeEvent(input$delPro, {
    load_ids <- input$delPro
    username <- pipelines$username
    withProgress(message = 'Delete the public data', value = 0.1, {
      for (i in 1:8) {
        # Each time through the loop, add another row of data. This a stand-in
        # for a long-running computation.
        
        # Increment the progress bar, and update the detail text.
        incProgress(0.1, detail = "This may take a while...")
        
        # Pause for 0.1 seconds to simulate a long computation.
        Sys.sleep(0.1)
      }
      
      prodir <- paste0("www/task/",username,"/project/",load_ids,".RData")
      message(prodir)
      #delete the file
      if(file.exists(prodir)){
        file.remove(prodir)
      }
      #change the prolist
      prolist <- readRDS("www/task/Projectlist.Rds")
      usernameTmp <- as.character(username)
      message(usernameTmp)
      message(load_ids)
      df <- prolist[!(prolist$ID == load_ids & prolist$username==usernameTmp),]
      rownames(df) <- NULL
      saveRDS(df,file = "www/task/Projectlist.Rds")
      VaryPro$total <- VaryPro$total+1
    })
    
    shinyalert("OK!", "The project has been deleted successfully", type = "success")
  })
  
  ## >> Export Project ----
  global <- reactiveValues(response = FALSE,load_ids=c(),ProRds=NULL)
  
  
  observeEvent(input$exPro, {
    username <- pipelines$username
    load_ids <- unlist(strsplit(input$exPro, "_"))[2]
    
    withProgress(message = 'Load the public data', value = 0.1, {
      
      prodir <- paste0("www/task/",username,"/project/",load_ids,".RData")
      
      incProgress(0.4, detail = "Load the project")
      load(prodir)
      
      incProgress(0.7, detail = "Export to a specific format")
      Result<-NULL
      if(input$exportFormat == "Seurat"){
        Result <- PISAS_pro@integrate$Combind
      }else if(input$exportFormat == "Monocle"){
        Result <- as.cell_data_set(PISAS_pro@integrate$Combind)
      }else if(input$exportFormat == "SCE"){
        Result <- as.SingleCellExperiment(PISAS_pro@integrate$Combind)
      }
      VaryPro$total <- VaryPro$total+1
      shinyalert("Confirmation",
                 "Do you want to download the data?",
                 type = "success",
                 callbackR = function(x) {
                   global$response <- x
                   global$load_ids <- load_ids
                   global$ProRds <- Result
                 },
                 showCancelButton = TRUE
      )
      
    })
  })
  
  observeEvent(global$response,{
    if(global$response){
      shinyjs::runjs("document.getElementById('download_export').click();")
      global$response <- FALSE
    }
  })
  
  output$download_export <- downloadHandler(
    
    filename = function() {
      paste0(global$load_ids,".Rds")
    },
    content = function(file) {
      withProgress(message = 'Saving the data...', value = 0.9, {
        saveRDS(global$ProRds,file)
      })
    }
  )
  
  
  
  #---------------> Upload New Files-------------------------------------
  
  observeEvent(input$data_format1, {
    daFormat <<- "exp"
  })
  
  observeEvent(input$data_format2, {
    daFormat <<- "10x"
  })
  
  observeEvent(input$way2, {
    
    if(input$way2 == "Local"){
      
      shinyjs::show("localWay")
      shinyjs::hide("ftpWay")
      shinyjs::hide("x10uuid")
      
    }else if(input$way2 == "Large"){
      shinyjs::hide("localWay")
      shinyjs::show("x10uuid")
      shinyjs::hide("ftpWay")
      
    }else{
      shinyjs::hide("localWay")
      shinyjs::show("ftpWay")
      shinyjs::hide("x10uuid")
    }
  })
  
  observeEvent(input$expFtp, {
    
    if(input$expFtp == "Local"){
      
      shinyjs::show("explocal")
      shinyjs::hide("expurl")
      shinyjs::hide("explarge")
    }else if(input$expFtp == "FTP"){
      
      shinyjs::hide("explocal")
      shinyjs::show("expurl")
      shinyjs::hide("explarge")
      
    }else{
      shinyjs::hide("explocal")
      shinyjs::show("explarge")
      shinyjs::hide("expurl")
    }
  })
  
  uploadFTP <- reactiveValues(status = 0,path="",exp="",meta="")
  
  output$fileMes <- renderText({
    Info_out$out_upload
  })
  
  #### >> Normal matrix ####
  #' button "Get file"
  observeEvent(input$btn_expftp, {
    tryCatch(
      {
        expFTPurl<- input$expFTPurl
        
        metaFTPurl<- input$metaFTPurl
        
        message(expFTPurl)
        
        Info_out$out_upload <- paste0("Getting file!<br>")
        
        withProgress(message = 'downloading the gene expression matrix...', value = 0.4, {
          
          if(grepl("https://|ftp://|http://",expFTPurl)){
            #FTP
            require("RCurl")
            mtx <- strsplit(expFTPurl,"/")[[1]]
            mtx <- mtx[length(mtx)]
            path1 <- paste0(taskdir,"public/",mtx)
            
            if(!file.exists(path1)){
              download.file(expFTPurl,path1)
              if(!file.exists(path1)){
                uploadFTP$status <- "2"
                Info_out$out_upload <- paste0(Info_out$out_upload ,"File download failed!<br>")
              }
            }else{
              uploadFTP$status <- "3"
              uploadFTP$exp <- path1
              setProgress(value = 1,message = paste0('Download ',mtx))
              Info_out$out_upload <- paste0(Info_out$out_upload ,paste0(mtx,' downloaded successfully!'))
            }
          }else if(grepl("UUID:",expFTPurl)){
            #UUID
            uuid <- gsub("UUID:","",expFTPurl)
            
            path1 <- paste0(uuiddir,uuid)
            message(path1)
            files <- list.files(path1)
            ms <- grepl("matrix",files)
            mtx <- c()
            if(sum(ms) == 1){
              mtx <- paste0(path1,"/",files[ms])
            }else{
              Info_out$out_upload <- paste0(Info_out$out_upload ,"Please make sure that the file group corresponding to the UUID contains a single matrix file!")
              uploadFTP$status <- 2
            }
            
            if(!file.exists(mtx)){
              Info_out$out_upload <- paste0(Info_out$out_upload ,"Please make sure the file has been uploaded successfully!")
              
              uploadFTP$status <- 2
            }else{
              uploadFTP$status <- "3"
              uploadFTP$exp <- mtx
              Info_out$out_upload <- paste0(Info_out$out_upload ,paste0('Meta downloaded successfully!'))
              
            }
          }else{
            #Local
            path1 <- expFTPurl
            if(!file.exists(path1)){
              Info_out$out_upload <- "Please ensure that the path contains single matrix file!"
              return(NULL)
            }else{
              uploadFTP$status <- "3"
              uploadFTP$exp <- path1
              Info_out$out_upload <- paste0(Info_out$out_upload ,paste0('Meta downloaded successfully!'))
              
            }
          }
        })
        
        if(metaFTPurl == ""){
          return(NULL)
        }
        withProgress(message = 'downloading the meta data...', value = 0.4, {
          if(grepl("https://|ftp://|http://",metaFTPurl)){
            # 1. FTP
            require("RCurl")
            mtx <- strsplit(metaFTPurl,"/")[[1]]
            mtx <- mtx[length(mtx)]
            path1 <- paste0(taskdir,"public/",mtx)
            
            if(!file.exists(path1)){
              download.file(metaFTPurl,path1)
              if(!file.exists(path1)){
                uploadFTP$status <- "2"
                Info_out$out_upload <- paste0(Info_out$out_upload ,"File download failed!<br>")
                
              }
            }else{
              uploadFTP$status <- "3"
              uploadFTP$meta <- path1
              setProgress(value =1,message = paste0('Download ',mtx))
              Info_out$out_upload <- paste0(Info_out$out_upload ,paste0(mtx,' downloaded successfully!'))
              
            }
          }else if(grepl("UUID:",metaFTPurl)){
            # 2. UUID
            uuid <- gsub("UUID:","",metaFTPurl)
            path1 <- paste0(uuiddir,uuid)
            files <- list.files(path1)
            ms <- grepl("meta",files)
            mtx <- c()
            if(sum(ms) == 1){
              mtx <- paste0(path1,"/",files[ms])
            }else{
              Info_out$out_upload <- paste0(Info_out$out_upload ,"Please ensure that the path contains single meta file!")
              return(NULL)
            }
            
            if(!file.exists(mtx)){
              return(NULL)
            }else{
              uploadFTP$status <- "3"
              uploadFTP$meta <- mtx
              Info_out$out_upload <- paste0(Info_out$out_upload ,paste0(mtx,' downloaded successfully!'))
              
            }
          }else{
            # 3. Local
            path1 <- metaFTPurl
            if(!file.exists(path1)){
              Info_out$out_upload <- paste0(Info_out$out_upload ,"Please ensure that the path contains single meta file!")
              
              return(NULL)
            }else{
              uploadFTP$status <- "3"
              uploadFTP$meta <- path1
              Info_out$out_upload <- paste0(Info_out$out_upload ,paste0(path1,' downloaded successfully!'))
              
            }
          }
        })
      },
      error = function(e) {
        Info_out$out_upload <- paste0(Info_out$out_upload,"<p style='color:red'>Stop<br>",safeError(e),"<br>","Failed!</p><br>")
        return(NULL)
        #stop(safeError(e))
      }
    )
  })
  
  ## >>> Preview ####
  
  exp.gal <- data.table()
  
  scRNA.table <- reactive({
    if(input$expFtp == "Local"){
      if(is.null(input$file1)) return(NULL)
      
      filepath <- input$file1$datapath
    }else{
      if(uploadFTP$exp!="" && uploadFTP$status == "3"){
        
        filepath <- uploadFTP$exp
      }else{
        
        return(NULL)
      }
    }
    scRNA.table <- NULL
    withProgress(message = 'Getting the data\n', value = 0.5, {
      
      tryCatch(
        {
          #scRNA.table<- read.csv(input$file1$datapath, header=input$header1, sep=input$sep1, quote=input$quote1,row.names = 1,check.names=FALSE)
          
          if(input$header1=="auto"){
            header1 <- "auto"
          }else if(input$header1=="TRUE"){
            header1 <- TRUE
          }else if(input$header1=="FALSE"){
            header1 <- FALSE
          }
          
          if(grepl("zip",filepath)){
            scRNA.table <- fread.zip(filepath,sep = input$sep1,header = header1, quote=input$quote1) 
            
          }else{
            scan1 <- function(what, ...) scan(filepath, nmax = 1, what = what, 
                                              quiet = TRUE, ...)
            if (scan1(character()) == "%%MatrixMarket"){
              Info_out$out_upload <- c("file is a Sparse Matrix file! Please go to 'Data Format 2: Sparse matrix' tab to load the data.")
              return(NULL)
            }else{
              scRNA.table<- fread(filepath, sep = input$sep1,header = header1, quote=input$quote1,check.names=FALSE)
              
            }
          }
        },
        error = function(e){
          # return a safeError if a parsing error occurs
          return(safeError('Reading data error! Make sure that the file you upload is in the correct format. After that, modify the values of the parameters on the right and wait to read the data again.'))
        })
    })
    
    exp.gal <<- scRNA.table
    
    return(scRNA.table)
  })
  
  meta.gal <- data.frame()
  
  primer.meta <- reactive({
    if(input$expFtp == "Local"){
      if(is.null(input$file2)) return(NULL)
      
      filepath <- input$file2$datapath
    }else{
      if(uploadFTP$status == "3" && uploadFTP$meta != ""){
        filepath <- uploadFTP$meta
      }else{
        return(NULL)
      }
    }
    primer.meta <- NULL
    tryCatch(
      {
        primer.meta<- read.csv(filepath, header=input$header2, sep=input$sep2, quote=input$quote2,row.names = 1,check.names=FALSE)
        #primer.meta <- fread(input$file2$datapath, sep = input$sep2,header = input$header2, quote=input$quote2,check.names=FALSE,na.strings = "NA")
      },
      error = function(e) {
        return(safeError(e))
      }
    )
    meta.gal <<- primer.meta
    return(primer.meta)
  })
  
  # matrix
  output$table1 <- DT::renderDataTable({
    req(scRNA.table())
    withProgress(message = 'Preview the data\n', value = 0.5, {
      validate(
        need(dim(scRNA.table())[2] > 1, "**Please reselect options(header/separator/quote) to load data correctly")
      )
      validate(
        need(sum(colnames(scRNA.table())[1:3]==c("V1","V2","V3"))!=3,"**Please reselect Header = TRUE")
      )
      
      scRNA <- head(scRNA.table())[,1:6]
    })
    return(DT::datatable(scRNA, options = list(orderClasses = TRUE,scrollX = TRUE)))
  })
  # barcode
  output$table2 <- DT::renderDataTable({
    req(primer.meta())
    validate(
      need(dim(primer.meta())[2] > 1, "**Please reselect options(header/separator/quote) to load data correctly")
    )
    
    primer <- head(primer.meta())
    return(DT::datatable(primer, options = list(orderClasses = TRUE,scrollX = TRUE)))
  })
  
  #### >> Sparse matrix ####
  # button "Get file"
  observeEvent(input$btn_ftp, {
    tryCatch(
      {
        url<- input$FTPurl
        
        Info_out$out_upload <- paste0("Connecting URL: ",url,"<br>")
        withProgress(message = 'downloading the files...', value = 0.4, {
          if(grepl("https://|ftp://|http://",url)){
            #FTP
            #res<- download_FTP(url,taskdir)
            
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
              
              if(file.exists(paste0(filepath,"/matrix.mtx.gz")) && 
                 file.exists(paste0(filepath,"/barcodes.tsv.gz")) &&
                 file.exists(paste0(filepath,"/genes.tsv.gz"))){
                
                uploadFTP$status <- 3
                uploadFTP$path <- filepath
              }else{
                #not exist-download
                if(!file.exists(filepath)){
                  dir.create(file.path(filepath),recursive = TRUE)
                }
                withProgress(message = 'downloading the matrix...', value = 0.2, {
                  
                  Info_out$out_upload <- paste0(Info_out$out_upload,"Downloading file: ",paste0(url,mtx),"<br>")
                  
                  download.file(paste0(url,mtx),paste0(filepath,"/matrix.mtx.gz"))
                })
                withProgress(message = 'downloading the barcodes...', value = 0.4, {
                  Info_out$out_upload <- paste0(Info_out$out_upload,"Downloading file: ",paste0(url,barcodes),"<br>")
                  
                  download.file(paste0(url,barcodes),paste0(filepath,"/barcodes.tsv.gz"))
                })
                withProgress(message = 'downloading the genes...', value = 0.8, {
                  Info_out$out_upload <- paste0(Info_out$out_upload,"Downloading file: ",paste0(url,genes),"<br>")
                  
                  download.file(paste0(url,genes),paste0(filepath,"/genes.tsv.gz"))
                })
                
                if(file.exists(paste0(filepath,"/matrix.mtx.gz")) && 
                   file.exists(paste0(filepath,"/barcodes.tsv.gz")) &&
                   file.exists(paste0(filepath,"/genes.tsv.gz"))){
                  
                  uploadFTP$status <- 3
                  uploadFTP$path <- filepath
                }else{
                  Info_out$out_upload <- paste0(Info_out$out_upload,"File download failed!<br>")
                  
                  uploadFTP$status <- 2
                }
              }
            }else{
              Info_out$out_upload <- "Please ensure that the ftp link contains single matrix.mtx(.gz),barcodes.tsv(.gz) and features.tsv(.gz) / genes.tsv(.gz)!"
              uploadFTP$status <- 1
            }
          }else if(grepl("UUID:",url)){
            #UUID
            uuid <- gsub("UUID:","",url)
            path1 <- paste0(uuiddir,uuid)
            files <- list.files(path1)
            uploadFTP$status <- 3
            uploadFTP$path <- path1
          }else{
            #Local
            path1 <- url
            if(!file.exists(path1)){
              Info_out$out_upload <- "Please ensure that the path contains single matrix.mtx(.gz),barcodes.tsv(.gz) and features.tsv(.gz) / genes.tsv(.gz)!"
              return(NULL)
            }else{
              uploadFTP$status <- 3
              uploadFTP$path <- path1
            }
          }
          
          setProgress(value =1,message = 'Finished! Wait to preview the data')
        })
      },
      error = function(e) {
        Info_out$out_upload <- paste0(Info_out$out_upload,"<p style='color:red'>Stop<br>",safeError(e),"<br>","Failed!</p><br>")
        return(NULL)
        #stop(safeError(e))
      }
    )
  })
  
  
  ## >>>Preview ####
  # mtx
  output$tab10x_matrix <- DT::renderDataTable({
    withProgress(message = 'Preview the data\n', value = 0.5, {
      
      if(input$way2 == "Local"){
        
        if(is.null(input$file10xm)){return(NULL)}
        
        matrix.mtx<- readMM(input$file10xm$datapath)
        
      }else{
        
        if(uploadFTP$path!="" && uploadFTP$status == "3"){
          files <- list.files(uploadFTP$path)
          ms <- grepl("matrix",files)
          matrix.mtx<- readMM(paste0(uploadFTP$path,"/",files[ms]))
          
        }else{
          return(NULL)
        }
      }
      return(DT::datatable(head(as.data.table(matrix.mtx[,1:15])), options = list(orderClasses = TRUE,scrollX = TRUE)))
      
      
    })
  })
  # genes
  output$tab10x_gene <- DT::renderDataTable({
    if(input$way2 == "Local"){
      if(is.null(input$file10xg)){return(NULL)}
      gene.csv<- read.csv(input$file10xg$datapath,header = FALSE,sep="\t")
      return(DT::datatable(gene.csv, options = list(orderClasses = TRUE,scrollX = TRUE)))
      
    }else{
      if(uploadFTP$path!="" && uploadFTP$status == "3"){
        files <- list.files(uploadFTP$path)
        ms <- grepl("gene|feature",files)
        gene.csv<- read.csv(paste0(uploadFTP$path,"/",files[ms]),header = FALSE,sep="\t")
        return(DT::datatable(gene.csv, options = list(orderClasses = TRUE,scrollX = TRUE)))
        
      }else{
        return(NULL)
      }
    }
  })
  # barcodes
  
  output$tab10x_barcode <- DT::renderDataTable({
    if(input$way2 == "Local"){
      if(is.null(input$file10xb)){return(NULL)}
      barcode.csv<- read.csv(input$file10xb$datapath,header = FALSE,sep="\t")
    }else{
      
      if(uploadFTP$path!="" && uploadFTP$status == "3"){
        files <- list.files(uploadFTP$path)
        ms <- grepl("barcode",files)
        barcode.csv<- read.csv(paste0(uploadFTP$path,"/",files[ms]),header = FALSE,sep="\t")
      }else{
        return(NULL)
      }
    }
    meta.gal <<- barcode.csv
    return(DT::datatable(as.data.frame(head(barcode.csv)), options = list(orderClasses = TRUE,scrollX = TRUE)))
    
  })
  
  #### >>Create DataSet ####
  #
  output$speciesSelect <- renderUI({
    
    db <- homologene::taxData
    group_lists <- as.list(c(db$name_txt,"Others"))
    selectInput("species", "Species:", group_lists,selected = "Homo sapiens")
  })
  
  # uplaod
  observeEvent(input$btn_upload, {
    if(input$DSName=="" | input$protocols == ""){
      shinyalert("Warning!", "Please enter required fields!", type = "warning")
      return(NULL)
    }
    username <- pipelines$username
    withProgress(message = 'Upload the dataset\n', value = 0.1, {
      
      list <- paste0(taskdir,"DSlist.Rds")
      create_time <- Sys.time()
      create.time <- format(create_time,format='%Y/%m/%d %H:%M:%S')
      if (!(file.exists(list))){
        dslist <- data.frame(name=input$DSName,username=username,species=input$species,
                             protocols = input$protocols,create.time= create.time)
        saveRDS(dslist,file=list)
      }else{
        dslist <- readRDS(list)
        if(sum(dslist$name == input$DSName & dslist$username == username) >0 ){
          
          shinyalert("Warning!", "Do not set the same dataset name!", type = "warning")
          return()
        }else{
          newlist <- data.frame(name=input$DSName,username=username,species=input$species,
                                protocols = input$protocols,create.time= create.time)
          dslist <- rbind(dslist,newlist)
          
        }
      }
      filepath <- paste0("www/task/",username)
      filepath10x <- paste0("www/task/",username,"/x10tmp/")
      
      if (!(file.exists(filepath))){
        dir.create(file.path(filepath),recursive = TRUE)
      }
      if (!(file.exists(filepath10x))){
        dir.create(file.path(filepath10x),recursive = TRUE)
      }
      
      setProgress(value = 0.4,detail="reading file ...")
      
      if(daFormat == "exp"){
        
        if(dim(exp.gal)[1]==0){
          shinyalert("Warning!", "Ensure upload the cell-gene expression matrix!", type = "warning")
          return(NULL)
        }
        scRNA.table <- As.DF(exp.gal)
      }else if(daFormat == "10x"){
        
        if (input$way2 == "Local"){
          if(is.null(input$file10xm) | is.null(input$file10xg) | is.null(input$file10xb)) return(NULL)
          f1 <- input$file10xm$datapath
          f2 <- input$file10xg$datapath
          f3 <- input$file10xb$datapath
          file.copy(f1, paste0(filepath10x,"/matrix.mtx"),overwrite = TRUE)
          file.copy(f2, paste0(filepath10x,"/genes.tsv"),overwrite = TRUE)
          file.copy(f3, paste0(filepath10x,"/barcodes.tsv"),overwrite = TRUE)
        }else{
          if(uploadFTP$path!="" && uploadFTP$status == "3"){
            filepath10x <- uploadFTP$path
            
          }else{
            return(NULL)
          }
        }
        message(filepath10x)
        
        #data10x <- Read10X(filepath10x,gene.column=input$geneColumn,unique.features=input$uniqueFeatures,strip.suffix = input$stripSuffix)
        
        data10x <- Read10X(filepath10x,gene.column=input$geneColumn)
        
        scRNA.table <- data10x
        #file.remove(filepath10x)
      }
      primer.meta <- meta.gal
      
      setProgress(value = 0.7,detail="Storing the dataset ...")
      
      saveRDS(dslist,file=list)
      VaryDT$total <- VaryDT$total +1
      
      # create the dataset
      pisasObj <- new("PISASDS", name = input$DSName, format=daFormat ,species = input$species,protocols = input$protocols,
                      description=input$description,create.time = create.time,count = scRNA.table,meta = primer.meta)
      saveRDS(pisasObj,file=paste0(filepath,"/",input$DSName,".Rds"))
      
      # add the dataset list
      setProgress(value = 0.9,detail="refreshing the database...")
      
      setProgress(value =1,message = 'Finished!')
      Sys.sleep(1)
    })
    
  })
  # reset
  observeEvent(input$btn_reset, {
    
    reset('file1')
    reset('file2')
    
    reset('file10xm')
    reset('file10xg')
    reset('file10xb')
    
    reset("DSName")
    reset("species")
    reset("protocols")
    reset("description")
  })
  
  ## >New Project-----------------------
  
  checkID <- c()
  #>>Tab list ----
  #refresh dataset lists
  TabList <- reactive({
    if(VaryDT$total > 0){
      message("Dataset list chenged")
    }
    list <- c("www/task/DSlist.Rds")
    username <- pipelines$username
    if (!(file.exists(list))){
      return(NULL)
    }else{
      withProgress(message = 'Refresh the datasets...', value = 0.5, {
        dslist <- readRDS(list)
        usernameTmp <- as.character(username)
        df <- subset(dslist,username == usernameTmp)
        len <- nrow(df)
        checkID <<- stri_rand_strings(1,3,'[a-z]')
        TabList <- data.frame(
          Checkbox = shinyInput(checkboxInput, len, checkID, value = FALSE),
          Dataset = df$name,
          Species = df$species,
          Protocol = df$protocols,
          Time = df$create.time,
          stringsAsFactors = FALSE
        )
        return(TabList)
      })
    }
  })
  
  #Tab List
  output$tabdslist <- DT::renderDataTable({
    TabList()
  }, server = FALSE, escape = FALSE, selection = 'none', options = list(
    preDrawCallback = JS('function() { Shiny.unbindAll(this.api().table().node()); }'),
    drawCallback = JS('function() { Shiny.bindAll(this.api().table().node()); } ')
  ))
  
  #print filetime in refreshdate1
  
  output$selectlist <- renderPrint({
    TabList <- TabList()
    len <-nrow(TabList)
    validate(
      need(len > 0, "Your dataset list is empty, please upload the dataset first!")
    )
    
    checkbox <- shinyValue(checkID, len)
    
    Tab_TRUE <- data.frame(checkbox = checkbox,Dataset=TabList$Dataset,
                           Species=TabList$Species,Protocol=TabList$Protocol,
                           Time=TabList$Time)
    selectList <<- subset(Tab_TRUE, checkbox == 'TRUE')
    s <- as.character.Array(selectList$Dataset)
    
    if (length(s)) {
      cat('These datasets were selected:\t')
      cat(s,sep = ", ")
    }
  })
  
  ## delete the selected datasets
  observeEvent(input$btn_delDT, {
    
    VaryDT$total <- VaryDT$total +1
    
    withProgress(message = 'Delete the selected datasets\n', value = 0.1, {
      
      if(nrow(selectList)==0){
        shinyalert("Warning!", "Please select the datasets to be deleted first!", type = "warning")
        return(NULL)
      }
      
      list <- c("www/task/DSlist.Rds")
      username <- pipelines$username
      if (!(file.exists(list))){
        return(NULL)
      }else{
        dslist <- readRDS(list)
        select <- FALSE
        usernameTmp <- as.character(username)
        for (i in 1:nrow(selectList)) {
          seltmp <- dslist$name == as.character(selectList$Dataset[i]) & dslist$username == usernameTmp
          select <- select | seltmp
          dtpath <- paste0(taskdir,username,"/",as.character(selectList$Dataset[i]),".Rds")
          message(dtpath)
          file.remove(dtpath)
        }
        
        dslist <- dslist[!select,]
        rownames(dslist) <- NULL
        saveRDS(dslist,file=list)
        
        selectList <<-  data.frame()
      }
      
      setProgress(value =1,message = 'Success !')
      Sys.sleep(2)
    })
    
  })
  
  ##>>Create Project ----
  observeEvent(input$btn_create, {
    
    len <- nrow(selectList)
    message(paste0("Select ",len," datasets"))
    username <- pipelines$username
    if( len == 0){
      shinyalert("Warning!", "Please select the datasets first!", type = "warning")
      
      return(NULL)
    }
    if((input$Type == "Single" && len == 1) || (input$Type == "Multiple" && len > 1)){
      withProgress(message = 'loading the datasets...', value = 0.5, {
        
        dslist <- list()
        dataset_info <- data.frame()
        species <- ""
        for (i in 1:nrow(selectList)) {
          name=selectList[i,]$Dataset
          list <- paste0(taskdir,username,"/",name,".Rds")
          message(list)
          ds <- readRDS(list)
          species <- ds@species[1]
          expdata <- ds@count
          metadata <- ds@meta
          name <- ds@name
          
          dsinfo.1 <- data.frame(ID = i,Name=name,Groups=ncol(metadata),Step="Raw",
                                 Cells=ncol(expdata),Genes=nrow(expdata),
                                 Status = "",Preprocessing ="",
                                 stringsAsFactors = FALSE
          )
          dsinfo.2 <- data.frame(ID = i,Name=name,Groups=ncol(metadata),Step="Preprocessed",
                                 Cells="-",Genes="-",
                                 Status = paste0('<div class="tag status-no">Non-preprocessed</div>'),
                                 Preprocessing = paste0('<a href="#preprocessing" data-toggle="tab" class="btn btn-md btn-info" onclick="Shiny.onInputChange(\'doPro\',',i,')">GO</a>'),
                                 stringsAsFactors = FALSE
          )
          dsinfo.3 <- data.frame(ID = i,Name=name,Groups=ncol(metadata),Step="Integrated",
                                 Cells="-",Genes="-",
                                 Status = "",Preprocessing ="",
                                 stringsAsFactors = FALSE
          )
          dataset_info<- rbind(dataset_info,dsinfo.1,dsinfo.2,dsinfo.3)
          
          ds@count <- expdata
          
          ds@meta <- metadata
          
          dslist <- c(dslist,ds)
        }
        
        create_time <- Sys.time()
        create.time <- format(create_time,format='%Y/%m/%d %H:%M:%S')
        PISAS_syn$PISAS_pro <<- new("PISAS", project = input$projectName, seed = input$randomSeed,species = species,
                                    create.time=create.time,dataset.list = dataset_info,
                                    array=dslist,type = input$Type,step =character(length = 0))
        
      })
      
    }else{
      shinyalert("Warning!", "The number of datasets you choose needs to be consistent with the \"type\" option! Please switch \"Single\" to \"Multiple\" or switch \"Multiple\" to \"Single\"", 
                 type = "warning",closeOnClickOutside=TRUE)
      return()
    }
    if(input$Type == "Single"){
      shinyjs::show("home_next")
      shinyjs::hide("btn_goMarge")
      
    }else{
      shinyjs::show("btn_goMarge")
      shinyjs::hide("home_next")
    }
    
    ispublic <<- FALSE
    
    setProgress(value=1,message = 'Success!')
    Sys.sleep(0.5)
    
  })
  
  ## Result: DataSet ----
  output$DataSetLists <- renderDataTable({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    if(VaryDS$total >0){
      message("processed success")
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    sumRes <- PISAS_pro@dataset.list
    tmp <- subset(sumRes,ID != "Integrate")
    
    dtable <- datatable(tmp, rownames = FALSE, escape = FALSE,
                        options = list(
                          pageLength=25,
                          rowsGroup = list(0,1,2) # merge cells of column 1
                          
                        ))
    
    path <- paste0(getwd(),"/www/js/") # folder containing dataTables.rowsGroup.js
    dep <- htmltools::htmlDependency(
      "RowsGroup", "2.0.0",
      path, script = "dataTables.rowsGroup.js")
    dtable$dependencies <- c(dtable$dependencies, list(dep))
    
    return(dtable)
  })
  
  ### Pipeline---------------------------------------
  
  ## 1 Preprocessing--------------------------------
  
  cur_proids <- reactive({
    input$doPro
  })
  
  ## preview load dataset
  observeEvent(input$doPro, {
    cur_proid <- cur_proids()
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    isy <- PISAS_pro@dataset.list$ID==cur_proid & PISAS_pro@dataset.list$Step =="Preprocessed"
    nameID <- PISAS_pro@dataset.list$Name[isy]
    
    cur_data <- PISAS_pro@array[[cur_proid]]
    gnames <- rownames(cur_data@count)[2:4]
    
    output$preview_gnames <- renderPrint({
      
      s <- paste(gnames, collapse = ",")
      cat(s)
      cat("...")
    })
  })
  
  output$info_pro<-renderText(Info_out$out_pro)
  
  output$campSelector_keytype <- renderUI({
    if(is.null(cur_proids()))
      return()
    cur_proid <- cur_proids()
    PISAS_pro <- PISAS_syn$PISAS_pro
    cur_data <- PISAS_pro@array[[cur_proid]]
    species <- cur_data@species
    if(species == "Homo sapiens")
    {
      library(org.Hs.eg.db)
      key_list <- keytypes(org.Hs.eg.db)
    }else if(species == "Mus musculus"){
      library(org.Mm.eg.db)
      key_list <-keytypes(org.Mm.eg.db)
    }else if(species == "Arabidopsis thaliana"){
      library(org.At.tair.db)
      key_list <-keytypes(org.At.tair.db)
    }else if(species == "Xenopus (Silurana) tropicalis"){
      library(org.Xl.eg.db)
      key_list <-keytypes(org.Xl.eg.db)
    }else if(species == "Rattus norvegicus"){
      library(org.Rn.eg.db)
      key_list <-keytypes(org.Rn.eg.db)
    }else if(species == "Saccharomyces cerevisiae"){
      library(org.Sc.sgd.db)
      key_list <-keytypes(org.Sc.sgd.db)
    }else if(species == "Anopheles gambiae"){
      library(org.Ag.eg.db)
      key_list <-keytypes(org.Ag.eg.db)
    }else if(species == "Gallus gallus"){
      library(org.Gg.eg.db)
      key_list <-keytypes(org.Gg.eg.db)
    }else if(species == "Pan troglodytes"){
      library(org.Pt.eg.db)
      key_list <-keytypes(org.Pt.eg.db)
    }else if(species == "Bos taurus"){
      library(org.Bt.eg.db)
      key_list <-keytypes(org.Bt.eg.db)
    }else if(species == "Danio rerio"){
      library(org.Dr.eg.db)
      key_list <-keytypes(org.Dr.eg.db)
    }
    
    group_lists <- as.list(key_list)
    selectInput("keytype", "Keytype From:", group_lists,selected = "ENSEMBL")
  })
  
  output$nfes_sel <- renderUI({
    if(is.null(cur_proids()))
      return()
    cur_proid <- cur_proids()
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    isy <- PISAS_pro@dataset.list$ID==cur_proid & PISAS_pro@dataset.list$Step =="Raw"
    nmax <- PISAS_pro@dataset.list$Genes[isy]
    n <- 2000
    if(nmax<2000) n <- nmax
    numericInput("nfeatures", "Top features:",n,min = 0, max = nmax)
  })
  
  observeEvent(input$btn_prepro, {
    
    tryCatch(
      {
        withProgress(message = 'Preprocessing ...\n', value = 0, {
          if(is.null(PISAS_syn$PISAS_pro)){
            return(NULL)
          }
          Info_out$out_pro <- c("<h3>Preprocessing</h3>Start<br>")
          
          if(is.null(cur_proids()))
            return()
          
          cur_proid <- cur_proids()
          
          PISAS_pro <- PISAS_syn$PISAS_pro
          
          isy <- PISAS_pro@dataset.list$ID==cur_proid & PISAS_pro@dataset.list$Step =="Preprocessed"
          nameID <- PISAS_pro@dataset.list$Name[isy]
          
          setProgress(value = 0.1, detail = "Reading Dataset")
          Info_out$out_pro <- paste0(Info_out$out_pro,"Reading Dataset...<br>")
          
          if(ispublic){
            cur_data <- readRDS(paste0("www/task/public/",nameID,".Rds"))
          }else{
            cur_data <- PISAS_pro@array[[cur_proid]]
          }
          scRNA.data.ds <- cur_data@count
          
          # if(sum(dim(scRNA.data.ds)==c(1,1))==2){
          #   #只可能出现在个人项目中,为后面改方案留下接口
          #   #读取与处理后的count
          #    
          #   if(cur_proid > length(PISAS_pro@pro.result)){
          #     cur_data <- readRDS(paste0("www/task/",username,"/",cur_data@name,".Rds"))
          #     scRNA.data.ds <- cur_data@count
          #   }else{
          #     scRNA.data.ds <- GetAssayData(PISAS_pro@pro.result[[1]],slot = "counts")
          #   }
          # }
          
          primer.meta <- cur_data@meta
          project <- cur_data@name
          species <- cur_data@species
          tech <- cur_data@protocols
          step <- PISAS_pro@step
          type <- PISAS_pro@type
          format <- cur_data@format
          
          set.seed(PISAS_pro@seed)
          epsilon <- 1e-05
          do.downsampling <- FALSE  # for downsampling
          nGene.low <- as.numeric(input$mingenes)  # 200
          nUMI.low <- 500
          nCell.low <- as.numeric(input$mincells) #3
          filter.low.exp.cells <- TRUE
          is.transfer <- input$transfer
          
          do.cpp <- TRUE  # for data scaling
          regress.model <- "linear"  # for suppressing batch effect, linear, poisson, negbinom
          pattern.mito <- "^MT-"
          
          setProgress(value = 0.2, detail = "Conversion gene annotation ID")
          
          Info_out$out_pro <- paste0(Info_out$out_pro,"Conversion gene annotation ID...<br>")
          
          gene.names <-  rownames(scRNA.data.ds)
          gene.symbol <- c()
          disable("btn_prepro")
          
          nk <- 0.3 *length(gene.names)
          #' 判断是否为symbol命名，如果是则不转换，否则转换；如果是单个数据集分析则也不转换
          #' 转换的目的是为了之后的整合。
          #' 此处应该判断基因名称的形式，如果为symbol则判断，如果为ensemble 则转换为symbol再判断
          #' 
          if(is.transfer){
            
            if(species == "Homo sapiens")
            {
              is.spike <- grepl("^ERCC", gene.names)
              require(org.Hs.eg.db)
              gene.ids <- mapIds(org.Hs.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
            }else if(species == "Mus musculus"){
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Mm.eg.db)
              gene.ids <- mapIds(org.Mm.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
            }else if(species == "Arabidopsis thaliana"){
              
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.At.tair.db)
              gene.ids <- mapIds(org.At.tair.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }else if(species == "Xenopus (Silurana) tropicalis"){
              
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Xl.eg.db)
              gene.ids <- mapIds(org.Xl.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }else if(species == "Rattus norvegicus"){
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Rn.eg.db)
              gene.ids <- mapIds(org.Rn.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }else if(species == "Saccharomyces cerevisiae"){
              
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Sc.sgd.db)
              gene.ids <- mapIds(org.Sc.sgd.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }else if(species == "Anopheles gambiae"){
              
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Ag.eg.db)
              gene.ids <- mapIds(org.Ag.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }else if(species == "Gallus gallus"){
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Gg.eg.db)
              gene.ids <- mapIds(org.Gg.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }else if(species == "Pan troglodytes"){
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Pt.eg.db)
              gene.ids <- mapIds(org.Pt.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }else if(species == "Bos taurus"){
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Bt.eg.db)
              gene.ids <- mapIds(org.Bt.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }else if(species == "Danio rerio"){
              is.spike <- grepl("^ERCC-", gene.names)
              require(org.Dr.eg.db)
              gene.ids <- mapIds(org.Dr.eg.db, keys=gene.names, keytype= input$keytype, column="SYMBOL")
              
            }
            
            is.MT <- grepl("^MT-", gene.names)
            
            
            #将两种特殊情况重新命名
            gene.ids[is.spike] <- gene.names[is.spike]
            gene.ids[is.MT] <- gene.names[is.MT]
            
            # Deduplication
            new.ids <- uniquifyFeatureNames(gene.names,gene.ids)
            rownames(scRNA.data.ds) <- as.character(new.ids)
            
          }else{
            is.spike <- grepl("^ERCC", gene.names) #结果返回大量的TRUE和FALSE
            is.MT <- grepl("^MT-", gene.names)
          }
          
          # withProgress calls can be nested, in which case the nested text appears
          # below, and a second bar is shown.
          setProgress(value = 0.35, detail = "Generating the object...")
          Info_out$out_pro <- paste0(Info_out$out_pro,"Generating the object...<br>")
          
          scRNA <- CreateSeuratObject(counts = scRNA.data.ds, project = project, min.cells = nCell.low, min.features = nGene.low)
          
          scRNA_old <- CreateSeuratObject(counts = scRNA.data.ds)
          
          if(length(levels(scRNA_old@meta.data$orig.ident)) == 1){
            scRNA_old@meta.data$orig.ident <- as.factor(project)
          }
          
          scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
          scRNA_old[["percent.mt"]] <- PercentageFeatureSet(scRNA_old, pattern = "^MT-")
          
          scRNA[["percent.ERCC"]] <- PercentageFeatureSet(scRNA, pattern = "^ERCC")
          scRNA_old[["percent.ERCC"]] <- PercentageFeatureSet(scRNA_old, pattern = "^ERCC")
          
          #update metadata
          if(format != "10x"){
            if(dim(primer.meta)[1] > 0){
              tmpdt <- primer.meta[match(rownames(scRNA@meta.data),rownames(primer.meta)),]
              scRNA@meta.data <- cbind(scRNA@meta.data,tmpdt)
            }
          }
          
          projects <- rep(project,nrow(scRNA@meta.data))
          names(projects) <- colnames(x = scRNA)
          scRNA <- AddMetaData(
            object = scRNA,
            metadata = projects,
            col.name = 'project.PIAS'
          )
          
          techs <- rep(tech,nrow(scRNA@meta.data))
          names(techs) <- colnames(x = scRNA)
          scRNA <- AddMetaData(
            object = scRNA,
            metadata = techs,
            col.name = 'technology.PIAS'
          )
          scRNA@meta.data$species.PIAS <- species
          
          # isOutlier
          libsize.drop <- FALSE
          feature.drop <- FALSE
          spike.drop <- FALSE
          mito.drop <- FALSE
          
          if(input$libsize.drop){
            libsize.drop <- isOutlier(scRNA@meta.data$nCount_RNA, nmads=3, type="lower", log=TRUE,batch=scRNA@meta.data$orig.ident)
          }
          if(input$feature.drop){
            feature.drop <- isOutlier(scRNA@meta.data$nFeature_RNA, nmads=3, type="lower",log=TRUE,batch=scRNA@meta.data$orig.ident)
          }
          if(input$spike.drop){
            spike.drop <- isOutlier(scRNA@meta.data$percent.ERCC, nmads=3, type="higher",batch=scRNA@meta.data$orig.ident)
          }
          if(input$mito.drop){
            mito.drop <- isOutlier(scRNA@meta.data$percent.mt,nmads=3,type="higher",batch=scRNA@meta.data$orig.ident)
          }
          
          is.keep <- !(libsize.drop | feature.drop | spike.drop | mito.drop)
          if(length(is.keep) >1){
            scRNA <- scRNA[,is.keep]
          }
          
          filters <- match(rownames(scRNA_old@meta.data),rownames(scRNA@meta.data))
          scRNA_old@meta.data$filter[is.na(filters)] <- "low quality"
          scRNA_old@meta.data$filter[!is.na(filters)] <- "normal"
          
          Info_out$out_pro <- paste0(Info_out$out_pro,"Adding meta data...<br>")
          
          #' normalization
          #' 
          if (input$isnormalize != "none") {
            setProgress(value = 0.5, detail = "Normalizing")
            scRNA <- NormalizeData(scRNA,normalization.method = input$isnormalize)
            Info_out$out_pro <- paste0(Info_out$out_pro,"Normalizing...<br>")
          }
          # re-calculate the library size and number of expressed genes
          scRNA[["nCount_RNA"]] <- apply(scRNA@assays$RNA@data,2,sum)
          scRNA[["nFeature_RNA"]] <-apply(scRNA@assays$RNA@data,2,function(x){sum(x>0)})
          
          setProgress(value = 0.6, detail = "Finding variable features")
          Info_out$out_pro <- paste0(Info_out$out_pro,"Finding variable features...<br>")
          
          scRNA <- FindVariableFeatures(scRNA, selection.method = input$hvgfunc,nfeatures=input$nfeatures)
          
          setProgress(value = 0.7, detail = "Removing batch effect")
          Info_out$out_pro <- paste0(Info_out$out_pro,"Removing batch effect...<br>")
          
          if(length(levels(scRNA$orig.ident)) > 1){
            
            scRNA_hvg <- subset(x = scRNA, features = VariableFeatures(object = scRNA))
            
            SCE_seurat <- as.SingleCellExperiment(scRNA_hvg)
            
            if (input$batchEffect == "Combat") {
              scRNA@assays$RNA@scale.data <- as.matrix(batch.normalise.comBat(data = SCE_seurat))
            }else if(input$batchEffect == "RUVg"){
              scRNA@assays$RNA@scale.data <- as.matrix(batch.normalise.ruvg(data = SCE_seurat))
            }else if(input$batchEffect == "GLM"){
              scRNA@assays$RNA@scale.data <- as.matrix(batch.normalise.glm(data=SCE_seurat,batch=SCE_seurat$orig.ident,indi=SCE_seurat$Donor))
            }else if(input$batchEffect == "ScaleData"){
              scRNA <- ScaleData(scRNA)
            }
          }else{
            scRNA <- ScaleData(scRNA)
          }
          
          setProgress(value = 0.9, detail = "Print result table")
          Info_out$out_pro <- paste0(Info_out$out_pro,"Print result table...<br>")
          
          #过滤信息
          pro_result <- filter.info(id=cur_proid,dataset.name = project,
                                    species=species,tech=tech,scRNA.data=scRNA.data.ds,
                                    scRNA = scRNA,is.spike = is.spike)
          
          PISAS_pro@dataset.list$Cells[isy] <- pro_result$Final.cells
          PISAS_pro@dataset.list$Genes[isy] <- pro_result$Final.genes
          PISAS_pro@dataset.list$Status[isy] <-  paste0('<div class="tag status-preprocessed">Preprocessed</div>')
          
          VaryDS$total <- VaryDS$total+1
          PISAS_pro@ori.result[[cur_proid]] <- list(meta=scRNA_old@meta.data,filter=pro_result)
          PISAS_pro@pro.result[[cur_proid]] <- scRNA
          
          if(PISAS_pro@type == "Single"){
            PISAS_pro@integrate[["Combind"]] <- scRNA
          }
          setProgress(value = 1, message = "Finished!",detail="")
          Sys.sleep(1)
          
          step <- c(step,"new_pro")
          PISAS_pro@step <- unique(c(step,"Preprocessing"))
          PISAS_syn$PISAS_pro <<- PISAS_pro
          enable("btn_prepro")
          Info_out$out_pro <- paste0(Info_out$out_pro,"<strong style='color:green'>Success!</strong><br>")
          
        })
      },
      error = function(e) {
        Info_out$out_pro <- paste0(Info_out$out_pro,"<p style='color:red'>Stop<br>",safeError(e),"<br>","Failed!</p><br>")
        enable("btn_prepro")
        return(NULL)
        #stop(safeError(e))
      }
    )
    
  })
  
  
  ## filter info----
  validata_table <- function() {
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    cur_proid <- cur_proids()
    if(cur_proid > length(PISAS_pro@ori.result)) return()
    set_res <- PISAS_pro@ori.result[[cur_proid]]$filter
    if(is.null(set_res)) return()
    else return(set_res)
  }
  output$dsName <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    message(as.character(set_res$Dataset.name))
    return(as.character(set_res$Dataset.name))
  })
  output$species <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    return(as.character(set_res$Species))
  })
  output$Protocols <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    return(as.character(set_res$Protocol))
  })
  output$Raw_cells <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    return(as.integer(set_res$Raw.cells))
  })
  
  output$Raw_genes <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    return(as.integer(set_res$Raw.genes))
  })
  
  output$Raw_sparsity <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    return(as.character(set_res$Raw.sparsity))
  })
  
  output$Pro_cells <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    return(as.integer(set_res$Final.cells))
  })
  
  output$Pro_genes <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    return(as.integer(set_res$Final.genes))
  })
  
  output$Pro_sparsity <- renderText({
    set_res <- validata_table()
    if(is.null(set_res)) return()
    return(as.character(set_res$Final.sparsity))
  })
  
  ## QC plot----
  validata_PISAS <- function(type) {
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("Preprocessing" %in% PISAS_pro@step, "")
    )
    cur_proid <- cur_proids()
    
    if(cur_proid > length(PISAS_pro@ori.result)) return()
    if(type =="ori.meta"){
      scRNA <- PISAS_pro@ori.result[[cur_proid]]$meta
    }else if(type =="ori.filter"){
      scRNA <- PISAS_pro@ori.result[[cur_proid]]$filter
    }else if(type =="pro"){
      scRNA <- PISAS_pro@pro.result[[cur_proid]]
    }
    return(scRNA)
  }
  
  output$table_QC <- DT::renderDataTable({
    scRNA <-validata_PISAS(type="pro")
    if (is.null(scRNA))
      return()
    return(DT::datatable(scRNA@meta.data, options = list(orderClasses = TRUE,scrollX = TRUE)))
  })
  
  output$ReFC <- renderPlotly({
    data <-validata_PISAS(type="ori.meta")
    if (is.null(data))
      return()
    
    plot_ly(
      data, x = ~nFeature_RNA, y = ~nCount_RNA, color = ~filter,
      colors = c("gray","#E64B35B2"),type='scatter',  mode='markers',text = ~rownames(data)
    ) %>% 
      layout(
        legend = list(
          title="Filter"
        ),
        xaxis = list(
          title = "Number of expressed genes"
        ),
        yaxis = list(
          title = "Library sizes",
          zeroline = F
        )
      )%>%
      config(displaylogo = FALSE)
    
  })
  
  output$ReFM <- renderPlotly({
    data <-validata_PISAS(type="ori.meta")
    if (is.null(data))
      return()
    plot_ly(
      data, x = ~nFeature_RNA, y = ~percent.mt, color = ~filter,
      colors = c("gray","#E64B35B2"),type='scatter',  mode='markers',text = ~rownames(data)
    ) %>% 
      layout(
        legend = list(
          title="Filter"
        ),
        xaxis = list(
          title = "Number of expressed genes"
        ),
        yaxis = list(
          title = "Mitochondrial proportions (%)",
          zeroline = F
        )
      )%>%
      config(displaylogo = FALSE)
    
  })
  output$ReFE <- renderPlotly({
    data <-validata_PISAS(type="ori.meta")
    if (is.null(data))
      return()
    plot_ly(
      data, x = ~nFeature_RNA, y = ~percent.ERCC, color = ~filter,
      colors = c("gray","#E64B35B2"),type='scatter',  mode='markers',text = ~rownames(data)
    ) %>% 
      layout(
        legend = list(
          title="Filter"
        ),
        xaxis = list(
          title = "Number of expressed genes"
        ),
        yaxis = list(
          title = "Spike-in proportions",
          zeroline = F
        )
      )%>%
      config(displaylogo = FALSE)
    
  })
  ## normazation----
  output$QCnGeneb <- renderPlotly({
    data <-validata_PISAS(type="ori.meta")
    if (is.null(data))
      return()
    
    p <- data %>%
      plot_ly(
        x = ~orig.ident,
        y = data[,input$QC_y],
        split = ~orig.ident,
        type = 'violin',
        text = ~rownames(data),
        box = list(
          visible = T
        ),
        meanline = list(
          visible = T
        )
      ) %>% 
      layout(
        xaxis = list(
          title = ""
        ),
        yaxis = list(
          title = input$QC_y,
          zeroline = F
        ),
        title = "Before"
      )%>%
      config(displaylogo = FALSE)
    return(p)
  })
  
  output$QCnGenea <- renderPlotly({
    SCRNA <-validata_PISAS(type="pro")
    
    if (is.null(SCRNA))
      return()
    
    p <- SCRNA@meta.data %>%
      plot_ly(
        x = ~orig.ident,
        y = SCRNA@meta.data[,input$QC_y],
        split = ~orig.ident,
        type = 'violin',
        text = ~rownames(SCRNA@meta.data),
        box = list(
          visible = T
        ),
        meanline = list(
          visible = T
        )
      ) %>% 
      layout(
        xaxis = list(
          title = ""
        ),
        yaxis = list(
          title = input$QC_y,
          zeroline = F
        ),
        title="After"
      )%>%
      config(displaylogo = FALSE)
    
    return(p)
  })
  
  
  
  output$histogram_beC <- renderPlotly({
    scRNA <-validata_PISAS(type="ori.meta")
    
    if (is.null(scRNA))
      return()
    
    fig <- plot_ly(
      x = scRNA[,input$QC_y], 
      type = "histogram",
      marker = list(color = "#F39B7FB2",line = list(width = 1,
                                                    color = 'rgb(0, 0, 0)'))) %>% 
      layout(
        xaxis = list(
          title = input$QC_y
        ),
        yaxis = list(
          title = "Number of cells",
          zeroline = F
        ),
        title="Before"
      )
    
    return(fig)
  })
  
  output$histogram_afC <- renderPlotly({
    scRNA <-validata_PISAS(type="pro")
    
    if (is.null(scRNA))
      return()
    fig <- plot_ly(
      x = scRNA@meta.data[,input$QC_y],
      type = "histogram",
      marker = list(color = "#F39B7FB2",line = list(width = 1,
                                                    color = 'rgb(0, 0, 0)'))) %>% 
      layout(
        xaxis = list(
          title = input$QC_y
        ),
        yaxis = list(
          title = "Number of cells",
          zeroline = F
        ),
        title="After"
      )
    return(fig)
  })
  
  ## hvg plot----
  output$table_HVFInfo <- DT::renderDataTable({
    scRNA <-validata_PISAS(type="pro")
    if (is.null(scRNA))
      return()
    return(DT::datatable(HVFInfo(object = scRNA), options = list(orderClasses = TRUE,scrollX = TRUE)))
  })
  
  output$plot_hvg <- renderPlotly({
    scRNA <-validata_PISAS(type="pro")
    if (is.null(scRNA))
      return()
    plot <- hvgPlotly(scRNA)
    return(plot)
  })
  
  ## batch----
  
  output$color_batch <- renderUI({
    
    scRNA <-validata_PISAS(type="pro")
    if (is.null(scRNA))
      return()
    group_names <- scRNA@meta.data
    
    group_lists <- as.list(colnames(group_names)[c(-2,-3,-4,-5)])
    
    selectInput("colorBy_batch", "Color by:", group_lists)
  })
  
  observeEvent(input$btn_batch, {
    
    scRNA <-validata_PISAS(type="pro")
    
    if (is.null(scRNA))
      return()
    
    
    tryCatch(
      {
        withProgress(message = 'Plot \n', value = 0.1, {
          disable("btn_batch")
          SCE_seurat <- as.SingleCellExperiment(scRNA)
          
          if("Normalized" %in% input$dobatch){
            
            output$batch_normalization <- renderPlot({
              withProgress(message = 'Normalized dataset\n', value = 0.5, {
                
                tmp <- scater::runTSNE(
                  SCE_seurat,
                  exprs_values = "logcounts"
                )
                p <- plotTSNE(tmp,colour_by = input$colorBy_batch,size_by = "nCount_RNA",
                              run_args=list(perplexity = 10))+ggtitle("normalization")
              })
              return(p)
            })
          }
          
          if("batch-effect" %in% input$dobatch){
            
            output$batch_scaledata <- renderPlot({
              withProgress(message = 'Dataset after batch-effect removal\n', value = 0.8, {
                assay(SCE_seurat, "scaledata") <- scRNA@assays$RNA@scale.data
                
                tmp <- scater::runTSNE(
                  SCE_seurat,
                  exprs_values = "scaledata"
                )
                p <- plotTSNE(tmp,colour_by = input$colorBy_batch,size_by = "nCount_RNA",
                              run_args=list(perplexity = 10))+ggtitle(input$batchEffect)
              })
              return(p)
            })
          }
          enable("btn_batch")
          
        })
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        enable("btn_batch")
        Info_out$out_pro <- paste0(Info_out$out_pro,safeError(e),"<br>","Failed!<br>")
        
        return(NULL)
      }
    )
    
  })
  
  ## 1.5 integrate ---------------------------------------------------------------
  
  output$info_int<-renderText(Info_out$out_int)
  
  output$homoSelect <- renderUI({
    
    if(is.null(PISAS_syn$PISAS_pro)) return()
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    db <- homologene::taxData
    group_lists <- as.list(db$name_txt)
    selectInput("TaxOut", "Target species:", group_lists,selected = "Homo sapiens")
  })
  
  observeEvent(input$btn_marge, {
    tryCatch({
      message("1.5 Integrating")
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      
      PISAS_pro <- PISAS_syn$PISAS_pro
      step <- PISAS_pro@step
      infolist <- PISAS_pro@dataset.list
      
      is_deal <- sum(grepl(">Non-preprocessed",infolist$Status))
      if(is_deal > 0){
        shinyalert("Warning!", "Please preprocess all datasets first", type = "warning")
        return()
      }
      Info_out$out_int <- c("<h3>Integrating</h3>")
      
      ori.scRNAs <- list()
      inte.genes <- list()
      withProgress(message = 'Translate gene homologs across species \n', value = 0.6, {
        
        Info_out$out_int <- paste0(Info_out$out_int,"Translate gene homologs across species...<br>")
        
        # Homolog
        db <- homologene::taxData
        outTax <- db$tax_id[db$name_txt == input$TaxOut]
        
        for (i in 1:length(PISAS_pro@pro.result)) {
          if(input$comb_range == "hvg"){
            curpro <- subset(PISAS_pro@pro.result[[i]],features = VariableFeatures(object = PISAS_pro@pro.result[[i]]))
          }else{
            curpro <- PISAS_pro@pro.result[[i]]
          }
          
          if(PISAS_pro@array[[i]]@species != "Homo sapiens"){
            
            inTax <- db$tax_id[db$name_txt == PISAS_pro@array[[i]]@species]
            ori.scRNAs[[i]] <- UpdateHomoSeurat(curpro,inTax,outTax)
          }else{
            ori.scRNAs[[i]] <- curpro
          }
          
          if(infolist[infolist$ID == i & infolist$Step == "Preprocessed",]$Cells == ncol(ori.scRNAs[[i]])){
            infolist[infolist$ID == i & infolist$Step == "Integrated",]$Cells <- ncol(ori.scRNAs[[i]])
          }else{
            infolist[infolist$ID == i & infolist$Step == "Integrated",]$Cells <- paste0("<span style='color:green;font-style:italic'>",ncol(ori.scRNAs[[i]]),"</span>")
          }
          
          if(infolist[infolist$ID == i & infolist$Step == "Preprocessed",]$Genes == nrow(ori.scRNAs[[i]])){
            infolist[infolist$ID == i & infolist$Step == "Integrated",]$Genes <- nrow(ori.scRNAs[[i]])
          }else{
            infolist[infolist$ID == i & infolist$Step == "Integrated",]$Genes <- paste0("<span style='color:green;font-style:italic'>",nrow(ori.scRNAs[[i]]),"</span>")
          }
          inte.genes[[PISAS_pro@array[[i]]@name]] <- rownames(ori.scRNAs[[i]])
          
        }
      })
      PISAS_pro@integrate$inte.genes<-inte.genes
      withProgress(message = 'Integrating \n', value = 0.1, {
        disable("btn_marge")
        
        #ori.scRNAs <- PISAS_pro@pro.result
        #PISAS_pro@pro.result <- list()
        
        t1 <- proc.time()[3][[1]]
        
        if(length(ori.scRNAs) == 2){
          #merge: no scale data
          mergeObj <- merge(x = ori.scRNAs[[1]], y = ori.scRNAs[[2]])
        }else if(length(ori.scRNAs) > 2){
          mergeObj <- merge(x = ori.scRNAs[[1]], y = ori.scRNAs[2:length(ori.scRNAs)])
        }
        
        if(input$comb_func == "Intersection"){
          common_genes <- list()
          for (i in 1:length(ori.scRNAs)) {
            cur_data <- ori.scRNAs[[i]]
            common_genes[[i]] <- rownames(cur_data)
          }
          keep_genes <- Reduce(intersect, common_genes)
          mergeObj <- subset(x = mergeObj, features  = keep_genes)
          
        }
        re_scale <- TRUE
        
        if("Combind" %in% names(PISAS_pro@integrate)){
          
          if(identical(rownames(mergeObj),rownames(PISAS_pro@integrate$Combind)) && identical(colnames(mergeObj),colnames(PISAS_pro@integrate$Combind))){
            re_scale <- FALSE
          }else{
            re_scale <- TRUE
          }
        }else{
          re_scale <- TRUE
        }
        if(re_scale){
          setProgress(value = 0.2, detail = "Combine multiple datasets(~1 min)...")
          
          mergeObj <- FindVariableFeatures(mergeObj,nfeatures=2000)
          #mergeObj <- ScaleData(mergeObj,vars.to.regress = "project.PIAS", do.center = FALSE,verbose = FALSE)
          mergeObj <- ScaleData(mergeObj,vars.to.regress = "project.PIAS", do.center = FALSE,verbose = FALSE)
          
          t2 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"Combine multiple datasets...(cost ",round(t2-t1,2)," secs)<br>")
          
          setProgress(value = 0.4, detail = "UMAP dimensionality reduction(~1 min)...")
          
          mergeObj <- RunPCA(mergeObj,verbose = FALSE)
          mergeObj <- RunUMAP(mergeObj, reduction = "pca", dims = 1:30,reduction.name = "combind_umap") 
          t3 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"UMAP dimensionality reduction...(cost ",round(t3-t2,2)," secs)<br>")
          
          
          infolist <- subset(infolist,Name!="Combind")
          infolist <- rbind(infolist,
                            data.frame(ID="Integrate",Name="Combind",Groups=0,Step="Integrated",
                                       Cells=paste0("<span class='badge badge-pill badge-warning'>",dim(mergeObj)[2],"</span>"),
                                       Genes=paste0("<span class='badge badge-pill badge-warning'>",dim(mergeObj)[1],"</span>"),
                                       Status="OK",
                                       Preprocessing=NA))
          
          step <- step[step != "new_pro"]
          PISAS_pro@integrate[["Combind"]] <- mergeObj
          
        }else{
          
          mergeObj <- PISAS_pro@integrate[["Combind"]]
        }
        
        #判断方法
        t3 <- proc.time()[3][[1]]
        
        if(input$integ_func == "Seurat3"){
          costabout <- round(1.5*length(ori.scRNAs),0)
          setProgress(value = 0.5, detail = paste0("Integrate data by Seurat(~",costabout," mins)..."))
          if(length(ori.scRNAs) < 4){
            data.anchors <- FindIntegrationAnchors(object.list = ori.scRNAs, dims = 1:30)
          }else{
            data.anchors <- FindIntegrationAnchors(object.list = ori.scRNAs, dims = 1:30,reference=1)
          }
          
          data.integrated <- IntegrateData(anchorset = data.anchors, dims = 1:30)
          DefaultAssay(data.integrated) <- "integrated"
          # Run the standard workflow for visualization and clustering
          setProgress(value = 0.7, detail = "Scale data...")
          Info_out$out_int <- paste0(Info_out$out_int,"Scale data...<br>")
          data.integrated <- ScaleData(data.integrated, verbose = FALSE)
          setProgress(value = 0.8, detail = "UMAP dimensionality reduction(~1 min)...")
          t4 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"Integrate data by Seurat3...(cost ",round(t4-t3,2)," secs)<br>")
          
          data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = FALSE)
          data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:10)
          Info_out$out_int <- paste0(Info_out$out_int,"UMAP...<br>")
          t5 <- proc.time()[3][[1]]
          
          Info_out$out_int <- paste0(Info_out$out_int,"UMAP dimensionality reduction...(cost ",round(t5-t4,2)," secs)<br>")
          
          PISAS_pro@integrate[["Combind"]]@reductions$Seurat3_umap <- data.integrated@reductions$umap
          PISAS_pro@integrate[["Combind"]]@reductions$Seurat3_pca <- data.integrated@reductions$pca
          
          PISAS_pro@integrate[["Combind"]]@assays$integrated <- data.integrated@assays$integrated
          
        }else if(input$integ_func == "mnnCorrect"){
          
          setProgress(value = 0.6, detail = "Integrate data by mnnCorrect(~4 mins)...")
          
          object.list <- CheckCellNames(object.list = ori.scRNAs)
          
          features <- SelectIntegrationFeatures(object.list = object.list,features = 2000)
          object.list <- ReOrder(object.list,features)
          
          mnns <- RunFastMNN(object.list = object.list,k=20)
          
          setProgress(value = 0.8, detail = "UMAP dimension reduction...")
          t4 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"Integrate data by mnnCorrect...(cost ",round(t4-t3,2)," secs)<br>")
          
          mnns <- RunUMAP(mnns, dims = 1:ncol(mnns[["mnn"]]), reduction = "mnn")
          t5 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"UMAP dimensionality reduction...(cost ",round(t5-t4,2)," secs)<br>")
          
          PISAS_pro@integrate[["Combind"]]@reductions$mnn <- mnns@reductions$mnn
          PISAS_pro@integrate[["Combind"]]@reductions$mnn_umap <- mnns@reductions$umap
          Key(PISAS_pro@integrate[["Combind"]][["mnn_umap"]])<- "UMAP_"
          DimPlot(mnns,reduction = "umap",group.by = "project.PIAS")
          setProgress(value = 0.9, detail = "Create integrate object...")
          
        }else if(input$integ_func == "Harmony"){
          require(harmony)
          
          setProgress(value = 0.6, detail = "Integrate data by Harmony(~2 mins)...")
          
          mergeObj <- mergeObj %>% 
            RunHarmony("project.PIAS", plot_convergence = FALSE)
          
          setProgress(value = 0.8, detail = "UMAP dimension reduction...")
          t4 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"Integrate data by Harmony...(cost ",round(t4-t3,2)," secs)<br>")
          
          mergeObj <- RunUMAP(mergeObj, reduction = "harmony",dims = 1:ncol(mergeObj[["harmony"]]),reduction.name = "harmony_umap")
          
          Key(mergeObj[["harmony_umap"]])<- "UMAP_"
          #DimPlot(mergeObj,reduction = "harmony_umap",group.by = "project.PIAS")
          setProgress(value = 0.9, detail = "Create integrate object...")
          t5 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"UMAP dimensionality reduction...(cost ",round(t5-t4,2)," secs)<br>")
          
          PISAS_pro@integrate[["Combind"]]@reductions$harmony <- mergeObj@reductions$harmony
          PISAS_pro@integrate[["Combind"]]@reductions$harmony_umap <- mergeObj@reductions$harmony_umap
          
        }else if(input$integ_func == "LIGER"){
          require(liger)
          costabout <- round(1.5*length(ori.scRNAs),0)
          
          setProgress(value = 0.5, detail = paste0("Integrate data by Liger(~",costabout," mins)..."))
          
          subObj <- subset(x = mergeObj, features = VariableFeatures(object = mergeObj))
          subObj@assays$RNA@scale.data <- as.matrix(subObj@assays$RNA@data)
          
          subObj <- RunOptimizeALS(subObj, k = 20, lambda = 5, split.by = "project.PIAS")
          
          subObj <- RunQuantileNorm(subObj, split.by = "project.PIAS")
          setProgress(value = 0.8, detail = "UMAP dimension reduction...")
          t4 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"Integrate data by LIGER...(cost ~",floor(t4-t3)%/%60," mins)<br>")
          
          subObj <- RunUMAP(subObj, dims = 1:ncol(subObj[["iNMF"]]), reduction = "iNMF",reduction.name = "liger_umap")
          t5 <- proc.time()[3][[1]]
          Info_out$out_int <- paste0(Info_out$out_int,"UMAP dimensionality reduction...(cost ",round(t5-t4,2)," secs)<br>")
          
          Key(subObj[["liger_umap"]])<- "UMAP_"
          DimPlot(subObj,reduction = "liger_umap",group.by = "project.PIAS")
          PISAS_pro@integrate[["Combind"]]@reductions$iNMF <- subObj@reductions$iNMF
          PISAS_pro@integrate[["Combind"]]@reductions$iNMF_raw <- subObj@reductions$iNMF_raw
          
          PISAS_pro@integrate[["Combind"]]@reductions$liger_umap <- subObj@reductions$liger_umap
          
          setProgress(value = 0.9, detail = "Create integrate object...")
          
          Info_out$out_int <- paste0(Info_out$out_int,"Create integrate object...<br>")
          
        }
        
        PISAS_pro@step <- unique(c(step,"Integrate","dimension","test"))
        PISAS_pro@dataset.list <- infolist
        PISAS_pro@species <- input$TaxOut
        
        #PISAS_pro@pro.result <- ori.scRNAs
        PISAS_syn$PISAS_pro <<- PISAS_pro
        
        setProgress(value = 1, message = "Finish!")
        Sys.sleep(1)
      })
      
      Info_out$out_int <- paste0(Info_out$out_int,"Total cost ",floor(t5-t1)%/%60," mins ",floor(t5-t1)%%60," secs<br>")
      
      Info_out$out_int <- paste0(Info_out$out_int,"<strong style='color:green'>Success!</strong><br>")
      
      enable("btn_marge")
    },
    error = function(e) {
      Info_out$out_int <- paste0(Info_out$out_int,"<p style='color:red'>Stop<br>",safeError(e),"<br>","Failed!</p><br>")
      enable("btn_marge")
      return(NULL)
    }
    )
  })
  
  output$int_info <- renderDataTable({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("Integrate" %in% PISAS_pro@step, "")
    )
    infolist<- PISAS_pro@dataset.list
    infolist <- infolist[infolist$Step != "Raw",]
    infolist <- infolist[,c(-1,-3,-7,-8)]
    las <- nrow(infolist)
    infolist <- infolist[c(las,1:las-1),]
    
    rownames(infolist) <- NULL
    return(DT::datatable(infolist,rownames = FALSE,escape = FALSE,  options = list(orderClasses = TRUE,scrollX = TRUE)))
    
    #return(DT::datatable(NULL,rownames = FALSE,escape = FALSE,  options = list(orderClasses = TRUE,scrollX = TRUE)))
    
  })
  
  output$plot_venn <- renderPlot({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("Integrate" %in% PISAS_pro@step, "")
    )
    
    lis <- PISAS_pro@integrate$inte.genes
    
    p1 <- venn::venn(lis,zcolor = "style",box=FALSE)
    
    return(p1)
  })
  
  output$int_exp <- renderDataTable({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("Integrate" %in% PISAS_pro@step, "")
    )
    tmp <- as.data.frame(PISAS_pro@integrate$Combind@assays$RNA@scale.data)
    subset_size <- 0.1 #subsample to 10% of the data
    batch <-  PISAS_pro@integrate$Combind$project.PIAS
    subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
    subset_idr <- sample.int(n = nrow(tmp), size = floor(subset_size * nrow(tmp)), replace=FALSE)
    tmp <-  tmp[subset_idr,subset_id]
    return(DT::datatable(tmp, options = list(orderClasses = TRUE,scrollX = TRUE)))
    
  })
  
  output$int_meta <- renderDataTable({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("Integrate" %in% PISAS_pro@step, "")
    )
    tmp <- PISAS_pro@integrate$Combind@meta.data
    return(DT::datatable(tmp, options = list(orderClasses = TRUE,scrollX = TRUE)))
  })
  
  output$mergePlot1 <- renderPlotly({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("Integrate" %in% PISAS_pro@step, "")
    )
    withProgress(message = 'Ploting', value = 0.9, {
      
      ori_obj <- PISAS_pro@integrate[["Combind"]]
      
      if("combind_umap" %in% names(ori_obj)){
        reduce_var <- ori_obj@reductions$combind_umap@cell.embeddings
      }else{
        return()
      }
      groups <- ori_obj@meta.data$project.PIAS
      p <- plot_umap(reduce_var,groups,"Combind")
      
      return(p)
      
      setProgress(1)
    })
    
  })
  
  Integs <- c("Seurat3_umap","mnn_umap","harmony_umap","liger_umap")
  InFuncs <- c("Seurat3","mnnCorrect","Harmony","LIGER")
  Indts <- c("Seurat3_pca","mnn","harmony","iNMF")
  
  output$mergePlot2 <- renderPlotly({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("Integrate" %in% PISAS_pro@step, "")
    )
    withProgress(message = 'Ploting', value = 0.9, {
      
      ori_obj <- PISAS_pro@integrate[["Combind"]]
      integ_func <- input$integ_func
      n <- grep(integ_func,InFuncs)
      
      validate(
        need(Integs[n] %in% names(ori_obj@reductions), paste0("Please run ",InFuncs[n]," first!"))
      )
      
      reduce_var <- Embeddings(ori_obj,reduction = Integs[n])
      groups <- ori_obj@meta.data$project.PIAS
      p <- plot_umap(reduce_var,groups,paste0("Integrated: ",input$integ_func))
      return(p)
      
      setProgress(1)
    })
    
  })
  
  com.score <- data.frame("Lisi"=rep(0,5),"kBET"= rep(0,5),"ASW"=rep(0,5),"CH"=rep(0,5),"ARI"=rep(0,5),"Score"=rep(0,5))
  rownames(com.score) <- c("Combind","Seurat3","mnnCorrect","Harmony","LIGER")
  
  VaryScore <- reactiveValues(total=0)
  # Score ----
  output$com_score <- renderTable({
    
    if(VaryScore$total >0) message("Score Changed")
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("Integrate" %in% PISAS_pro@step, "")
    )
    ori_obj <- PISAS_pro@integrate[["Combind"]]
    
    sels <- c("Combind",InFuncs[Integs %in% names(ori_obj@reductions)])
    
    com.scores <- com.score[rownames(com.score) %in% sels,]
    com.scores[,2] <- c(1,0.97,1,0.99,0.94)
    for (i in colnames(com.scores)) {
      x <- com.scores[,i]
      if(max(x) == 0){
        com.scores[,i] <- 0
      }else{
        if(i %in% c("kBET", "ASW", "ARI","CH")){
          com.scores[,i] <- round((max(x)-x)/(max(x)-min(x)),2)
        }else{
          com.scores[,i] <- round((x-min(x))/(max(x)-min(x)),2)
        }
      }
    }
    p<-apply(com.scores,2,Percentage)
    
    if(input$set_weight == "EWM"){
      e <-apply(p,2,entropy)
      w <- eweight(e)
    }else{
      w <- as.numeric(c(input$Lisi_w,input$kBET_w,input$ASW_w,input$CH_w,input$ARI_w,100))
      w <- w/100
    }
    
    com.scores <- rbind(com.scores,Weight=w)
    com.scores[,"Score"] <- c(round(p %*% w,4)*100,1)
    
    maxval <- max(com.scores[rownames(com.scores) != "Weight",]$Score)
    best <- which(com.scores$Score == maxval)
    vals <- rownames(com.scores)[best]
    vals2 <-paste0(vals,collapse =",")
    
    output$MaxScore <- renderText({
      if(is.nan(maxval)) return("NULL")
      maxval
    })
    
    output$BestFunc <- renderText({
      if(maxval == 0 || is.nan(maxval)) return("NULL")
      vals2
    })
    
    return(com.scores)
  },striped = TRUE,hover = TRUE,align = 'c',rownames = TRUE)
  
  #Lisi----
  observeEvent(input$btn_Lisi, {
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    withProgress(message = 'LISI', value = 0.5, {
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      ori_obj <- PISAS_pro@integrate[["Combind"]]
      
      ori_var <- ori_obj@reductions$combind_umap@cell.embeddings
      
      ori_groups <- data.frame(project.PIAS = ori_obj@meta.data$project.PIAS)
      rownames(ori_groups) <- rownames(ori_obj@meta.data)
      ori_lisi <- lisi::compute_lisi(ori_var, ori_groups, c('project.PIAS'))
      df_lisi <- data.frame(mixing=ori_lisi$project.PIAS,methods=rep("Combind",nrow(ori_lisi)))
      com.score["Combind","Lisi"] <<- mean(ori_lisi$project.PIAS)
      
      p1 <- ori_var %>% 
        cbind(ori_lisi) %>% 
        dplyr::sample_frac(1L, FALSE) %>% 
        tidyr::gather(key, ori_var, project.PIAS) %>% 
        ggplot(aes(UMAP_1,  UMAP_2, color = ori_var)) + geom_point(size=1) + 
        labs(title = "Combind")+
        facet_wrap(~key)
      lisi_plot <- list(p1)
      
      for (i in 1:length(Integs)) {
        
        if(Integs[i] %in% names(ori_obj)){
          
          iLISI <- Embeddings(ori_obj,reduction = Integs[i])
          
          lisi <- lisi::compute_lisi(iLISI, ori_groups, c('project.PIAS'))
          p2 <- iLISI %>% 
            cbind(lisi) %>% 
            dplyr::sample_frac(1L, FALSE) %>% 
            tidyr::gather(key, iLISI, project.PIAS) %>% 
            ggplot(aes(UMAP_1,  UMAP_2, color = iLISI)) + geom_point(size=1) + 
            labs(title = paste0("Integrated: ",InFuncs[i]))+
            facet_wrap(~key)
          
          lisi_plot <- c(lisi_plot,list(p2))
          df_tmp <- data.frame(mixing=lisi$project.PIAS,methods=rep(InFuncs[i],nrow(lisi)))
          df_lisi <- rbind(df_lisi,df_tmp)
          com.score[InFuncs[i],"Lisi"] <<- mean(lisi$project.PIAS)
        }
      }
      
      output$lisi_plot <- renderPlot({
        p <- do.call(multiplot,c(lisi_plot,cols=2))
        return(p)
      })
      
      VaryScore$total <- VaryScore$total+1
      
      output$lisi_violin <- renderPlotly({
        
        fig <- df_lisi %>%
          plot_ly(
            x = ~methods,
            y = ~mixing,
            split = ~methods,
            type = 'violin',
            box = list(
              visible = T
            ),
            meanline = list(
              visible = T
            )
          ) %>%
          layout(
            title="LISI Results",
            yaxis = list(
              title = "iLISI batch",
              zeroline = F,
              showgrid=F
            ),
            xaxis = list(
              title = "Methods",
              zeroline = F,
              showgrid=F
            ) %>% config(displaylogo = FALSE)
          )
        
        
        return(fig)
      })
      
    })
  })
  #kBET----
  observeEvent(input$btn_kBET, {
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    withProgress(message = 'kBET(~3 mins)...', value = 0.1, {
      require(kBET)
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      ori_obj <- PISAS_pro@integrate[["Combind"]]
      
      ori_var <- ori_obj@reductions$combind_umap@cell.embeddings
      
      batch<- ori_obj$project.PIAS
      
      batch.estimate <- estimate.kBET(ori_var,batch)
      
      plot.data <- data.frame(class=rep(c('Combind'), 
                                        each=length(batch.estimate$stats$kBET.observed)), 
                              data =  c(batch.estimate$stats$kBET.observed))
      com.score["Combind","kBET"] <<- batch.estimate$summary$kBET.observed[1]
      
      for (i in 1:length(Integs)) {
        
        if(Integs[i] %in% names(ori_obj)){
          
          incProgress(0.2,  detail = InFuncs[i])
          var <-Embeddings(ori_obj,reduction = Integs[i])
          #var <- t(ori_obj@assays$integrated@scale.data)
          batch.estimate <- estimate.kBET(var,batch)
          tmpdata <- data.frame(class=rep(c(InFuncs[i]), 
                                          each=length(batch.estimate$stats$kBET.observed)), 
                                data =  c(batch.estimate$stats$kBET.observed))
          
          plot.data <- rbind(plot.data,tmpdata)
          
          com.score[InFuncs[i],"kBET"] <<-  batch.estimate$summary$kBET.observed[1]
          
        }
      }
      
      VaryScore$total <- VaryScore$total+1
      
      output$kBET_violin <- renderPlotly({
        
        fig <- plot.data %>%
          plot_ly(
            x = ~class,
            y = ~data,
            split = ~class,
            type = 'violin',
            box = list(
              visible = T
            ),
            meanline = list(
              visible = T
            )
          ) %>%
          layout(
            title="kBET Results",
            xaxis = list(
              title = "Methods",
              zeroline = F,
              showgrid=F
            ),
            yaxis = list(
              title = "Rejection rate",
              zeroline = F,
              showgrid=F
            )  %>% config(displaylogo = FALSE)
          )
        
        return(fig)
      })
      
    })
    
  })
  #ASW----
  observeEvent(input$btn_ASW, {
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    withProgress(message = 'ASW...', value = 0.5, {
      
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      ori_obj <- PISAS_pro@integrate[["Combind"]]
      
      ori_var <- ori_obj@reductions$combind_umap@cell.embeddings
      batch <- ori_obj@meta.data$project.PIAS
      
      batch.estimate <- estimate.ASW(ori_var,batch)
      
      
      plot.data <- data.frame(class=rep(c('Combind'), 
                                        each=length(batch.estimate)),
                              data =  c(batch.estimate))
      com.score["Combind","ASW"] <<- mean(batch.estimate)
      
      for (i in 1:length(Integs)) {
        
        if(Integs[i] %in% names(ori_obj)){
          evar <- Embeddings(ori_obj,reduction = Integs[i])
          
          batch.estimate <- estimate.ASW(evar,batch)
          tmpdata <- data.frame(class=rep(InFuncs[i], 
                                          each=length(batch.estimate)), 
                                data =  c(batch.estimate))
          plot.data <- rbind(plot.data,tmpdata)
          com.score[InFuncs[i],"ASW"] <<- mean(batch.estimate)
        }
      }
      
      VaryScore$total <- VaryScore$total+1
      
      output$ASW_violin <- renderPlotly({
        
        fig <- plot.data %>%
          plot_ly(
            x = ~class,
            y = ~data,
            split = ~class,
            type = 'violin',
            box = list(
              visible = T
            ),
            meanline = list(
              visible = T
            )
          ) %>%
          layout(
            title="ASW Results",
            xaxis = list(
              title = "Methods",
              zeroline = F,showgrid=F
            ),
            yaxis = list(
              title = "ASW batch",
              zeroline = F,showgrid=F
            )  %>% config(displaylogo = FALSE)
          )
        return(fig)
      })
    })
  })
  
  #CH----
  observeEvent(input$btn_CH, {
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    withProgress(message = 'Calinski-Harabasz index...', value = 0.5, {
      
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      ori_obj <- PISAS_pro@integrate[["Combind"]]
      
      ori_var <- ori_obj@reductions$combind_umap@cell.embeddings
      
      batch<- ori_obj@meta.data
      
      batch.estimate <- estimate.CH(ori_var,batch)
      
      
      plot.data <- data.frame(class=rep(c('Combind'), 
                                        each=length(batch.estimate)),
                              data =  c(batch.estimate))
      com.score["Combind","CH"] <<- batch.estimate
      
      for (i in 1:length(Integs)) {
        
        if(Integs[i] %in% names(ori_obj)){
          var <- Embeddings(ori_obj,reduction = Integs[i])
          
          batch.estimate <- estimate.CH(var,batch)
          
          tmpdata <- data.frame(class=rep(InFuncs[i], 
                                          each=length(batch.estimate)), 
                                data =  c(batch.estimate))
          plot.data <- rbind(plot.data,tmpdata)
          com.score[InFuncs[i],"CH"] <<- batch.estimate
        }
      }
      
      VaryScore$total <- VaryScore$total+1
      
      output$CH_violin <- renderPlotly({
        
        fig <- plot_ly(data = plot.data, x = ~class, y = ~data,type = "scatter",symbol = ~class,
                       marker = list(size = 10,
                                     color = 'rgba(255, 182, 193, .9)',
                                     line = list(color = 'rgba(152, 0, 0, .8)',
                                                 width = 2)))
        
        fig <- fig %>% layout(title = 'Calinski-Harabasz Results',
                              yaxis = list(title = "CH batch",zeroline = FALSE,showgrid=F),
                              xaxis = list(title = "Methods",zeroline = FALSE,showgrid=F)) %>% config(displaylogo = FALSE)
        
        return(fig)
      })
    })
  })
  
  #ARI----
  observeEvent(input$btn_ARI, {
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    withProgress(message = 'ARI...', value = 0.5, {
      
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      ori_obj <- PISAS_pro@integrate[["Combind"]]
      
      ori_var <- ori_obj@reductions$combind_umap@cell.embeddings
      
      batch<- ori_obj@meta.data
      
      batch.estimate <- estimate.ARI(ori_var,batch)
      
      
      plot.data <- data.frame(class=rep(c('Combind'), 
                                        each=length(batch.estimate)),
                              data =  c(batch.estimate))
      com.score["Combind","ARI"] <<- batch.estimate
      
      for (i in 1:length(Integs)) {
        
        if(Integs[i] %in% names(ori_obj)){
          var <- Embeddings(ori_obj,reduction = Integs[i])
          
          batch.estimate <- estimate.ARI(var,batch)
          
          tmpdata <- data.frame(class=rep(InFuncs[i], 
                                          each=length(batch.estimate)), 
                                data =  c(batch.estimate))
          plot.data <- rbind(plot.data,tmpdata)
          com.score[InFuncs[i],"ARI"] <<- batch.estimate
        }
      }
      
      VaryScore$total <- VaryScore$total+1
      
      output$ARI_violin <- renderPlotly({
        
        fig <- plot_ly(data = plot.data, x = ~class, y = ~data,type = "scatter",symbol = ~class,
                       marker = list(size = 10,
                                     color = 'rgba(255, 182, 193, .9)',
                                     line = list(color = 'rgba(152, 0, 0, .8)',
                                                 width = 2)))
        
        fig <- fig %>% layout(title = 'ARI Results',
                              yaxis = list(title = "Average ARI batch",zeroline = FALSE,showgrid=F),
                              xaxis = list(title = "Methods",zeroline = FALSE,showgrid=F)) %>% config(displaylogo = FALSE)
        
        return(fig)
      })
    })
  })
  
  ##
  
  ## 2 dimension ---------------------------------------------------------------
  #   define the selector of metadata
  
  output$info_dim<-renderText(Info_out$out_dim)
  
  output$choose_intfunc2 <- renderUI({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    if(PISAS_pro@type=="Single"){
      group_lists <- list("Single")
    }else{
      group_lists <- list()
      for (i in 1:length(Integs)) {
        if(Integs[i] %in% names(PISAS_pro@integrate$Combind)){
          group_lists <- c(group_lists,InFuncs[i])
        }
      }
    }
    selectInput("cho_mes2", "Integrated methods:",group_lists,selected = sel_intfunc)
  })
  
  output$campSelector2 <- renderUI({
    
    if(is.null(PISAS_syn$PISAS_pro)) return()
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    validate(
      need("Preprocessing" %in% PISAS_pro@step, "")
    )
    group_names <- PISAS_pro@integrate[["Combind"]]@meta.data
    group_lists <- as.list(colnames(group_names))
    selectInput("groupBy", "Color by:", group_lists,selected = "nCount_RNA")
  })
  #
  observeEvent(input$btn_dim, {
    
    tryCatch({
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      disable("btn_dim")
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      validate(
        need("Preprocessing" %in% PISAS_pro@step, "Please run preprocessing first....")
      )
      
      scRNA <- PISAS_pro@integrate[["Combind"]]
      step <- PISAS_pro@step
      withProgress(message = 'Calculating the data', value = 0.1, {
        setProgress(0.4)
        withProgress(message = 'Dimension reduction:PCA', detail = "This may take a while...", value = 7, {
          pcs.use <- 1:input$pcs.compute
          Info_out$out_dim <- c("<h3>Dimension reduction</h3>")
          
          type <- PISAS_pro@type
          
          if(type == "Multiple"){
            n <- grep(input$cho_mes2,InFuncs)
            nsid <- Indts[n]
            redt <- Embeddings(scRNA,reduction = nsid)
            if(nsid == "Seurat3_pca") nsid = "Seurat3"
            #is tsne
            if(input$dimfunc == "tSNE"){
              Info_out$out_dim <- paste0(Info_out$out_dim,"Dimension reduction by TSNE...<br>")
              
              redt1 <- RunTSNE(object = redt, dims.use = pcs.use, do.fast = TRUE,
                               perplexity = input$perplexity, max_iter = input$max_iter,dim.embed=2)
              
              reid <- paste0(nsid,"_tsne")
              scRNA@reductions[[reid]] <- redt1
              
              redt2 <- RunTSNE(object = redt, dims.use = pcs.use, do.fast = TRUE,
                               perplexity = input$perplexity, max_iter = input$max_iter,dim.embed=3)
              
              reid <- paste0(nsid,"_tsne_3d")
              scRNA@reductions[[reid]] <- redt2
              
            }else if(input$dimfunc == "UMAP"){
              Info_out$out_dim <- paste0(Info_out$out_dim,"Dimension reduction by UMAP...<br>")
              
              # redt1 <- RunUMAP(redt, dims = pcs.use,n.components = 2)
              # reid <- paste0(nsid,"_umap")
              # scRNA@reductions[[reid]] <- redt1
              # 
              redt2 <- RunUMAP(redt, dims = pcs.use,n.components = 3,reduction.name = "umap_3d")
              reid <- paste0(nsid,"_umap_3d")
              scRNA@reductions[[reid]] <- redt2
            }
          }else{
            
            scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA),verbose = F)
            
            Info_out$out_dim <- paste0(Info_out$out_dim,"Dimension reduction by PCA...<br>")
            
            #is tsne
            if(input$dimfunc == "tSNE"){
              Info_out$out_dim <- paste0(Info_out$out_dim,"Dimension reduction by TSNE...<br>")
              
              scRNA <- RunTSNE(object = scRNA, dims.use = pcs.use, do.fast = TRUE,
                               perplexity = input$perplexity, max_iter = input$max_iter,dim.embed=2)
              
              scRNA <- RunTSNE(object = scRNA, dims.use = pcs.use, do.fast = TRUE,
                               perplexity = input$perplexity, max_iter = input$max_iter,dim.embed=3,reduction.name = "tsne_3d")
              
            }else if(input$dimfunc == "UMAP"){
              Info_out$out_dim <- paste0(Info_out$out_dim,"Dimension reduction by UMAP...<br>")
              scRNA <- RunUMAP(scRNA, dims = pcs.use,n.components = 2)
              scRNA <- RunUMAP(scRNA, dims = pcs.use,n.components = 3,reduction.name = "umap_3d")
            }
          }
          PISAS_pro@integrate[["Combind"]] <- scRNA
          PISAS_pro@step <- unique(c(step,"dimension"))
          PISAS_syn$PISAS_pro <<- PISAS_pro
        })
        setProgress(value = 1, message = "Finish!")
        Sys.sleep(1)
      })
      enable("btn_dim")
      Info_out$out_dim <- paste0(Info_out$out_dim,"Finished!")
      
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      enable("btn_dim")
      
      Info_out$out_dim <- paste0(Info_out$out_dim,safeError(e),"<br>Failed!")
      
      return(NULL)
    })
  })
  
  output$Dimplot2D <- renderPlotly({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("dimension" %in% PISAS_pro@step, "Please run dimension reduction first....")
    )
    validate(
      need(!is.null(input$groupBy), "Loading....")
    )
    scRNA <- PISAS_pro@integrate[["Combind"]]
    
    type <- PISAS_pro@type
    
    if(type == "Multiple"){
      n <- grep(input$cho_mes2,InFuncs)
      nsid <- Indts[n]
      
      if(input$dimfunc == "PCA"){
        recname <- nsid
      }else{
        if(nsid == "Seurat3_pca") nsid = "Seurat3"
        recname <- paste0(nsid,"_",tolower(input$dimfunc))
      }
      
      validate(
        need(recname %in% names(scRNA@reductions), paste0("Please run ",input$dimfunc," first...."))
      )
      
      message("recname:",recname)
      reduce_var <- scRNA[[recname]]@cell.embeddings[,1:2]
      
    }else{
      recname <- tolower(input$dimfunc)
      validate(
        need(recname %in% names(scRNA@reductions), paste0("Please run ",toupper(recname)," first...."))
      )
      reduce_var <- scRNA[[recname]]@cell.embeddings[,1:2]
      
    }
    
    groupBy <- input$groupBy
    orders <- which(colnames(scRNA@meta.data) == groupBy)
    groups <- scRNA@meta.data[[orders]]
    
    fig <- plot_umap(reduce_var,groups,paste0(input$dimfunc,"_2D"),colorbar_name= groupBy)
    
    return(fig)
  })
  
  output$Dimplot3D <- renderPlotly({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("dimension" %in% PISAS_pro@step, "")
    )
    validate(
      need(!is.null(input$groupBy), "Loading....")
    )
    
    scRNA <- PISAS_pro@integrate[["Combind"]]
    type <- PISAS_pro@type
    if(type == "Multiple"){
      n <- grep(input$cho_mes2,InFuncs)
      nsid <- Indts[n]
      if(input$dimfunc == "PCA"){
        recname <- nsid
      }else{
        if(nsid == "Seurat3_pca") nsid = "Seurat3"
        recname <- paste0(nsid,"_",tolower(input$dimfunc),"_3d")
      }
      
      validate(
        need(recname %in% names(scRNA@reductions), paste0("Please run ",input$dimfunc," 3D first...."))
      )
      
      message("recname:",recname)
      
      reduce_var <- scRNA[[recname]]@cell.embeddings[,1:3]
      
    }else{
      recname <- paste0(tolower(input$dimfunc),"_3d")
      if(input$dimfunc == "PCA"){
        recname <- "pca"
      }
      validate(
        need(recname %in% names(scRNA@reductions), paste0("Please run ",toupper(recname)," first...."))
      )
      reduce_var <- scRNA[[recname]]@cell.embeddings[,1:3]
      
    }
    
    groupBy <- input$groupBy
    orders <- which(colnames(scRNA@meta.data) == groupBy)
    groups <- scRNA@meta.data[[orders]]
    
    fig <- plotly_3d(reduce_var,groups,paste0(input$dimfunc,"_3D"),colorbar_name= groupBy)
    
    return(fig)
  })
  
  ## 
  
  ## 3 cluster ---------------------------------------------------------------
  
  output$info_clu<-renderText(Info_out$out_clu)
  sel_intfunc <- "Seurat3"
  output$choose_intfunc <- renderUI({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    if(PISAS_pro@type=="Single"){
      group_lists <- list("Single")
    }else{
      group_lists <- list()
      for (i in 1:length(Integs)) {
        if(Integs[i] %in% names(PISAS_pro@integrate$Combind)){
          group_lists <- c(group_lists,InFuncs[i])
        }
      }
    }
    selectInput("cho_mes", "Integrated methods:",group_lists,selected = sel_intfunc)
  })
  
  
  
  sel_g3 <- ""
  output$campSelector3 <- renderUI({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    group_lists <- list()
    meta <- PISAS_pro@integrate$Combind@meta.data
    j <- grep("_clusters", colnames(meta))
    group_lists <- c(group_lists,as.list(colnames(meta)[j]))
    
    selectInput("groupBy3", "Group by:",group_lists,selected = sel_g3)
  })
  
  observeEvent(input$btn_cluster, {
    
    tryCatch({
      
      disable("btn_cluster")
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      type <- PISAS_pro@type
      
      validate(
        need("dimension" %in% PISAS_pro@step, "")
      )
      
      step <- PISAS_pro@step
      
      sel_intfunc <<- input$cho_mes
      
      withProgress(message = 'Clustering...', value = 0.1, {
        
        setProgress(0.4)
        Info_out$out_clu <- c("<h3>Cluster</h3>")
        
        use.phenograph <- FALSE
        scRNA <- PISAS_pro@integrate[["Combind"]]
        
        if(type == "Single"){
          Info_out$out_clu <- paste0(Info_out$out_clu,paste0(input$dimchos,"...<br>"))
          
          if(!(tolower(input$dimchos) %in% names(scRNA@reductions))){
            redps <- paste0("Please return to the previous tab for ",input$dimchos," dimensionality reduction")
            shinyalert("Sorry!",redps, type = "warning")
            return(NULL)
          }
          
          if(input$dimchos == "PCA"){
            input.dat <- scRNA@reductions[["pca"]]@cell.embeddings
            
          }else if(input$dimchos == "tSNE"){
            input.dat <- scRNA@reductions[["tsne"]]@cell.embeddings
            
          }else if(input$dimchos == "UMAP"){
            input.dat <- scRNA@reductions[["umap"]]@cell.embeddings
            
          }
          
        }else{
          Info_out$out_clu <- paste0(Info_out$out_clu,paste0(input$cho_mes,".",input$dimchos,"...<br>"))
          
          n <- grep(input$cho_mes,InFuncs)
          nsid <- Indts[n]
          if(input$dimchos == "PCA"){
            recname <- nsid
          }else{
            if(nsid == "Seurat3_pca") nsid = "Seurat3"
            recname <- paste0(nsid,"_",tolower(input$dimchos))
          }
          message("recname:",recname)
          
          if(recname %in% names(scRNA@reductions)){
            input.dat <- Embeddings(scRNA,reduction = recname)
          }else{
            Info_out$out_clu <- paste0(Info_out$out_clu,"Please run ",recname," first...<br>")
          }
          
        }
        
        Info_out$out_clu <- paste0(Info_out$out_clu,paste0(input$clusterfunc,"...<br>"))
        
        col.name <- character()
        if(input$clusterfunc == "hierarchical"){
          
          dat <- apply(t(input.dat), 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
          
          dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
          
          hc <- hclust(dd, method = input$method.hc)
          
          result <- cutree(hc, k = input$k.hc)
          scRNA@active.ident <- factor(result)
          col.name <- paste0(input$cho_mes,".",input$dimchos,'.hierarchical_clusters')
          PISAS_pro@integrate$Combind@meta.data[,col.name] <- factor(result)
          
        }else if(input$clusterfunc == "kmeans"){
          cluster <- as.character(kmeans(input.dat, centers = input$centers.ks)$clust)
          scRNA@active.ident <- factor(cluster)
          col.name <-  paste0(input$cho_mes,".",input$dimchos,'.kmeans_clusters')
          PISAS_pro@integrate$Combind@meta.data[,col.name] <- factor(cluster)
          
        }else{
          
          message(input$cho_mes)
          
          if(input$cho_mes == "Seurat3"){
            DefaultAssay(scRNA) <- "integrated"
          }
          n <- grep(input$cho_mes,InFuncs)
          if(input$cho_mes == "Single"){
            scRNA <- FindNeighbors(scRNA, reduction = "pca")
          }else{
            scRNA <- FindNeighbors(scRNA, reduction = Indts[n])
          }
          
          if(input$clusterfunc == "Louvain"){
            algorithm = 1
          }else if(input$clusterfunc == "SLM"){
            algorithm = 3
          }else if(input$clusterfunc =="Leiden"){
            library(leiden)
            algorithm = 4
          }
          
          #scRNA <- FindClusters(object = scRNA, resolution = 0.1,algorithm = 1, verbose = FALSE)
          #DimPlot(scRNA,reduction = "liger_umap")+labs(title="LIGER+Louvain")
          scRNA <- FindClusters(object = scRNA, resolution = input$resolution,algorithm = algorithm, verbose = FALSE)
          col.name <-  paste0(input$cho_mes,".",input$dimchos,".",input$clusterfunc,".res",input$resolution,"_clusters")
          PISAS_pro@integrate$Combind@meta.data[,col.name] <- scRNA@meta.data$seurat_clusters
        }
        
        #PISAS_pro@integrate$Combind <- scRNA
        
        sel_g3 <<- col.name
        
        PISAS_pro@step <- unique(c(step,"cluster"))
        
        PISAS_syn$PISAS_pro <<- PISAS_pro
        
        setProgress(1,message = "Finished!")
        Sys.sleep(1)
      })
      enable("btn_cluster")
      Info_out$out_clu <- paste0(Info_out$out_clu,c("Finished!<br>"))
      
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      enable("btn_cluster")
      Info_out$out_clu <- paste0(Info_out$out_clu,safeError(e),"<br>Failed!<br>")
      return(NULL)
      #stop(safeError(e))
      
    })
    
  })
  
  
  observeEvent(input$btn_k, {
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    type <- PISAS_pro@type
    
    validate(
      need("dimension" %in% PISAS_pro@step, "")
    )
    
    disable("btn_k")
    scRNA <- PISAS_pro@integrate[["Combind"]]
    
    if(type == "Single"){
      if(input$dimchos == "PCA"){
        input.dat <- scRNA@reductions[["pca"]]@cell.embeddings
        
      }else if(input$dimchos == "tSNE"){
        input.dat <- scRNA@reductions[["tsne"]]@cell.embeddings
        
      }else if(input$dimchos == "UMAP"){
        input.dat <- scRNA@reductions[["umap"]]@cell.embeddings
      }
    }else{
      
      Info_out$out_clu <- paste0(Info_out$out_clu,paste0(input$cho_mes,".",input$dimchos,"...<br>"))
      
      n <- grep(input$cho_mes,InFuncs)
      nsid <- Indts[n]
      if(input$dimchos == "PCA"){
        recname <- nsid
      }else{
        if(nsid == "Seurat3_pca") nsid = "Seurat3"
        recname <- paste0(nsid,"_",tolower(input$dimchos))
      }
      
      if(recname %in% names(scRNA@reductions)){
        input.dat <- Embeddings(scRNA,reduction = recname)
      }else{
        Info_out$out_clu <- paste0(Info_out$out_clu,"Please run ",recname," first...<br>")
      }
      
    }
    
    library(factoextra)
    library(cluster)
    output$elbowPlot <- renderPlot({
      withProgress(message = 'Determining Optimal Clusters...', value = 0.4, {
        fviz_nbclust(input.dat, kmeans, method = "wss")
      })
    })
    
    output$silPlot <- renderPlot({
      withProgress(message = 'Determining Optimal Clusters...', value = 0.8, {
        fviz_nbclust(input.dat, kmeans, method = "silhouette")
      })
    })
    
    enable("btn_k")
  })
  
  output$clusterPlot <- renderPlotly({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("cluster" %in% PISAS_pro@step, "Please click \"Run\" to cluster samples first....")
    )
    
    type <- PISAS_pro@type
    groupBy <- input$groupBy3
    message("groupBy",groupBy=="")
    validate(
      need(groupBy!="", "")
    )
    mes <- unlist(strsplit(groupBy, split = "[.]"))[1]
    scRNA <- PISAS_pro@integrate[["Combind"]]
    if(type == "Single"){
      if(input$dimchos == "PCA"){
        reduce_var <- scRNA@reductions[["pca"]]@cell.embeddings
        reduce_umap <- scRNA@reductions[["pca"]]@cell.embeddings
        
      }else if(input$dimchos == "tSNE"){
        reduce_var <- scRNA@reductions[["tsne"]]@cell.embeddings
        reduce_umap <- scRNA@reductions[["tsne"]]@cell.embeddings
        
      }else if(input$dimchos == "UMAP"){
        reduce_var <- scRNA@reductions[["umap"]]@cell.embeddings
        reduce_umap <- scRNA@reductions[["umap"]]@cell.embeddings
        
      }
      
    }else{
      n <- grep(mes,InFuncs)
      
      nsid <- Indts[n]
      if(input$dimchos == "PCA"){
        recname <- nsid
      }else{
        if(nsid == "Seurat3_pca") nsid = "Seurat3"
        recname <- paste0(nsid,"_",tolower(input$dimchos))
      }
      
      validate(
        need(recname %in% names(scRNA@reductions), paste0("Please run ",mes,"_",input$dimchos," first!"))
      )
      
      reduce_var <- Embeddings(scRNA,reduction = Indts[n])
      reduce_umap <- Embeddings(scRNA,reduction = recname)
      
    }
    
    #orders <- which(colnames(scRNA@meta.data) == groupBy)
    groups <- scRNA@meta.data[[groupBy]]
    
    # len_batch <- length(levels(factor(groups)))
    # validate(
    #   need(len_batch>1, "Please select another group!")
    # )
    
    ASW.groups <- estimate.ASW(reduce_var,groups)
    
    output$ASWPlot <- renderPlotly({
      validate(
        need("cluster" %in% PISAS_pro@step, "Please click \"Run\" to cluster samples first....")
      )
      plot.data <- data.frame(class=rep(mes, each=length(ASW.groups)),
                              data = c(ASW.groups))
      fig <- plot.data %>%
        plot_ly(
          x = ~class,
          y = ~data,
          split = ~class,
          type = 'violin',
          box = list(
            visible = T
          ),
          meanline = list(
            visible = T
          )
        ) %>%
        layout(
          title="",
          xaxis = list(
            title = "Methods",
            zeroline = F,showgrid=F
          ),
          yaxis = list(
            title = "ASW cluster",
            zeroline = F,showgrid=F
          )  %>% config(displaylogo = FALSE)
        )
      return(fig)
    })
    
    p <- plot_umap(reduce_umap,groups,groupBy)
    return(p)
  })
  
  cluster_results <- data.frame()
  
  output$table_cluster <- DT::renderDataTable({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    type <- PISAS_pro@type
    validate(
      need("cluster" %in% PISAS_pro@step, "")
    )
    meta <- PISAS_pro@integrate$Combind@meta.data
    j <- grep("_clusters", colnames(meta))
    cluster_results <<- meta
    return(DT::datatable(meta[,j], options = list(orderClasses = TRUE,scrollX = TRUE)))
  })
  
  output$download_cluster <- downloadHandler(
    filename = "cluster_result.txt",
    content = function(file) {
      if(nrow(cluster_results) == 0) return()
      new_meta <- cluster_results
      write.table (new_meta, file, sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)
    }
  )
  
  ## 
  
  ## 4 differential ----------------------------------------------------------
  
  output$info_diff<-renderText(Info_out$out_diff)
  
  output$choose_cluster <- renderUI({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    group_lists <- list()
    meta <- PISAS_pro@integrate$Combind@meta.data
    j <- grep("_clusters", colnames(meta))
    group_lists <- c(group_lists,as.list(colnames(meta)[j]))
    if("annotations" %in% colnames(meta)) group_lists <- c(group_lists,"annotations")
    selectInput("cho_clu", "Cluster results:",group_lists)
  })
  
  output$choose_set1 <- renderUI({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    type <- PISAS_pro@type
    meta <- PISAS_pro@integrate$Combind@meta.data
    s <- levels(factor(meta[,input$cho_clu]))
    group_lists <- as.list(s)
    selectInput("set1", "experimental group:",group_lists)
  })
  
  output$choose_set2 <- renderUI({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    type <- PISAS_pro@type
    meta <- PISAS_pro@integrate$Combind@meta.data
    s <- levels(factor(meta[,input$cho_clu]))
    group_lists <- as.list(s)
    selectInput("set2", "control group:",group_lists,multiple = T)
  })
  
  observeEvent(input$btn_diff, {
    tryCatch({
      disable("btn_diff")
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      type <- PISAS_pro@type
      
      validate(
        need("cluster" %in% PISAS_pro@step, "")
      )
      step <- PISAS_pro@step
      scRNA <- PISAS_pro@integrate$Combind
      #Idents(scRNA) <- factor(scRNA@meta.data[,input$cho_clu])
      Idents(scRNA) <- input$cho_clu
      
      #names(PISAS_pro@integrate[["Combind"]]@active.ident) <- rownames(scRNA@meta.data)
      
      step <- PISAS_pro@step
      withProgress(message = 'Differential expression analysis', detail = "This may take a while...", value = 0.3, {
        Info_out$out_diff <- c("<h3>Differential expression analysis</h3>")
        # judge the genesets
        if(input$deSet == "all"){
          Info_out$out_diff <- paste0(Info_out$out_diff,c("Finds markers(differentially expressed genes) for each of cluster...<br>"))
          
          scRNA.markers <- FindAllMarkers(object = scRNA, only.pos = input$only_pos,
                                          logfc.threshold = input$logfc.threshold,min.pct=input$min.pct,test.use = input$DEfunc)
          
        }else{
          
          if(input$set1 %in% input$set2){
            shinyalert("Sorry!", "Please do not choose gene sets that have intersection!", type = "warning")
            return(NULL)
          }else{
            
            Info_out$out_diff <- paste0(Info_out$out_diff,c("Finds markers(differentially expressed genes) for a cluster...<br>"))
            
            scRNA.markers <- FindMarkers(scRNA, ident.1 = as.numeric(input$set1), ident.2 = input$set2,only.pos = input$only_pos,
                                         logfc.threshold = input$logfc.threshold,min.pct=input$min.pct,test.use = input$DEfunc)
            scRNA.markers$cluster <- input$set1
            scRNA.markers$gene <- rownames(scRNA.markers)
          }
        }
        PISAS_pro@integrate[["Markers"]] <- scRNA.markers
        PISAS_pro@integrate[["DE.info"]] <- c(input$cho_clu)
        
        if(PISAS_pro@type == "Single"){
          PISAS_pro@integrate[["assay.map"]] <- tolower(unlist(strsplit(PISAS_pro@integrate[["DE.info"]],"[.]"))[2])
        }else{
          n <- grep(unlist(strsplit(input$cho_clu, split = "[.]"))[1],InFuncs)
          PISAS_pro@integrate[["assay.map"]] <-  Integs[n]
        }
        message("assay.map: ",PISAS_pro@integrate[["assay.map"]])
        PISAS_pro@integrate$Combind <- scRNA
        PISAS_pro@step <- c(step,"DE")
        PISAS_syn$PISAS_pro <<- PISAS_pro
        message("DE... finish!")
        
        updateSelectInput(session, "cho_clu",
                          selected = input$cho_clu)
        updateSelectInput(session, "set1",
                          selected = input$set1)
        updateSelectInput(session, "set2",
                          selected = input$set2)
        
        setProgress(value = 1, message = "Finished!",detail="")
        Sys.sleep(1)
      })
      enable("btn_diff")
      Info_out$out_diff <- paste0(Info_out$out_diff,c("Finished!<br>"))
      
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      enable("btn_diff")
      Info_out$out_diff <- paste0(Info_out$out_diff,safeError(e),"<br>Failed!<br>")
      return(NULL)
    })
    
  })
  
  marker.genes <- reactive({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("DE" %in% PISAS_pro@step, "")
    )
    scRNA.markers <- PISAS_pro@integrate[["Markers"]]
    topMarkers <- NULL
    if(input$DEfunc %in% c('wilcox', 'bimod', 't', "LR")){
      if("avg_log2FC" %in% colnames(scRNA.markers)){
        topMarkers <- scRNA.markers %>% group_by(cluster) %>% top_n(input$topNum,wt= avg_log2FC)
      }else{
        topMarkers <- scRNA.markers %>% group_by(cluster) %>% top_n(input$topNum,wt= avg_logFC)
      }
      
    }else if(input$DEfunc == "roc"){
      
      if("power" %in% colnames(scRNA.markers)){
        topMarkers <- scRNA.markers %>% group_by(cluster) %>% top_n(input$topNum,wt= power)
      }
      
    }
    return(topMarkers)
  })
  
  output$gene_Lists <- renderUI({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    if(is.null(marker.genes())){
      return(NULL)
    }
    topMarkers <- marker.genes()
    subMarkers <- topMarkers[,c("cluster","gene")]
    group_lists <- split(paste(subMarkers$cluster,":",subMarkers$gene), subMarkers$cluster)
    selectInput("selgenes", "Interesting genes:",group_lists,multiple = T)
  })
  
  output$diffPlot <- renderPlot({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("DE" %in% PISAS_pro@step, "")
    )
    type <- PISAS_pro@type
    scRNA <- PISAS_pro@integrate$Combind
    
    assay.map <- PISAS_pro@integrate[["assay.map"]]
    message(assay.map)
    validate(
      need(!is.null(input$selgenes), "** Please select the genes of interest to display!")
    )
    markers.to.plot <- input$selgenes
    #markers.to.plot <- gsub("[0-9]+ : ","",markers.to.plot)
    markers.to.plot <- gsub(".*: ","",markers.to.plot)
    
    if(!is.null(markers.to.plot)){
      if(input$deplots == "Feature"){
        p <- FeaturePlot(scRNA, features = markers.to.plot)
        
        if(length(assay.map)!=0){
          for (i in 1:length(markers.to.plot)) {
            p[[i]]$data[,1:2] <- Embeddings(scRNA,reduction = assay.map) 
          }
        }
      }else if(input$deplots == "Ridge"){
        p <- RidgePlot(scRNA, features = markers.to.plot, ncol = 2)
        
      }else if(input$deplots == "Dot"){
        p <- DotPlot(scRNA, features = markers.to.plot) + RotatedAxis()
        
      }else if(input$deplots == "Violin"){
        p <- VlnPlot(scRNA, features = markers.to.plot,ncol=2,pt.size = 0)
        
      }else if(input$deplots == "Heatmap"){
        p <- DoHeatmap(scRNA, features = markers.to.plot, size = 3)
        
      }
    }else{
      return(NULL)
    }
    
    return(p)
  })
  
  
  output$table_DE <- DT::renderDataTable({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("DE" %in% PISAS_pro@step, "Please run cell differential expression analysis reduction first....")
    )
    
    scRNA.markers <- PISAS_pro@integrate[["Markers"]]
    
    return(DT::datatable(scRNA.markers, options = list(orderClasses = TRUE,scrollX = TRUE)))
    
  })
  
  output$download_DE <- downloadHandler(
    filename = "DE_result.txt",
    content = function(file) {
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      validate(
        need("DE" %in% PISAS_pro@step, "Please run cell differential expression analysis reduction first....")
      )
      
      scRNA.markers <- PISAS_pro@integrate[["Markers"]]
      write.table (scRNA.markers, file, sep ="\t", row.names =TRUE, col.names =TRUE, quote =FALSE)
    }
  )
  
  
  
  ## 5 enrichment analysis ------------------------------------------------------------
  output$info_rich<-renderText(Info_out$out_rich)
  
  output$clusterSelector1 <- renderUI({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("DE" %in% PISAS_pro@step, "Please run cell differential expression analysis reduction first....")
    )
    scRNA.markers <- PISAS_pro@integrate[["Markers"]]
    
    cluster_level <- levels(scRNA.markers$cluster)
    selectInput("gcluster1", "Gene list:", as.list(cluster_level))
  })
  
  observeEvent(input$btn_go, {
    tryCatch({
      disable("btn_go")
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      
      validate(
        need("DE" %in% PISAS_pro@step, "Please run cell differential expression analysis reduction first....")
      )
      scRNA.markers <- PISAS_pro@integrate[["Markers"]]
      step <- PISAS_pro@step
      withProgress(message = 'Enrichment analysis...', detail = "This may take about 5 minutes...",value = 0.3, {
        
        Info_out$out_rich <- c("<h3>Enrichment analysis</h3>")
        
        geneList <- scRNA.markers[(scRNA.markers$p_val_adj < 0.05) & (scRNA.markers$cluster == input$gcluster1),]
        
        if("avg_logFC" %in% colnames(geneList)){
          geneSet <- geneList$avg_logFC
        }else if("avg_log2FC" %in% colnames(geneList)){
          geneSet <- geneList$avg_log2FC
        }else if("power" %in% colnames(geneList)){
          geneSet <- geneList$power
        }
        
        geneSet <- sort(geneSet, decreasing = TRUE)
        
        species <- PISAS_pro@species
        Org.Db <- LoadOrgDB(species)
        if(species == "Xenopus (Silurana) tropicalis") species <- "Xenopus tropicalis"
        
        organisms <- search_kegg_organism(species, by='scientific_name')
        organisms <- organisms$kegg_code[1]
        
        ENTREZIDs <- mapIds(Org.Db, keys=geneList$gene,  keytype= "SYMBOL", column="ENTREZID")
        names(geneSet) <- ENTREZIDs
        
        gores <- NULL
        keggres <- NULL
        
        if("GO" %in% input$Enmeas){
          if(input$enrichFc == "ORA"){
            Info_out$out_rich <- paste0(Info_out$out_rich,"GO analysis by ORA...<br>")
            
            gores <- enrichGO(names(geneSet),OrgDb = Org.Db, ont = input$ont, keyType = "ENTREZID",pAdjustMethod = input$pAdjustMethod, pvalueCutoff = input$pvalueCutoff, qvalueCutoff = input$qvalueCutoff)
            #gores <- enrichGO(geneList$gene,OrgDb = org.Hs.eg.db, ont = "BP", keyType = "ENSEMBL",pAdjustMethod = "fdr", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
            
          }else{
            Info_out$out_rich <- paste0(Info_out$out_rich,"GO analysis by FCS...<br>")
            gores <- gseGO(
              geneList  = geneSet,
              OrgDb  = Org.Db,
              ont  = input$ont,
              keyType = "ENTREZID",
              pvalueCutoff = input$pvalueCutoff)
          }
          
          gores <- setReadable(gores, OrgDb = Org.Db)
          step <- unique(c(step,"GO"))
          PISAS_pro@integrate[["GO"]] <- gores
        }
        
        if("KEGG" %in% input$Enmeas){
          
          if(input$enrichFc == "ORA"){
            Info_out$out_rich <- paste0(Info_out$out_rich,"KEGG pathway analysis by ORA...<br>")
            
            keggres <- enrichKEGG(
              gene = names(geneSet),
              keyType = "kegg",
              organism  = organisms,
              pvalueCutoff  = input$pvalueCutoff,
              pAdjustMethod  = input$pAdjustMethod,
              qvalueCutoff  = input$qvalueCutoff
            )
            
          }else{
            Info_out$out_rich <- paste0(Info_out$out_rich,"KEGG pathway analysis by FCS...<br>")
            
            keggres <- gseKEGG(
              geneList  = geneSet,
              keyType  = 'kegg',
              organism = organisms,
              pvalueCutoff = input$pvalueCutoff
            )
          }
          step <- unique(c(step,"KEGG"))
          PISAS_pro@integrate[["KEGG"]] <- keggres
        }
        message(step)
        PISAS_pro@step <- step
        PISAS_syn$PISAS_pro <<- PISAS_pro
        setProgress(1)
      })
      enable("btn_go")
      Info_out$out_rich <- paste0(Info_out$out_rich,"Finished<br>")
      
    },
    error = function(e) {
      # return a safeError if a parsing error occurs
      enable("btn_go")
      Info_out$out_rich <- paste0(Info_out$out_rich,safeError(e),"<br>Failed<br>")
      return(NULL)
    })
  })
  
  #GO绘图
  output$goPlot <- renderPlot({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("GO" %in% PISAS_pro@step, "Please run cell GO enrichment analysis reduction first....")
    )
    gores <- PISAS_pro@integrate[["GO"]]
    
    validate(
      need(nrow(gores@result) > 0, "No gene can be mapped or No term enriched under specific pvalueCutoff! Please increase the value of p-value cutoff.")
    )
    
    if(input$goplots == "dotplot"){
      
      goplot <- enrichplot::dotplot(gores,showCategory=input$showCategory,color=input$gocolorBy,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free")
      
    }else if(input$goplots == "cnetplot"){
      ## remove redundent GO terms
      #gores2 <- simplify(gores) #take too long 
      goplot <-  enrichplot::cnetplot(gores,showCategory = input$showCategory,categorySize=input$gocolorBy,circular = TRUE, colorEdge = TRUE)
      #upsetplot(gores)
      
    }else if(input$goplots == "emapplot"){
      goplot <- enrichplot::emapplot(gores, showCategory = input$showCategory, color = input$gocolorBy)
      
    }else if(input$goplots == "heatplot"){
      if("geneList" %in%  slotNames(gores)){
        goplot <- enrichplot::heatplot(gores, showCategory = input$showCategory, foldChange = gores@geneList)+ggplot2::coord_flip()
      }else{
        goplot <- enrichplot::heatplot(gores, showCategory = input$showCategory)+ggplot2::coord_flip()
      }
      
    }else if(input$goplots == "ridgeplot"){
      
      validate(
        need("geneList" %in%  slotNames(gores), "Only the results obtained by FCS can be visualized by this method!")
      )
      goplot <- enrichplot::ridgeplot(gores, showCategory = input$showCategory, core_enrichment = TRUE)
      
    }
    
    return(goplot)
  })
  
  #表格
  output$table_gores <- DT::renderDataTable({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("GO" %in% PISAS_pro@step, "Please run cell GO enrichment analysis reduction first....")
    )
    gores <- PISAS_pro@integrate[["GO"]] 
    
    validate(
      need(nrow(gores@result) > 0, "No gene can be mapped ... or No term enriched under specific pvalueCutoff...")
    )
    
    golist <- gores@result
    
    golist$ID <- paste0("<a href='http://amigo.geneontology.org/amigo/term/",golist$ID,"' target='_blank'>",golist$ID,"</a>")
    
    return(DT::datatable(golist, options = list(orderClasses = TRUE,scrollX = TRUE),escape=FALSE))
    
  })
  
  output$download_golist <- downloadHandler(
    
    filename = "golist.csv",
    
    content = function(file) {
      
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      validate(
        need("GO" %in% PISAS_pro@step, "Please run GO enrichment analysis first....")
      )
      
      gores <- PISAS_pro@integrate[["GO"]] 
      
      gores <- gores@result
      
      write.csv(gores, file)
    }
  )
  
  
  output$gsea_geneSetID <- renderUI({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("GO" %in% PISAS_pro@step, "")
    )
    gores <- PISAS_pro@integrate[["GO"]] 
    
    selectInput("geneSetID", "GeneSet ID:", as.list(gores@result$ID),multiple=TRUE)
  })
  
  output$GSEA_plot <- renderPlot({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("GO" %in% PISAS_pro@step, "Please run cell differential expression analysis reduction first....")
    )
    gores <- PISAS_pro@integrate[["GO"]]
    
    validate(
      need("geneList" %in%  slotNames(gores), "Only the results obtained by FCS can be visualized by this method!")
    )
    indx <- which(gores@result$ID %in% input$geneSetID)
    message(indx)
    n <- length(indx)
    message(n)
    if(n>0){
      color_pairs <- hue_pal()(n)
      plot <- enrichplot::gseaplot2(gores, geneSetID = indx, pvalue_table = TRUE,
                                    color = color_pairs, ES_geom = input$go_style)
      print(plot)
      return(plot)
    }else{
      return(NULL)
    }
    
  })
  
  #KEGG绘图
  output$keggPlot <- renderPlot({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("KEGG" %in% PISAS_pro@step, "Please run KEGG pathway enrichment analysis first....")
    )
    keggres <- PISAS_pro@integrate[["KEGG"]]
    
    validate(
      need(nrow(keggres@result) > 0, "No gene can be mapped ... or No term enriched under specific pvalueCutoff...")
    )
    keggplot <- NULL
    if(input$keggplots == "dotplot"){
      keggplot <- enrichplot::dotplot(keggres,showCategory=input$showCategory2,color=input$gocolorBy)
      
    }else if(input$keggplots == "cnetplot"){
      keggplot <-  cnetplot(keggres,showCategory = input$showCategory2,categorySize=input$gocolorBy,circular = TRUE, colorEdge = TRUE)
      
    }else if(input$keggplots == "emapplot"){
      keggplot <- emapplot(keggres, showCategory = input$showCategory2, color = input$gocolorBy)
      
    }else if(input$keggplots == "heatplot"){
      if("geneList" %in%  slotNames(keggres)){
        keggplot <- heatplot(keggres, showCategory = input$showCategory2, foldChange = keggres@geneList)
      }else{
        keggplot <- heatplot(keggres, showCategory = input$showCategory2)
      }
    }else if(input$keggplots == "ridgeplot"){
      validate(
        need("geneList" %in%  slotNames(keggres), "Only the results obtained by FCS can be visualized by this method!")
      )
      keggplot <- ridgeplot(keggres, showCategory = input$showCategory, fill = input$gocolorBy, core_enrichment = TRUE)
    }
    
    return(keggplot)
  })
  
  
  output$kegg_pathview <- renderUI({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("KEGG" %in% PISAS_pro@step, "")
    )
    keggres <- PISAS_pro@integrate[["KEGG"]] 
    
    selectInput("pathwayid", "Pathway ID:", as.list(keggres@result$ID))
  })
  
  output$pathwayView <- renderImage({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("KEGG" %in% PISAS_pro@step, "")
    )
    keggres <- PISAS_pro@integrate[["KEGG"]]
    validate(
      need("geneList" %in%  slotNames(keggres), "Only the results obtained by FCS can be visualized by this method!")
    )
    library("pathview")
    species <- PISAS_pro@species
    organisms <- search_kegg_organism(species, by='scientific_name')
    organisms <- organisms$kegg_code[1]
    
    pathview(gene.data  = keggres@geneList,
             pathway.id = input$pathwayid,
             species    = organisms,
             limit      = list(gene=max(abs(keggres@geneList)), cpd=1))
    
    dir <- getwd()
    outfile <- paste0(dir,"/",input$pathwayid,".png")
    file.remove(paste0(dir,"/",input$pathwayid,".pathview.png"))
    file.remove(paste0(dir,"/",input$pathwayid,".xml"))
    list(src = outfile,
         contentType = 'image/png',
         width = 800,
         alt = "This is alternate text")
  }, deleteFile = TRUE)
  
  
  #表格
  output$table_keggres <- DT::renderDataTable({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("KEGG" %in% PISAS_pro@step, "Please run cell KEGG pathway enrichment analysis reduction first....")
    )
    keggres <- PISAS_pro@integrate[["KEGG"]]
    
    validate(
      need(nrow(keggres@result) > 0, "No gene can be mapped ... or No term enriched under specific pvalueCutoff...")
    )
    
    kegglist <-  keggres@result
    kegg_urls <- paste0("http://www.kegg.jp/kegg-bin/show_pathway?", kegglist$ID, "/", keggres[kegglist$ID, "geneID"])
    
    kegglist$ID <- paste0("<a href='", kegg_urls,"' target='_blank'>",kegglist$ID,"</a>")
    
    return(DT::datatable(kegglist, options = list(orderClasses = TRUE,scrollX = TRUE),escape=FALSE))
    
  })
  
  output$download_kegglist <- downloadHandler(
    
    filename = "kegg.csv",
    
    content = function(file) {
      
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      validate(
        need("KEGG" %in% PISAS_pro@step, "")
      )
      
      keggres <- PISAS_pro@integrate[["KEGG"]]
      
      kegglist <-  keggres@result
      
      write.csv(kegglist, file)
    }
  )
  
  
  ## 
  
  ## 6 pseudotime ------------------------------------------------------------
  
  output$info_pseudo<-renderText(Info_out$out_pseudo)
  
  pseudo_vals <- shiny::reactiveValues(keeprows = NULL)
  ica_space_df <- data.frame()
  edge_df <- data.frame()
  
  observeEvent(input$btn_pseudotime, {
    tryCatch({
      disable("btn_pseudotime")
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      validate(
        need("DE" %in% PISAS_pro@step, "")
      )
      step <- PISAS_pro@step
      type <- PISAS_pro@type
      scRNA <- PISAS_pro@integrate$Combind
      assay.map <- PISAS_pro@integrate[["assay.map"]]
      
      withProgress(message = 'Trajectory analysis...', value = 0.3, {
        Info_out$out_pseudo <- c("<h3>Trajectory analysis</h3>")
        
        require(monocle3)
        setProgress(value = 0.4,detail="UMAP dimensionality reduction")
        names(scRNA@reductions)[names(scRNA@reductions) == assay.map ] <- "UMAP"
        #scRNA <- RunUMAP(scRNA, dims = 1:5, reduction.name = "UMAP")
        cds <- as.cell_data_set(scRNA)
        cds <- cluster_cells(cds)
        Info_out$out_pseudo <-paste0(Info_out$out_pseudo,c("Learn principal graph...<br>"))
        
        cds <- learn_graph(cds)
        
        # cds <- order_cells(cds)
        reduction_method <- "UMAP"
        reduced_dim_coords <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst)
        
        ica_space_df <- as.data.frame(reduced_dim_coords)
        
        num_reduced_dim <- ncol(ica_space_df)
        if (num_reduced_dim >= 3) {
          use_3d = TRUE
        }else {
          use_3d = FALSE
        }
        
        colnames(ica_space_df) <- vapply(seq_along(ica_space_df), 
                                         function(i) {
                                           paste0("prin_graph_dim_", i)
                                         }, c("a"))
        
        ica_space_df$sample_name <- row.names(ica_space_df)
        ica_space_df$sample_state <- row.names(ica_space_df)
        
        ica_space_df <<- ica_space_df
        
        dp_mst <- principal_graph(cds)[[reduction_method]]
        if(use_3d){
          edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select_(source = "from", 
                                                                           target = "to") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                                 dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                                source_prin_graph_dim_2 = "prin_graph_dim_2", 
                                                                                                                                source_prin_graph_dim_3 = "prin_graph_dim_3"), 
                                                                                                               by = "source") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                                                                     dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                                                                    target_prin_graph_dim_2 = "prin_graph_dim_2", 
                                                                                                                                                                    target_prin_graph_dim_3 = "prin_graph_dim_3"), 
                                                                                                                                                   by = "target")
        }else {
          edge_df <- dp_mst %>% igraph::as_data_frame() %>% dplyr::select_(source = "from", 
                                                                           target = "to") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                                 dplyr::select_(source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                                source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                               by = "source") %>% dplyr::left_join(ica_space_df %>% 
                                                                                                                                                     dplyr::select_(target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                                                                                                                                    target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                                                                                                                                                   by = "target")
        }
        
        num_roots <- nrow(ica_space_df)
        sel <- rep(FALSE, nrow(ica_space_df))
        
        pseudo_vals$keeprows <- rep(TRUE,nrow(ica_space_df))
        
        #
        # cds <- order_cells(cds,root_pr_nodes = "Y_41")
        # vv <- plot_cells(cds, color_cells_by = "partition", label_cell_groups = FALSE, label_leaves = FALSE, 
        #            label_branch_points = FALSE)+
        #   geom_point(data = ica_space_df, mapping = aes(x = prin_graph_dim_1, y = prin_graph_dim_2,label=sample_name))+
        #   geom_text(data = ica_space_df,mapping = aes(prin_graph_dim_1,prin_graph_dim_2,label=sample_name),hjust=0, vjust=0)
        # 
        # ggplotly(vv,tooltip="sample_name")
        
        
        
        
        # pseudo <- toMonocle(scRNA, import_all = FALSE)
        # pseudo <- estimateSizeFactors(pseudo)
        # #pseudo <- estimateDispersions(pseudo)
        # pseudo <- detectGenes(pseudo, min_expr = 0.1)
        # pseudo@featureData@data$use_for_ordering <- pseudo@featureData@data$num_cells_expressed > 0.05 * ncol(pseudo)
        # pseudo <- setOrderingFilter(pseudo, ordering_genes  = difMarker$gene)
        # pseudo <- reduceDimension(pseudo, method = "SimplePPT",norm_method="none")
        # pseudo <- orderCells(pseudo)
        
        PISAS_pro@integrate[["pseudo"]] <- cds
        
        # parts <- cds@clusters$UMAP$partitions
        # 
        # addpart <- pseudo$Pseudotime
        # 
        # names(addpart) <- colnames(pseudo)
        # 
        # addpart[addpart == min(addpart[names(parts[parts==1])])]
        # 
        PISAS_syn$PISAS_pro <<- PISAS_pro
        
        setProgress(1)
      })
      enable("btn_pseudotime")
    },
    error=function(e){
      enable("btn_pseudotime")
      message(safeError(e))
      return(NULL)
    })
    
  })
  
  output$nodeplot <- renderPlot({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("pseudo" %in% names(PISAS_pro@integrate), "Please click [RUN] on the left first!")
    )
    cds <- PISAS_pro@integrate[["pseudo"]]
    if(is.null(pseudo_vals$keeprows)) return()
    keep <- ica_space_df[pseudo_vals$keeprows, , drop = FALSE]
    exclude <- ica_space_df[!pseudo_vals$keeprows, , drop = FALSE]
    reduction_method = "UMAP"
    reduced_dims <- as.data.frame(reducedDims(cds)[[reduction_method]])
    colnames(reduced_dims) <- vapply(seq_along(colnames(reduced_dims)), 
                                     function(i) {
                                       paste0("V", i)
                                     }, c("a"))
    
    ggplot(keep, aes(prin_graph_dim_1, prin_graph_dim_2)) + 
      geom_point(data = reduced_dims, aes(x = V1, 
                                          y = V2), size = 0.5, color = "gray", alpha = 0.3) + 
      geom_point(alpha = 0.7) + geom_point(data = exclude, 
                                           shape = 21, fill = "red", color = "red") 
  }, height = function() {
    session$clientData$output_nodeplot_width
  })
  
  shiny::observeEvent(input$nodeplot_click, {
    tryCatch({
      Info_out$out_pseudo <-paste0(Info_out$out_pseudo,c("Click root nodes...<br>"))
      
      res <- shiny::nearPoints(ica_space_df, xvar = "prin_graph_dim_1", 
                               yvar = "prin_graph_dim_2", input$nodeplot_click, 
                               allRows = TRUE)
      pseudo_vals$keeprows <- xor(pseudo_vals$keeprows, res$selected_)
    },
    error=function(e){
      
      message(safeError(e))
      return(NULL)
    })
  })
  
  shiny::observeEvent(input$choose_toggle, {
    tryCatch({
      if(is.null(input$nodeplot_brush)) return()
      res <- shiny::brushedPoints(ica_space_df, input$nodeplot_brush, 
                                  xvar = "prin_graph_dim_1", yvar = "prin_graph_dim_2", 
                                  allRows = TRUE)
      pseudo_vals$keeprows <- xor(pseudo_vals$keeprows, res$selected_)
    },
    error=function(e){
      
      message(safeError(e))
      return(NULL)
    })
  })
  shiny::observeEvent(input$pseudo_reset, {
    pseudo_vals$keeprows <- rep(TRUE, nrow(ica_space_df))
  })
  
  shiny::observeEvent(input$pseudo_done, {
    tryCatch({
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      step <- PISAS_pro@step
      validate(
        need("pseudo" %in% names(PISAS_pro@integrate), "")
      )
      cds <- PISAS_pro@integrate[["pseudo"]]
      
      if(is.null(pseudo_vals$keeprows)) return()
      
      sel <- pseudo_vals$keeprows
      
      validate(
        need(sum(!sel) > 0, "Please click the 'Choose' button first to confirm the selected root cells")
      )
      select_pr_nodes <- as.character(ica_space_df$sample_name[which(!sel)])
      
      cds <- order_cells(cds,root_pr_nodes = select_pr_nodes)
      
      PISAS_pro@integrate[["pseudo"]] <- cds
      reduction_method = "UMAP"
      print(cds@principal_graph_aux[[reduction_method]]$pseudotime)
      PISAS_pro@step <- unique(c(step,"pseudo"))
      PISAS_syn$PISAS_pro <<- PISAS_pro
      Info_out$out_pseudo <-paste0(Info_out$out_pseudo,c("Finished!<br>"))
    },
    error=function(e){
      message(safeError(e))
      return(NULL)
    })
  })
  
  output$pseu_color <- renderUI({
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    if(is.null(PISAS_pro)) return()
    
    meta <- PISAS_pro@integrate$Combind@meta.data
    group_lists <- c("pseudotime","partition",as.list(colnames(meta)))
    selectInput("pse_col", "Color cells by:",group_lists, selected = "pseudotime")
  })
  
  output$pseudotimePlot1 <- renderPlotly({
    
    validate(
      need("monocle3" %in% installed.packages()[, "Package"], "The monocle3 package has not been installed. If you want to use this module, please follow https://cole-trapnell-lab.github.io/monocle3/ to install monocle3")
    )
    require(monocle3)
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    
    validate(
      need("pseudo" %in% PISAS_pro@step, "Note:Please select trajectory roots first....")
    )
    cds <-  PISAS_pro@integrate[["pseudo"]]
    
    validate(
      need("pseudotime" %in% names(cds@principal_graph_aux$UMAP), "Note: No psedotime for UMAP calculated. Please select the trajectory roots first to run the order_cells!")
    )
    
    ply <- monocle3::plot_cells(cds,
                                color_cells_by = input$pse_col,
                                label_cell_groups=TRUE,
                                label_leaves=TRUE,
                                label_branch_points=TRUE,
                                graph_label_size=2.5)
    
    ply <- ggplotly(ply,autosize = T) %>% config(displaylogo = FALSE)
    
    ggiris <- qplot(Petal.Width, Sepal.Length, data = iris, color = Species)
    ggplotly(ggiris)
    
    return(ply)
  })
  
  
  output$table_pheno <- DT::renderDataTable({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("pseudo" %in% names(PISAS_pro@integrate), "Please run this step first....")
    )
    cds <-  PISAS_pro@integrate[["pseudo"]]
    
    pdata <- data.frame(partitions = cds@clusters$UMAP$clusters, pseudotime =cds@principal_graph_aux$UMAP$pseudotime)
    
    return(DT::datatable(pdata, options = list(orderClasses = TRUE,scrollX = TRUE)))
    
  })
  
  output$download_pheno <- downloadHandler(
    
    filename = "pheno.csv",
    
    content = function(file) {
      
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      if (!("pseudo" %in% names(PISAS_pro@integrate)))
        return()
      
      pseudo <-  PISAS_pro@integrate[["pseudo"]]
      
      write.csv(pData(pseudo), file)
    }
  )
  
  
  ##  
  ### 7 cluster annotation ----
  
  ## listen to methods_anno
  
  output$info_anno<-renderText(Info_out$out_anno)
  
  markers.set <- reactive({
    
    markers.set <- NULL
    
    tryCatch(
      {
        if(is.null(input$file_markers)) return(NULL)
        filepath <- input$file_markers$datapath
        markers.set<- read.csv(filepath, header=TRUE, sep="\t",row.names = 1,check.names=FALSE,stringsAsFactors=FALSE)
      },
      error = function(e) {
        return(safeError(e))
      }
    )
    return(markers.set)
  })
  
  anno.set <- reactive({
    
    anno.set <- NULL
    
    tryCatch(
      {
        if(is.null(input$file_annotation)) return(NULL)
        filepath <- input$file_annotation$datapath
        anno.set<- read.csv(filepath, header=TRUE, sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
        
      },
      error = function(e) {
        return(safeError(e))
      }
    )
    return(anno.set)
  })
  
  
  observeEvent(input$btn_annotation, {
    tryCatch({
      disable("btn_annotation")
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      validate(
        need("DE" %in% PISAS_pro@step, "")
      )
      step <- PISAS_pro@step
      
      difMarker <- PISAS_pro@integrate$Markers
      de.info <- PISAS_pro@integrate$DE.info
      
      type <- PISAS_pro@type
      scRNA <- PISAS_pro@integrate$Combind
      Idents(scRNA) <- de.info
      withProgress(message = 'Cluster annotation ...', value = 0.3, {
        
        Info_out$out_anno <- c("<h3>Cluster annotation</h3>")
        
        library(scCATCH)
        
        
        if(input$methods_anno == "scCATCH"){
          message(input$species_anno)
          if(input$species_anno != "Others"){
            
            res <- AutoSingleType(difMarker,species=input$species_anno,tissue= input$tissue.anno,auto=input$autoMarkers)
            PISAS_pro@integrate[["annotation"]] <- res[[1]]
            new.labels <- res$annotation$cell_type
          }else{
            Info_out$out_anno <- paste0(Info_out$out_anno,c("Please select the AUCell method to annotate the clusters. Failed!<br>"))
            return(NULL)
          }
        }else if(input$methods_anno == "AUCell"){
          Info_out$out_anno <- paste0(Info_out$out_anno,c("Evidence-based score and annotation for each cluster...<br>"))
          
          library(AUCell)
          req(markers.set())
          markers.set <- markers.set()
          xy.list <- list()
          for (i in 1:nrow(markers.set)) {
            xy.list[[rownames(markers.set)[i]]] <- gsub(" ","",unlist(strsplit(markers.set$cellMarker[i],",")))
          }
          
          afm <- scRNA@assays$RNA@data
          pooled <- matrix(nrow=nrow(afm), ncol = 0)
          for (i in levels(scRNA@active.ident)) {
            
            m <- afm[,which(scRNA@active.ident==i)]
            pooled <- cbind(pooled, rowSums(m)/ncol(m))
          }
          colnames(pooled) <- levels(scRNA@active.ident)
          
          rankings <- AUCell_buildRankings(pooled,plotStats=FALSE, verbose=FALSE)
          cell.aucs <- AUCell_calcAUC(xy.list, rankings,aucMaxRank=nrow(rankings)*0.05,verbose=FALSE)
          results <- t(cell.aucs@assays$data$AUC)
          PISAS_pro@integrate[["annotation"]] <- results
          new.labels <- colnames(results)[max.col(results)]
          
        }else{
          Info_out$out_anno <- paste0(Info_out$out_anno,c("Add user's own annotation for each cell...<br>"))
          req(anno.set())
          anno.set <- anno.set()
          anno.set <- arrange(anno.set,cluster)
          PISAS_pro@integrate[["annotation"]] <- anno.set
          new.labels <- anno.set$cellType
          
        }
        scRNA@meta.data$annotations <- FoldNames(new.labels,scRNA)
        PISAS_pro@integrate$Combind <- scRNA
        PISAS_pro@integrate[["anno_func"]] <- input$methods_anno
        Info_out$out_anno <- paste0(Info_out$out_anno,c("Adding annotation meta-data...<br>"))
        
        PISAS_pro@step <- unique(c(step,"annotation"))
        PISAS_syn$PISAS_pro <<- PISAS_pro
        
        setProgress(1)
      })
      enable("btn_annotation")
      Info_out$out_anno <- paste0(Info_out$out_anno,c("Finished!<br>"))
      
    },
    error=function(e){
      enable("btn_annotation")
      Info_out$out_anno <- paste0(Info_out$out_anno,safeError(e),"<br>Failed!<br>")
      return(NULL)
    })
  })
  
  output$annotation_plot <- renderPlotly({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("annotation" %in% PISAS_pro@step, "")
    )
    
    type <- PISAS_pro@type
    scRNA <- PISAS_pro@integrate$Combind
    
    if(type == "Single"){
      if("umap" %in% names(scRNA@reductions)){
        input.dat <- "umap"
      }else if("tsne" %in% names(scRNA@reductions)){
        input.dat <- "tsne"
      }else{
        input.dat <- "pca"
      }
    }else{
      
      input.dat <- PISAS_pro@integrate[["assay.map"]]
    }
    Idents(scRNA) <- "annotations"
    DimPlotly(scRNA,reduction = input.dat,label=TRUE,pt.size = 3,label.size = 16)
  })
  
  output$table_anno <- DT::renderDataTable({
    
    if(is.null(PISAS_syn$PISAS_pro)){
      return(NULL)
    }
    PISAS_pro <- PISAS_syn$PISAS_pro
    validate(
      need("annotation" %in% names(PISAS_pro@integrate), "Please run this step first....")
    )
    annotation <-  PISAS_pro@integrate[["annotation"]]
    func <- PISAS_pro@integrate[["anno_func"]]
    if(func=="scCATCH") annotation <- annotation[,-2]
    return(DT::datatable(annotation, options = list(orderClasses = TRUE,scrollX = TRUE)))
    
  })
  
  output$download_anno <- downloadHandler(
    
    filename = "annotation.csv",
    
    content = function(file) {
      if(is.null(PISAS_syn$PISAS_pro)){
        return(NULL)
      }
      PISAS_pro <- PISAS_syn$PISAS_pro
      if (!("annotation" %in% names(PISAS_pro@integrate)))
        return()
      annotation <-  PISAS_pro@integrate[["annotation"]]
      
      write.csv(annotation, file)
    }
  )
  
  
  
  
  ### Save: save the project------------------------------
  
  observeEvent(input$savetest, {
    if(is.null(PISAS_syn$PISAS_pro)){
      shinyalert("Warning!", "Please load the project first!", type = "error")
      return(NULL)
    }
    username <- as.character(pipelines$username)
    message("save by: ",username)
    withProgress(message = paste0('Saving the project by ',username), value = 0.1, {
      for (i in 1:7) {
        incProgress(0.1, detail = "This will take a long time.Storage time depends on your data volume and operating steps.")
        Sys.sleep(0.1)
      }
      #-----update the project list
      
      PISAS_pro <- PISAS_syn$PISAS_pro
      len <- length(PISAS_pro@array)
      
      list <- paste0(taskdir,"/Projectlist.Rds")
      create_time <- Sys.time()
      
      create.time <- as.character(format(create_time,format='%Y/%m/%d %H:%M:%S'))
      
      Size <- as.character(format(object.size(PISAS_pro), units = "auto"))
      
      message("Project_ID:",Project_ID)
      
      if(is.null(Project_ID)){
        #if project id does not exist, create info
        ids <- stri_rand_strings(1,9)
        
        if (!(file.exists(list))){
          #If list does not exists, create
          dslist <- data.frame(ID= ids,Project=input$projectName,username=username,Type=input$Type,Size= Size,
                               Number = len ,Seed=input$randomSeed,create.time= create.time,update.time= create.time,stringsAsFactors = FALSE)
          saveRDS(dslist,file=list)
        }else{
          #If list exists, update
          dslist <- readRDS(list)
          while(ids %in% dslist$ID){
            #avoid duplication
            ids <- stri_rand_strings(1,9)
          }
          
          newlist <- data.frame(ID= ids,Project=input$projectName,username=username,Type=input$Type,Size= Size,
                                Number = len,Seed=input$randomSeed,create.time= create.time,update.time= create.time)
          
          dslist <- rbind(dslist,newlist)
          saveRDS(dslist,file=list)
        }
      }else{
        #update project info
        message("update project info")
        ids <- Project_ID
        dslist <- readRDS(list)
        message("nns",nrow(dslist))
        if(ids %in% dslist$ID){
          is <- dslist$username == username & dslist$ID == ids
          
          # dslist[is,] <- data.frame(ID= ids,Project=input$projectName,username=username,Type=input$Type,Size= Size,
          #                           Number = len,Seed=input$randomSeed,create.time= dslist[is,]$create.time,update.time= create.time)
          #
          
          dslist$Size <- as.character(dslist$Size)
          dslist[is,]$Size <- Size
          dslist[is,]$Number <- as.character(len)
          dslist[is,]$Seed <- as.character(input$randomSeed)
          dslist$update.time <- as.character(dslist$update.time)
          dslist[is,]$update.time <- as.character(create.time)
        }
        saveRDS(dslist,file=list)
      }
      #-----save the Rdata
      filepath <- paste0(taskdir,username,"/project")
      if (!(file.exists(filepath))){
        dir.create(file.path(filepath),recursive = TRUE)
      }
      
      RDataName <- paste0(ids,".RData")
      save(PISAS_pro,file=paste0(filepath,"/",RDataName))
      shinyalert("OK!", "The project has been successfully saved", type = "success")
      VaryPro$total <- VaryPro$total+1
      Project_ID <<- ids
    })
  })
  
}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
