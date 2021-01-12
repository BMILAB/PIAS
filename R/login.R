
loginUI <- function(id, title = "", user_title = "User Name", pass_title = "Password",
                    login_title = "Login", error_message = "Invalid username or password!",inactivated_message="Your account has not been activated!") {
  ns <- shiny::NS(id)

  shiny::div(id = ns("panel"), style = "width: 500px; max-width: 100%;margin-bottom: 110px; padding: 20px;border-top: 3px solid #0071b9;box-shadow: 0 2px 8px rgba(0, 0, 0, 0.2);",
      shiny::wellPanel(
        style = "    background-color: #fff;border: none",
        shiny::tags$h2(title, class = "text-center", style = "padding-top: 0;"),

        shiny::textInput(ns("user_name"), shiny::tagList(shiny::icon("user"), user_title)),

        shiny::passwordInput(ns("password"), shiny::tagList(shiny::icon("unlock-alt"), pass_title)),

        shiny::div(
          style = "text-align: center;",
          shiny::actionButton(ns("button"), login_title, class = "btn-primary btn-block", style = "color: white;"),
          shiny::tags$a("Register",class="btn btn-warning btn-block",href="http://219.229.80.154/pscap/index.php",target="_blank")
          
        ),
        shiny::br(),
        shiny::div(
          style = "text-align: right;font-size:12px",
          shiny::tags$a("Forgot username / password?",href="http://219.229.80.154/pscap/reset.php",target="_blank")
        ),
		

        shinyjs::hidden(
          shiny::div(id = ns("error"),
                     shiny::tags$p(error_message,
                     style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
        ),
        shinyjs::hidden(
          shiny::div(id = ns("inactivated"),
                     shiny::tags$p(inactivated_message,
                                   style = "color: red; font-weight: bold; padding-top: 5px;", class = "text-center"))
        )
      )
  )
}

#' md5 + salt
#' 
md5hash <- function(pswd){
  psw <- md5(pswd)
  salt <- substr(psw,2,22)
  pswsalt <- paste0(salt,psw)
  passwd <- md5(pswsalt)
  return(passwd)
}

login <- function(input, output, session, data, user_col, pwd_col,
                  hashed = FALSE, algo = c("md5", "sha1", "crc32", "sha256", "sha512", "xxhash32", "xxhash64", "murmur32"), 
                  log_out = NULL) {
  
  algo <- match.arg(algo, several.ok = FALSE)

  credentials <- shiny::reactiveValues(user_auth = FALSE, info = NULL)

  shiny::observeEvent(log_out(), {
    credentials$user_auth <- FALSE
    credentials$info <- NULL
    shiny::updateTextInput(session, "password", value = "")
  })

  shiny::observeEvent(credentials$user_auth, ignoreInit = TRUE, {
    shinyjs::toggle(id = "panel")
  })

  users <- dplyr::enquo(user_col)
  pwds <- dplyr::enquo(pwd_col)
  
  # ensure all text columns are character class
  
  #data <- dplyr::mutate_if(user_base, is.factor, as.character)
  
  shiny::observeEvent(input$button, {
    
    # check for match of input username to username column in data
    # 
    #row_username <- which(dplyr::pull(data, !! users) == input$user_name)
    #
    con <- dbConnect(RMariaDB::MariaDB(), username = "root", password = "bmilab",dbname = "pscap")
    kk <- input$user_name
    
    sql1 <- paste0("SELECT * FROM members WHERE username ='", kk,"';")
    
    res <- dbSendQuery(con,sql1)
    
    result <- dbFetch(res)
    if(nrow(result) == 1){
      #password is true
      userpass <- md5hash(input$password)
      
      if(result$password == userpass){
        if(result$active == "Yes"){
          credentials$user_auth <- TRUE
          credentials$info <- result$username
        }else{
          shinyjs::toggle(id = "inactivated", anim = TRUE, time = 1, animType = "fade")
          shinyjs::delay(5000, shinyjs::toggle(id = "inactivated", anim = TRUE, time = 1, animType = "fade"))
        }
      }else{
        shinyjs::toggle(id = "error", anim = TRUE, time = 1, animType = "fade")
        shinyjs::delay(5000, shinyjs::toggle(id = "error", anim = TRUE, time = 1, animType = "fade"))
      }
      
    }else{
      shinyjs::toggle(id = "error", anim = TRUE, time = 1, animType = "fade")
      shinyjs::delay(5000, shinyjs::toggle(id = "error", anim = TRUE, time = 1, animType = "fade"))
    }
    
    dbClearResult(res)
    
    # Disconnect from the database
    dbDisconnect(con)
  })

  # return reactive list containing auth boolean and user information
  shiny::reactive({
    shiny::reactiveValuesToList(credentials)
  })

}
