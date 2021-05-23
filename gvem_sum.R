library(shiny)
library(psych)
library(Rcpp)
library(testit)
library(Matrix)
library(gtools)
source("initialization.r")
source("cfa_2pl.r")
source("efa_2pl.r")
source("cfa_3pl.r")
source("efa_3pl.r")
source("lassoc1_2pl.r")
source("lassoc2_2pl.r")
source("lassoc1_3pl.r")
source("lassoc2_3pl.r")
source("adptc1_2pl.r")
source("adptc2_2pl.r")
source("adptc1_3pl.r")
source("adptc2_3pl.r")
options(shiny.maxRequestSize = 30*1024^2)
# Define UI for data upload app ----
ui <- navbarPage("VEMIRT Shiny App",
                 
                 # App title ----
                 tabPanel("Parallel Analysis",
                          
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              # Input: Select a file ----
                              fileInput("file", "Choose Data CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              # Input: Select quotes ----
                              radioButtons("quote", "Quote",
                                           choices = c(None = "",
                                                       "Double Quote" = '"',
                                                       "Single Quote" = "'"),
                                           selected = '"'),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Select number of rows to display ----
                              radioButtons("disp", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              
                              #Input: Select correlation methods ----
                              selectInput("cor", 
                                          label = "Choose a correlation method",
                                          choices = list("Pearson"='cor', 
                                                         "Polychoric"='poly'),
                                          selected = "poly"),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              
                              # Input: Select the number of simulated analyses to perform ----
                              sliderInput("niter",
                                          label = "Number of simulated analyses to perform",
                                          min = 0, max = 100, value = 20),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Display the eigenvalue plot  ----
                              radioButtons("plot",
                                           label = "Display the eigenvalue plot",
                                           choices = c("TRUE" = "T",
                                                       "FALSE" = "F"),
                                           selected = "T"),
                              # Horizontal line ----
                              tags$hr(),
                              actionButton("go", "Run")
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              h2("Data"),
                              # Output: Data file ----
                              tableOutput("contents"),
                              
                              h1(textOutput("nfa")),
                              
                              #strong(textOutput("nfa")),
                              
                              #strong("Parallel analysis suggests that the number of factors =  "),
                              
                              h2("Eigenvalues"),
                              #Output: EigenValues,
                              tableOutput("eigen"),
                              
                              h2("Eigenvalue Plot"),
                              #Output:Scree Plot,
                              plotOutput("graph")
                              
                            )
                            
                          )
                 ),
                 tabPanel("EFA",
                          
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              # Input: Select a file ----
                              fileInput("file1", "Choose Data CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header1", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep1", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Select number of rows to display ----
                              radioButtons("disp1", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              # Horizontal line ----
                              tags$hr(),
                              #Input: Select EFA methods ----
                              selectInput("method", 
                                          label = "Choose an EFA method",
                                          choices = list("Rotation"='Rot', 
                                                         "Adaptive Lasso"='AL',
                                                         "Lasso"="La"),
                                          selected = "Rot"),
                              
                              # Horizontal line ----
                              tags$hr(),
                              #Choose Constraint 1
                              conditionalPanel(
                                condition = "input.method == 'AL' || input.method == 'La'",
                                radioButtons("cons", "Constraint",
                                             choices = c("Constraint 1" = "c1",
                                                         "Constraint 2" = "c2"),
                                             selected = "c1"),
                                helpText("Note: For constraint 1, you need to set a K by K submatrix of the indicator 
                                         matrix to be an identity matrix. For constraint 2, you need to set a K by K submatrix of the indicator 
                                         matrix to be a triangular matrix. ") 
                              ),
                              
                              # Horizontal line ----
                              tags$hr(),
                              conditionalPanel(
                                condition = "(input.method == 'AL' || input.method == 'La') && input.cons == 'c1' ",
                                # Input: Select the indicator dile ----
                                fileInput("file5", "Choose Loading Indicator CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                helpText("Note:  The loading indicator file should be a binary Item by Factor matirx. 
                                         0 refers to the item is irrelevant with this factor, 1 otherwise.")
                              ),
                              #input the domain
                              tags$hr(),
                              numericInput("dom", "domain", min=1, max=10, value=5),
                              
                              # Input: Select the model to perform ----
                              tags$hr(),
                              selectInput("mod1", "Choose the model to fit the data",
                                          c("M2PL" = "2pl",
                                            "M3PL" = "3pl")),
                              
                              
                              
                              conditionalPanel(
                                condition = "input.mod1 == '3pl'",
                                sliderInput("fr1", "Step size", min = 0.5, max = 1.0, value = 0.51),
                                numericInput("sample1", "Subsample", min=1, max=1000, value=10),
                                numericInput("bmu2", "Prior Distribution for b (Mean)", min=-3, max=3, value=0),
                                numericInput("bsig2", "Prior Distribution for b (Variance)", min=0.1, max=20, value=1),
                                numericInput("alpha2", "Prior Distribution for g (Alpha)", min=0, max=50, value=2),
                                numericInput("beta2", "Prior Distribution for g (Beta)", min=0, max=50, value=5)
                              ),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              actionButton("go1", "Run"),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              radioButtons("checkGroup1", "Download Results",
                                           choices = list("All Results" = "all",
                                                          "Item Parameters" = "item",
                                                          "Covariance Matrix" = "cov"),
                                           selected = "all"),
                              # Button
                              downloadButton("downloadData", "Download Results")
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              h2("Data"),
                              # Output: Data file ----
                              tableOutput("contents1"),
                              
                              #warning
                              h1(textOutput("warn")),
                              
                              h2("Item Parameter Results"),
                              tableOutput("par1"),
                              
                              h2("Covariance Matrix"),
                              tableOutput("cov1"),
                              
                              h2("Q Matrix"),
                              tableOutput("fl1"),
                              
                              h2("Model Fit"),
                              textOutput("gic1"),
                              
                            )
                            
                          )
                 ),
                 tabPanel("CFA",
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              
                              # Input: Select a file ----
                              fileInput("file3", "Choose Data CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header3", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep3", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # # Input: Select number of rows to display ----
                              radioButtons("disp3", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              #Horizontal line ----
                              tags$hr(),
                              
                              # Input: Select the indicator dile ----
                              fileInput("file4", "Choose Loading Indicator CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              helpText("Note: The loading indicator file should be a binary Item by Factor matirx. 
               0 refers to the item is irrelevant with this factor, 1 otherwise."),
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header4", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep4", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              #Input: Select number of rows to display ----
                              radioButtons("disp4", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              # Input: Select the model to perform ----
                              tags$hr(),
                              selectInput("mod2", "Choose the model to fit the data",
                                          c("M2PL" = "2pl",
                                            "M3PL" = "3pl")),
                              
                              
                              # radioButtons("mod",
                              #              label = "Choose the model to fit the data",
                              #              choices = c("M2PL" = "2pl",
                              #                          "M3PL" = "3pl"),
                              #              selected = "2pl"),
                              
                              conditionalPanel(
                                condition = "input.mod2 == '3pl'",
                                sliderInput("fr2", "Step size", min = 0.5, max = 1.0, value = 0.6),
                                
                                numericInput("sample2", "Subsample", min=1, max=1000, value=10),
                                numericInput("bmu1", "Prior Distribution for b (Mean)", min=-3, max=3, value=0),
                                numericInput("bsig1", "Prior Distribution for b (Variance)", min=0.1, max=20, value=1),
                                numericInput("alpha1", "Prior Distribution for g (Alpha)", min=0, max=50, value=2),
                                numericInput("beta1", "Prior Distribution for g (Beta)", min=0, max=50, value=5)
                              ),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              actionButton("go3", "Run"),
                              
                              tags$hr(),
                              radioButtons("checkGroup2", "Download Results",
                                                 choices = list("All Results" = "all",
                                                                "Item Parameters" = "item",
                                                                "Covariance Matrix" = "cov"),
                                                 selected = "all"),
                              #Button
                              downloadButton("downloadDataall2", "Download Results")
                              #downloadButton("downloadDataall2", "Download All Results")
                            # downloadButton("downloadDataitem2", "Download Item Parameter Results"),
                            # downloadButton("downloadDatasigma2", "Download Covariance Results")
                            
                              
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              h2("Data"),
                              # Output: Data file ----
                              tableOutput("contents3"),
                              
                              
                              h2("Loading Indicator Matrix"),
                              #Output: indicator matrix,
                              tableOutput("contents4"),
                              
                              
                              h2("Item Parameter Results"),
                              tableOutput("par2"),
                              
                              h2("Covariance Matrix"),
                              tableOutput("cov2"),
                              
                              h2("Model Fit"),
                              textOutput("gic2")
                              
                              
                            )
                            
                          )
                 ))

# Define server logic to read selected file ----
server <- function(input, output,session) {
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, head of that data file by default,
  # or all rows if selected, will be shown.
  
  #parallel analysis results
  data<-reactive({req(input$file)
    df <- read.csv(input$file$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    return(df)
  })
  
  output$contents <- renderTable({
    input$go
    isolate({
      if(input$disp == "head") {
        return(data()[1:6,])
      }
      else {
        return(data())
      }
    })
  })
  
  
  result1<-reactive({fa.parallel(data(),cor=input$cor,fa='fa',n.iter=input$niter,plot=F)
  })
  
  output$nfa<-renderText({
    input$go
    isolate({
      paste("Parallel analysis suggests that the number of factors =", result1()$nfact)
    })
  })
  
  output$eigen<-renderTable({
    input$go
    isolate({
      result1()
      m<-cbind(result1()$fa.values,result1()$fa.sim)
      colnames(m)<-c("Actual Data",
                     "Simulated Data")
      return(m)
    })
  })
  
  output$graph<-renderPlot({
    input$go
    isolate({
      if(input$plot == "T") {
        d<-result1()
        # get the range for the x and y axis
        xrange<-range(1,length(d$fa.values))
        yrange <- range(c(d$fa.values,d$fa.sim))
        
        # set up the plot
        plot(xrange, yrange, type="n", xlab="Factor Number",
             ylab="Eigenvalues" )
        colors <- rainbow(2)
        linetype <- c(1:2)
        plotchar <- seq(18,18+2,1)
        
        # add lines
        lines(1:length(d$fa.values), d$fa.values, type="b", lwd=1.5,
              lty=linetype[1], col=colors[1], pch=plotchar[1])
        lines(1:length(d$fa.values), d$fa.sim, type="b", lwd=1.5,
              lty=linetype[2], col=colors[2], pch=plotchar[2])
        abline(h=1)
        
        # add a title and subtitle
        title("Parallel Analysis Scree Plots")
        
        # add a legend
        legend("topright",legend=c("Actual Data", "Simulated Data"),cex=0.8, col=colors,
               pch=plotchar, lty=linetype)
      }
      else {
        return(NULL)
      }
    })
  })
  
  #EFA
  u<-reactive({req(input$file1)
    df <- data.matrix(read.csv(input$file1$datapath,
                               header = input$header1,
                               sep = input$sep1))
    return(df)
  })
  output$contents1 <- renderTable({
    if(input$disp1 == "head") {
      return(u()[1:6,])
    }
    else {
      return(u())
    }
    
  })
  domain<-reactive({req(input$dom)})
  
  result0<-reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Iteration times", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 10
      }
      progress$set(value = value, detail = detail)
    }
    if (input$method == "Rot" && input$mod1 == '2pl') {
    result<-vem_2PLEFA(u(),domain(),updateProgress)}
    else if (input$method == "Rot" && input$mod1 == '3pl') {
      result<-sgvem_3PLEFA(u(),domain(),input$sample1,input$fr1,input$bmu2,
                           input$bsig2,input$alpha2,
                           input$beta2,updateProgress)}
    else if(input$method == "AL" && input$cons=="c1" && input$mod1 == '2pl'){
      indic<-reactive({req(input$file5)
          df2 <- data.matrix(read.csv(input$file5$datapath))
          return(df2)
        })
      result<-vem_2PLEFA_adaptive_const1_all(u(),domain(),indic(),updateProgress)}
    else if(input$method == "AL" && input$cons=="c1" && input$mod1 == '3pl'){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      result<-sgvem_3PLEFA_adaptive_const1_all(u(),domain(),indic(),
                                             input$sample1,input$fr1,input$bmu2,
                                             input$bsig2,input$alpha2,
                                             input$beta2,updateProgress)}
    else if(input$method == "AL" && input$cons=="c2" && input$mod1 == '2pl'){
      result<-vem_2PLEFA_adaptive_const2_all(u(),domain(),updateProgress)}
    else if(input$method == "AL" && input$cons=="c2" && input$mod1 == '3pl'){
      result<-sgvem_3PLEFA_adaptive_const2_all(u(),domain(),
                                             input$sample1,input$fr1,input$bmu2,
                                             input$bsig2,input$alpha2,
                                             input$beta2,updateProgress)}
    else if(input$method == "La" && input$cons=="c1" && input$mod1 == '2pl'){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      result<-vem_2PLEFA_L1_const1_all(u(),domain(),indic(),updateProgress)}
    else if(input$method == "La" && input$cons=="c1" && input$mod1 == '3pl'){
      indic<-reactive({req(input$file5)
        df2 <- data.matrix(read.csv(input$file5$datapath))
        return(df2)
      })
      result<-sgvem_3PLEFA_L1_const1_all(u(),domain(),indic(),
                                         input$sample1,input$fr1,input$bmu2,
                                         input$bsig2,input$alpha2,
                                         input$beta2,updateProgress)}
    else if(input$method == "La" && input$cons=="c2" && input$mod1 == '2pl'){
      result<-vem_2PLEFA_L1_const2_all(u(),domain(),updateProgress)}
    else if(input$method == "La" && input$cons=="c2" && input$mod1 == '3pl'){
      result<-sgvem_3PLEFA_L1_const2_all(u(),domain(),
                                         input$sample1,input$fr1,input$bmu2,
                                         input$bsig2,input$alpha2,
                                         input$beta2,updateProgress)}
    return(result)
  }
  )
  
  output$warn<-renderText({
    input$go1
    isolate({
      if(input$method != "Rot" && (result0 ()$lbd == 0.1 || result0 ()$lbd == 40)){
        return("Warning: The optimal penalty parameter may be out of range")
      }else{
        return(NULL)
    }})
  })
  
  output$par1<-renderTable({
    input$go1
    isolate({
      if(input$mod1=='2pl'){
      m<-cbind(result0()$ra,result0()$rb)
      colnames(m)<-c(paste("a",1:domain(),sep=""),"b")}else{
        m<-cbind(result0()$ra,result0()$rb,result0()$rc)
        colnames(m)<-c(paste("a",1:domain(),sep=""),"b","g") 
      }
      return(m)
    })
    
  })
  
  
  output$cov1<-renderTable({
    input$go1
    isolate({
      m<-result0()$rsigma
      rownames(m)<-colnames(m)<-paste("dimension",1:domain(),sep="")
      return(m)
    })
    
  })
  output$fl1<-renderTable({
    input$go1
    isolate({
      m<-result0()$Q_mat
      colnames(m)<-paste("dimension",1:input$dom,sep="")
      #rownames(m)<-paste("Item",1:dim(m)[1],sep="")
      return(m)
    })
    
  })
  
  output$gic1<-renderText({
    input$go1
    isolate({
      paste("GIC =", round(result0()$GIC,2))
    })
  })
  
  #Downloadable csv of selected dataset ----
  # output$downloadData <- downloadHandler(
  #   # filename = function() {
  #   #   paste(input$mod, "GVEMresults.rds", sep = "")
  #   # },
  #   filename = "GVEMEFAresults.rds",
  #   content = function(file) {
  #     saveRDS(result0(), file)
  #   }
  # )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if(input$checkGroup1=="all"){
        paste("GVEMEFA",input$checkGroup1 ,"results.rds",sep="")}else{
          paste("GVEMEFA",input$checkGroup1 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup1=="all"){
        saveRDS(result0(), file)}else if(input$checkGroup1=="item"){
          if(input$mod1=='2pl'){
            m<-cbind(result0()$ra,result0()$rb)
            colnames(m)<-c(paste("a",1:domain(),sep=""),"b")}else{
              m<-cbind(result0()$ra,result0()$rb,result0()$rc)
              colnames(m)<-c(paste("a",1:domain(),sep=""),"b","g") 
            }
          write.csv(m,file,row.names = F)
        }else if(input$checkGroup1=="cov"){
          m<-result0()$rsigma
          rownames(m)<-colnames(m)<-paste("dimension",1:domain(),sep="")
          write.csv(m,file,row.names = F)
        }
    }
  )
  
  #CFA output
  u1<-reactive({req(input$file3)
    df <- data.matrix(read.csv(input$file3$datapath,
                               header = input$header3,
                               sep = input$sep3))
    return(df)
  })
  output$contents3 <- renderTable({
    if(input$disp3 == "head") {
      return(u1()[1:6,])
    }
    else {
      return(u1())
    }
    
  })
  
  indic1<-reactive({req(input$file4)
    df2 <- data.matrix(read.csv(input$file4$datapath,
                                header = input$header4,
                                sep = input$sep4))
    return(df2)
  })
  
  output$contents4 <- renderTable({
    
    if(input$disp4 == "head") {
      return(indic1()[1:6,])
    }
    else {
      return(indic1())
    }
    
  })
  domain1<-reactive(dim(indic1())[2])
  
  result<-reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Iteration times", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 10
      }
      progress$set(value = value, detail = detail)
    }
    if(input$mod2=='2pl'){
      result<-vem_2PLCFA(u1(),domain1(),indic1(),updateProgress) 
    }else{
      result<-saem_3PLCFA(u1(),domain1(),indic1(),input$sample2,input$fr2,input$bmu1,
                          input$bsig1,input$alpha1,
                          input$beta1,updateProgress) 
    }
    return(result)
  })
  
  output$par2<-renderTable({
    input$go3
    isolate({
      if(input$mod2=='2pl'){
      m<-cbind(result()$ra,result()$rb)
      colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")}else{
        m<-cbind(result()$ra,result()$rb,result()$rc)
        colnames(m)<-c(paste("a",1:domain1(),sep=""),"b","g") 
      }
      return(m)
    })
    
  })
  
  output$cov2<-renderTable({
    input$go3
    isolate({
      m<-result()$rsigma
      rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
      return(m)
    })
    
  })
  
  output$gic2<-renderText({
    input$go3
    isolate({
      paste("GIC =", round(result()$GIC,2))
    })
  })
  
  #Downloadable csv of selected dataset ----
  output$downloadDataall2 <- downloadHandler(
     filename = function() {
       if(input$checkGroup2=="all"){
       paste("GVEMCFA",input$checkGroup2 ,"results.rds",sep="")}else{
        paste("GVEMCFA",input$checkGroup2 ,"results.csv",sep="")}
          
     },
    content = function(file) {
      if(input$checkGroup2=="all"){
      saveRDS(result(), file)}else if(input$checkGroup2=="item"){
        if(input$mod2=='2pl'){
          m<-cbind(result()$ra,result()$rb)
          colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")}else{
            m<-cbind(result()$ra,result()$rb,result()$rc)
            colnames(m)<-c(paste("a",1:domain1(),sep=""),"b","g") 
          }
        write.csv(m,file,row.names = F)
      }else if(input$checkGroup2=="cov"){
        m<-result()$rsigma
        rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
        write.csv(m,file,row.names = F)
      }
    }
  )
}
# Run the app ----
shinyApp(ui, server)