library(shiny)
library(psych)
options(shiny.maxRequestSize = 30*1024^2)
# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Parallel Analysis"),
  
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
      
      #Input: Select results to display ----
      # selectInput("matrix",
      #             label = "Display the eigenvalues for",
      #              choices = list("Factor Analysis"='fa',
      #                          "Principal Components"='pc',
      #                          "Both"='both'),
      #              selected = "Factor Analysis"),

      # Horizontal line ----
      #tags$hr(),
      
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
      uiOutput("tab"),
      uiOutput("tab2"),
      
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
)

# Define server logic to read selected file ----
server <- function(input, output) {
  url <- a("M2PL estimation", href="https://www.google.com/")
  output$tab <- renderUI({
    tagList("Other Functions:", url)
  })
  url2 <- a("M3PL estimation", href="https://www.google.com/",)
  output$tab2 <- renderUI({
    tagList(url2)
  })
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
}
# Run the app ----
shinyApp(ui, server)