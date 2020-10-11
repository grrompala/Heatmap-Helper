# Load all packages here
library(shiny)
library(pheatmap)
library(dplyr)
library(tibble)
library(viridis)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(shinyWidgets)

COLORS <- rownames(brewer.pal.info)

options(shiny.maxRequestSize = 30 * 1024 ^ 2)  # max-csv upload set to 30 MB

ui <- fluidPage(
  h1(style="font-family:impact;font-size:300%","The Heaping Heatmap Helper"),
                
              
                sidebarLayout(
                  sidebarPanel(
                   
                    h3("Provide Data Here"),
                     wellPanel(style = "background:lightblue;border-width:thick;border-color:black",
                    # Counts file
                    fileInput(
                      inputId = "file1",
                      label = "Upload counts .csv",
                      accept = c(".csv")
                    ),
                    
                    # Metadata file
                    fileInput(
                      inputId = "file2",
                      label = "Upload meta .csv",
                      accept = c(".csv")
                    )
                  ,
                    
                    # Gene List
                    textAreaInput(
                      inputId = "Anno.genes",
                      label = "Choose a subset of genes",
                      height = 100,
                    )),
                  # Start the heatmapping
                  actionButton("action2", "Start Ungrouped Heatmap",icon=icon("play-circle"),
                               style='padding:16px; font-size:125%'),
                  actionButton("action", "Start Grouped Heatmap",icon=icon("play-circle"),
                               style='padding:16px; font-size:125%'),
                  
                  h3("Annotation and Filtering Options"),
                    
                    wellPanel(style = "background:lightgreen;border-width:thick;border-color:black",
                    # Reactive dropdown to loaded metadata
                    selectInput(
                      'mydropdown',
                      label = 'Annotation or Grouping Variables',
                      choices = 'No choices here yet',
                      multiple = T
                    ),
                  
                  
                      # Filter 1
                      selectInput(
                        inputId='filter1',
                        label = 'Select if you want to filter for a trait',
                        choices = '',
                        multiple = F
                      ),
                      # Filter 2
                      selectInput(
                        'filter2',
                        label = "Filter for subjects with this trait value",
                        choices = 'No choices here yet',
                        multiple = F
                      ),
                      # Filter 3
                      selectInput(
                        'filter3',
                        label = "Go ahead and filter subjects by criteria?",
                        choices = c("YES", "NO"),
                        selected = "NO"
                      ),
                     # Checkbox for rownames
                      checkboxInput(inputId = 'Rows',
                                    label="Row Names?",
                                    value = FALSE, 
                                    width = NULL),
                    # Checkbox for colnames
                    checkboxInput(inputId="Col",
                                  label="Column Names?",
                                  value=F,
                                  width=NULL)),
                   
                   
                  h3("Style and Formatting"),
                  wellPanel(style = "background:#ffcccb;border-width:thick;border-color:black",
                     # Width of heatmap
                    setSliderColor(sliderId =c(1,2,3,4),color = c("Black","Black","Black","Black")),
                    sliderInput(inputId="groupWidth",
                                label="Heatmap Width",min = 10,
                                
                                max=1200,
                                value=1000
                                ),
                    # height of heatmap
                    sliderInput(inputId="Height",
                                label="Heatmap Height",
                                min=10,
                                max=500,
                                value=400,
                                
                                ),
                    
                    selectInput(inputId="Color",
                                label="Choose color palette",
                                choices=c("Default",COLORS),
                                selected = "Default"
                    ),
                    
                    sliderInput(inputId="Breaks",
                                label="Choose Color Palette to Adjust Color Breaks",
                               min=3,max=4,value=1)
                  ),
                    
                    # Download heatmap
                    downloadButton("downloadData2", "Download Ungrouped Heatmap",
                                   style='padding:16px; font-size:100%'),
                    downloadButton("downloadData", "Download Grouped Heatmap",
                                   style='padding:16px; font-size:100%')
                    
                    
                  ),
                  
                  mainPanel(tabsetPanel(
                    # Heatmap panel
                    tabPanel(
                      "Heatmaps",
                      br(),
                      h2("Ungrouped Heatmap"),
                      wellPanel(style = "background:white;border-width:thick;border-color:black", plotOutput("pheatmap.full"))
                      ,
                      br(),
                      h2("Grouped Heatmap"),
                      wellPanel(style = "background:white;border-width:thick;border-color:black", plotOutput("pheatmap"))
                      
                    ),
                    # Table of heatmap column values
                    tabPanel("Prac", DT::dataTableOutput("prac")),
                    tabPanel("Area", textOutput("Area"))
                  ))
                  
                ))


server <- function(input, output, session) {
  # Reactive counts variable
  
  starter <-
    reactive({
      read.csv(input$file1$datapath,
               header = T,
               row.names = 1)
    })
  
  # filter counts for selected genes in
  
  Genes <-
    reactive({
      readLines(textConnection(input$Anno.genes))
    })
  
  f.counts <-
    reactive({
      if (!input$Anno.genes == "") {
        starter()[rownames(starter()) %in% Genes(), ]
      } else{
        starter()
      }
    })
  
  # For colors
  col.pal <- reactive({if(input$Color=="Default"){colorRampPalette(rev(brewer.pal(n = input$Breaks, name =
        "RdYlBu")))(100)}else{colorRampPalette(rev(brewer.pal(n=input$Breaks
                                                ,name= input$Color)))(100)}
    })
  
  observeEvent(input$Color,{if(input$Color=="Default"){
    updateSliderInput(session, "Breaks", label = "Number of color breaks in heatmap",
                      min=3,
                      max=brewer.pal.info %>% filter(rownames(brewer.pal.info) %in% "RdYlBu") 
                      %>% select(maxcolors) %>% .$maxcolors,
                      value=7
    )}
    else{
    
    
    
  
              updateSliderInput(session, "Breaks", label = "Number of color breaks in heatmap",
                                min=3,
                                max=brewer.pal.info %>% filter(rownames(brewer.pal.info) %in% input$Color) 
                               %>% select(maxcolors) %>% .$maxcolors,
                                value=brewer.pal.info %>% filter(rownames(brewer.pal.info) %in% input$Color) 
                               %>% select(maxcolors) %>% .$maxcolors
                             )}
})           
    
  
  # Reactive metadata variable
  
  meta <-
    reactive({
      read.csv(input$file2$datapath,
               header = T,
               row.names = 1)
    })
  
  # Observe upload of meta file to update filter selection
  
  observeEvent(input$file2, {
    meta <- read.csv(input$file2$datapath,
                     header = T,
                     row.names = 1)
    
    updateSelectInput(session, "mydropdown", choices = colnames(meta))
    updateSelectInput(session, "filter1", choices = colnames(meta))
    
  })
  
  observeEvent(input$filter1, {
    updateSelectInput(session, "filter2", label = "Filter for subjects with this trait value",
                      if (input$filter1 == "") {
                        "None yet"
                      } else{
                        unique(meta() %>% select(input$filter1))
                      })
  })
  
  # Meta filtered for annotation track in full subject heatmap
  
  metaA <- reactive({
    meta() %>% select(input$mydropdown)
  })
  
  # Filter counts for specific metadata
  
  
  counts <- reactive({
    if (input$filter3 == "NO") {
      f.counts()
    } else{
      META.A <-
        meta() %>% select(input$filter1) %>% filter(.[1] == input$filter2)
      
      con <- f.counts()
      f.counts()[, colnames(con) %in% rownames(META.A)]
    }
  })
  
  # Action-Button for ungrouped heatmap
  UNGROUPED.HEAT <- eventReactive(input$action2, {
    counts()
  })
  


  
  # Heatmap columns function
  DATA <- eventReactive(input$action, {  
  
    metaZ <- meta() %>% select(input$mydropdown)
    
    total <-
      unite(metaZ, "total",sep=".") # Use for labelling grouped variables
    
    collapse <-
      paste(input$mydropdown, collapse = ".") # grouped variable IDs
    
    metaZZ <- cbind(metaZ, total)
    
    metaZZ$total <- paste(metaZZ$total, collapse, sep = ".")
    
  

        
    df <- data.frame(blank = "")
    
    names <- rownames(meta())
    
    split <-
      metaZZ %>% mutate(Names = names) # makes column for Subject IDs
    
    con <- counts()
    
    fun <-
      unique(metaZZ$total) # Iterate through each filtering and get means for each group variable
    
    
    Anno.Grouped <- data.frame()   # Annotation container
    
    A <- length(input$mydropdown)
  
   for (i in fun) {
      temp <- split %>% filter(total == i)
      
    
     x <- unlist(strsplit(i,"[.]"))
     y <- x[1:A]
     z <- x[(A+1):length(x)]
     Anno.Grouped <- rbind(Anno.Grouped,y)
     
      
      mine <-  con[, colnames(con) %in% temp$Names] # filter for correct subjects
      
      GroupMean <- rowMeans(mine) # get Row Means
      
      df <- cbind(df, GroupMean)
    }
   
    
    colnames(Anno.Grouped) <- z
    rownames(Anno.Grouped) <- fun
    
    df <- df[, -1]
    
    colnames(df) <- fun
    
    #Genes <- readLines(textConnection(input$Anno.genes))
    
    #df <- if(!input$Anno.genes==""){df[rownames(df) %in% Genes,]}else{df} # Filter dataframe for annotation genes
    
    df <- list(df,Anno.Grouped)
    df
    
  })
  
  # Functions for saving plots
    # Save grouped plot
  plotInput <-
    function() {
      pheatmap(
        DATA()[[1]],
        annotation_col=DATA()[[2]],
        cluster_rows = T,
        scale = "row",
        show_rownames = input$Rows,
        show_colnames = input$Col,
        color=col.pal()
        cluster_cols = T
      )
    }
    # Save ungrouped plot
  plotInput.full <- function() {
      pheatmap(
        UNGROUPED.HEAT(),
        cluster_rows = T,
        scale = "row",
        annotation_col=if(is.null(input$mydropdown)){NA}else{metaA()},
        show_rownames =input$Rows, 
        show_colnames = input$Col,
        color=col.pal(),
        cluster_cols = T
      )
    }
   
  # Grouped plot in UI
  output$pheatmap = renderPlot({
    pheatmap(
      DATA()[[1]],
      annotation_col = DATA()[[2]],
      cluster_rows = T,
      scale = "row",
      show_rownames = input$Rows,
      show_colnames = input$Col,
      color=col.pal(),
      angle_col = 45,
      cluster_cols = T
    )
    
  },width=function(){input$groupWidth},height=function(){input$Height})
  
  # Ungrouped plot in UI
  output$pheatmap.full = renderPlot({
    req(input$file1)
   # if () {
      pheatmap(
        UNGROUPED.HEAT(),
        cluster_rows = T,
        scale = "row",
        annotation_col= if(is.null(input$mydropdown)){NA}else{metaA()},
        show_rownames=input$Rows,
        show_colnames = input$Col,
        color = col.pal(),
        cluster_cols = T
      )
  },width=function(){input$groupWidth},height=function(){input$Height})
  
  
  output$prac <- DT::renderDataTable({
    DATA()[[2]]
  })
  
  
  
  output$Area <-
    renderText({
      readLines(textConnection(input$Anno.genes))
    })
  
  
  
  output$downloadData <-
    downloadHandler(
      filename = "Shinyplot.png",
      content = function(file) {
        png(file,height=input$Height,width=input$groupWidth) # can add height and width
        plotInput()
        dev.off()
      }
    )
  
  output$downloadData2 <-
    downloadHandler(
      filename = "Shinyplot.full.png",
      content = function(file) {
        png(file,height=input$Height,width=input$groupWidth) # can add height and width)
        plotInput.full()
        dev.off()
      }
    )
  
  
  
  
}







shinyApp(ui, server)




```
