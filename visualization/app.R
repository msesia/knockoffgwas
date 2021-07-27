suppressMessages(library(tidyverse))
suppressMessages(library(gridExtra))
suppressMessages(library(shiny))
#suppressMessages(library(dqshiny))
source("utils_clumping.R")
source("utils_plotting.R")
source("utils_shiny.R")
source("utils_manhattan.R")

data_dir <- "../data"
res_dir <- "../results"
lmm_dir <- "../data/lmm"
chromosomes <- 1:22
phenotypes <- c("example")

annotations <- load_annotations(data_dir)
genes <- unique(annotations$Exons.canonical$name2)

# create user interface
ui <- fluidPage(theme = "theme.css",

  # title
  headerPanel('KnockoffGWAS discovery view'),
  
  # side panel with inputs, error messages, and information box
  fluidRow(
    column(3,
      h4("Step 1: Select a phenotype."),
      wellPanel(
        selectInput(inputId = 'phenotype', label = 'Phenotype', choices = c(phenotypes)),    
        actionButton("load.results", "Load results")
      ),
      h4("Step 2: Type a chromosome or gene."),
      wellPanel(
        fluidRow( 
          column(6,
                 verticalLayout(uiOutput("ui_chr_select"),
                                actionButton("zoom.chr", "Zoom to chromosome")
                                )
                 ),
          column(6,
                    verticalLayout(uiOutput("ui_gene_select"),
                                   actionButton("zoom.gene", "Zoom to gene")
                                   )
                 )
        )
      ),
      h4("Step 3: Refine the locus."),
      wellPanel(
        fluidRow(
          column(6,actionButton("zoom.in", "Zoom in (slider)")),
          column(6,actionButton("zoom.out", "Zoom out (x10)"))          
        )
      ),
      actionButton("info", "?"),
      span(textOutput("message"), style="color:red"),
      absolutePanel(
        bottom = 0, left = 0, width = "24%",
        fixed = TRUE,        
        div(
          style="padding: 8px; border: 5px solid #CCC; background: #FFFFFF;", 
          HTML("<font size=\"3\">This website presents the results of the KnockoffGWAS methodology
            applied to a toy dataset. See [bioRxiv link] for the manuscript, and <a href=\"https://msesia.github.io/knockoffgwas\" target=\"_blank\"/>this webpage</a> 
            for more information.</font>")
        )
      )
  ),

  # main panel displaying results
  column(9,
    tabsetPanel(type = "tabs", id = "Tabset",
                tabPanel(title = "Genome", value = "manhattan", 
                         h4(textOutput("placeholder.manhattan")),
                         plotOutput('plot.manhattan', width = "100%", height = "700px")),
                tabPanel(title = "Locus", value = "manhattan.chr", 
                         h4(textOutput("placeholder.locus")),
                         splitLayout(width = "100%", cellWidths = c("4.6%", "76.1%", "19.3%"), 
                                     {}, uiOutput("slider"), {}),
                         plotOutput('plot.annotations', width = "100%", height = "700px")
                         )
    )
  )
)
)


# back-end code
server <- function(input, output, session) {

  output$placeholder.manhattan <- renderText({"[Select a phenotype to get started.]"})
  output$placeholder.locus <- renderText({"[Select a chromosome or gene to produce locus view.]"})
  # all parameters required to describe state of the app  
  state <- reactiveValues(phenotype = NULL,
                         association_results = NULL,
                         chr = NULL,
                         max.BP = NULL,
                         window.left = NULL, 
                         window.right = NULL,
                         slider.left = NULL,
                         slider.right = NULL,
                         highlight.gene = NULL)
  
  # what to do if "Info" button is pressed
  observeEvent(input$info, {
    if(input$Tabset == "manhattan"){
      box.title <- "Information on low-resolution results"
      box.message.1 <- "Top: Manhattan plot with BOLT-LMM p-values."
      box.message.2 <- "Bottom: Manhattan plot with KnockoffGWAS test statistics at low-resolution."
      box.message <- sprintf("%s<br>%s", box.message.1, box.message.2)
    } else{
      box.title <- "Information on high-resolution results"
      box.message.1 <- "Top: Manhattan plot with BOLT-LMM p-values."
      box.message.2 <- "Middle: Chicago plot with KnockoffGWAS test statistics at multiple resolutions."
      box.message.3 <- "Bottom: Functional annotations and gene positions."
      box.message <- sprintf("%s<br>%s<br>%s", box.message.1, box.message.2, box.message.3)
    }

    showModal(modalDialog(
      title = box.title,
      HTML(box.message)
    ))
  })
  
  # what to do if "Load association results" button is pressed
  observeEvent(input$load.results, {
    # clear error message
    output$message <- NULL
    # switch to the appropriate tab
    updateTabsetPanel(session, inputId = "Tabset", selected = "manhattan")
    # clear placeholder
    output$placeholder.manhattan <- NULL
    # check if phenotype has changed; if so, clear lower-level variables and 
    # load association results for new phenotype
    if(!is.null(state$phenotype)){
      if(input$phenotype != state$phenotype){
        state$chr <- NULL
        state$highlight.gene <- NULL
        state$max.BP <- NULL
        state$window.left <- NULL
        state$window.right <- NULL
        state$slider.left <- NULL
        state$slider.right <- NULL
        output$plot.annotations <- NULL # important: clear Chicago plot
        output$placeholder.locus <- renderText({"[Select a chromosome or gene to produce locus view.]"})
        state$phenotype <- input$phenotype 
        state$association_results <- load_association_results(res_dir, lmm_dir, input$phenotype)
      }
    } else{
        state$phenotype <- input$phenotype
        withProgress(message = 'Loading results...', value = 0, {
          state$association_results <- load_association_results(res_dir, lmm_dir, input$phenotype)
        })
    }
    # produce plot
    output$plot.manhattan <- renderPlot({
      withProgress(message = 'Rendering plot...', value = 0, {
        plot_manhattan_knockoffs(state$association_results$LMM,
                                 state$association_results$Pvalues,
                                 ytrans="identity")
      })
    })
  })

  # what to do if "Zoom to chromosome" button is pressed
  observeEvent(input$zoom.chr,{
      error <- TRUE
      # check if association data are loaded
      if(is.null(state$phenotype) | is.null(state$association_results)){
        output$message <- renderText({"Before clicking this button, first select a 
          phenotype and load association results."})
      } else{
        # check if valid chromosome number was entered
        chr <- as.integer(input$chr)
        if(is.na(chr) | is.null(chr)){
          error <- TRUE
        } else{
          if(!(chr %in% 1:22)){
            error <- TRUE
          } else{
            error <- FALSE
          }
        } 
        if(error){
          output$message <- renderText({"Type a chromosome number between 1 and 22."})
        } else{
          # clear error message
          output$message <- NULL
          # clear placeholder
          output$placeholder.locus <- NULL
          # clear highlighted gene if chromosome has changed
          if(!is.null(state$chr)){
            if(state$chr != chr){
              state$chr <- chr
              state$highlight.gene <- NULL
            }
          } else{
            state$chr <- chr
            state$highlight.gene <- NULL
          }
          # set window parameters to show whole chromosome
          chr.boundaries <- find_chr_boundaries(state$association_results, state$chr)
          state$min.BP <- chr.boundaries$min.BP
          state$max.BP <- chr.boundaries$max.BP
          state$chr <- chr
          state$window.left <- state$min.BP
          state$window.right <- state$max.BP
          state$slider.left <- state$window.left
          state$slider.right <- state$window.right
          updateTabsetPanel(session, inputId = "Tabset", selected = "manhattan.chr")
          output$plot.annotations <- renderPlot({
            withProgress(message = 'Rendering plot...', value = 0, {
              plot_combined_state(state, annotations)})
            })
        }
      }
  }
  )
  
  # what to do if "Zoom to gene" button is pressed
  observeEvent(input$zoom.gene,{
    # check if association data are loaded
    if(is.null(state$phenotype) | is.null(state$association_results)){
      output$message <- renderText({"Before clicking this button, first select a phenotype and
        load association results."})
    } else{
       if(input$gene %in% genes){ # check if valid gene name was entered
         withProgress(message = 'Finding location of gene...', value = 0, {
           # clear placeholder
           output$placeholder.locus <- NULL
           # clear error message
           output$message <- NULL
           # set chromosome appropriately
           filtered_exons <- filter(annotations$Exons.canonical, name2==input$gene)
           state$chr <- filtered_exons$chrom[1]
           state$max.BP <- max(filter(state$association_results$LMM, CHR==state$chr)$BP)
           # set center of gene to be center of window
           gene_min <- min(filtered_exons$txStart)
           gene_max <- max(filtered_exons$txEnd)
           window.center <- (gene_min + gene_max)/2
           # choose a window of width 1Mb
           state$window.left <- max(0, window.center - 0.25e6)
           state$window.right <- min(window.center + 0.25e6, state$max.BP)
           # adjust slider appropriately
           state$slider.left <- state$window.left
           state$slider.right <- state$window.right
           # set highlighted gene
           state$highlight.gene <- input$gene
           # switch to the appropriate tab
           updateTabsetPanel(session, inputId = "Tabset", selected = "manhattan.chr")
         })
         # produce Chicago plot
         output$plot.annotations <- renderPlot({
           withProgress(message = 'Rendering plot...', value = 0, {
             plot_combined_state(state, annotations)})
         })
       } else{
         output$message <- renderText({"Type a valid gene name."})
       }
    }
  })
  
  # what to do if "Zoom in" button is pressed
  observeEvent(input$zoom.in,{
    # check if association results are loaded
    if(is.null(state$phenotype) | is.null(state$association_results)){
      output$message <- renderText({"Before clicking this button, load association results for
        a phenotype and then choose a chromosome or gene."})
    } else{
      # check if window is chosen
      if(is.null(state$window.left) | is.null(state$window.right)){
        output$message <- renderText({"Before clicking this button, choose a chromosome or gene."})
      } else{
        # switch to appropriate tab
        updateTabsetPanel(session, inputId = "Tabset", selected = "manhattan.chr")
        output$message <- NULL
        # reset window based on slider
        state$window.left <- input$window[1]*1e6
        state$window.right <- input$window[2]*1e6
        state$slider.left <- state$window.left
        state$slider.right <- state$window.right
        # produce plot
        output$plot.annotations <- renderPlot({
          withProgress(message = 'Rendering plot...', value = 0, {
            plot_combined_state(state, annotations)})
        })
      }
    }
  })

  # what to do if "Zoom out" button is pressed
  observeEvent(input$zoom.out, {
    # check if association results are loaded
    if(is.null(state$phenotype) | is.null(state$association_results)){
      output$message <- renderText({"Before clicking this button, load association results for
        a phenotype and then choose a chromosome or gene."})
    } else{
      if(is.null(state$window.left) | is.null(state$window.right)){
        output$message <- renderText({"Before clicking this button, choose a chromosome or gene."})
      } else{
        # switch to appropriate tab
        updateTabsetPanel(session, inputId = "Tabset", selected = "manhattan.chr")
        # clear error message
        output$message <- NULL
        # set new window parameters
        window.center <- 0.5*(state$window.left + state$window.right)
        window.width <- 10*(state$window.right - state$window.left)
        state$window.left <- max(window.center - 0.5*window.width, state$min.BP)
        state$window.right <- min(window.center + 0.5*window.width, state$max.BP)
        # adjust slider
        state$slider.left <- state$window.left
        state$slider.right <- state$window.right
        # produce plot
        output$plot.annotations <- renderPlot({
          withProgress(message = 'Rendering plot...', value = 0, {
            plot_combined_state(state, annotations)})
        })
      }
    }
  })
  
  # produce slider UI element
  output$slider <- renderUI({
    # make sure window parameters are set
    if(!is.null(state$chr) & 
       !is.null(state$window.left) & 
       !is.null(state$window.right) & 
       !is.null(state$slider.left) &
       !is.null(state$slider.right)){
      
      # Make sure that the two ends of the slider do not overlap
      if(state$window.right<=state$window.left+0.2e6) {
        state.center <- max(0,(state$window.left+state$window.right)/2)
        state$window.left <- max(0,state.center-0.1e6)
        state$window.right <- state.center+0.1e6
        state$slider.left <- state$window.left
        state$slider.right <- state$window.right
      }
      
      # convert from BP to Mb
      window.left.Mb <- round(1e-6*state$window.left, 1)
      window.right.Mb <- round(1e-6*state$window.right, 1)
      slider.left.Mb <- round(1e-6*state$slider.left, 1)
      slider.right.Mb <- round(1e-6*state$slider.right, 1)
      
      # create slider
      sliderInput("window", label = NULL,
                  width = "99%", min = window.left.Mb, max = window.right.Mb,
                  step = min(0.5, window.right.Mb - window.left.Mb)/100,
                  value = c(slider.left.Mb, slider.right.Mb))
    } else{
      NULL
    }
  })
  
  # produce gene selection UI element
  output$ui_gene_select <- renderUI({
   autocomplete_input("gene", label = "Gene", value = state$highlight.gene,
                      options = genes, width = "90%")
  })
  
  # produce chromosome selection UI element
  output$ui_chr_select <- renderUI({textInput(inputId = 'chr', width = "90%", 
                                           value = state$chr, label = 'Chromosome')
  })
}

shinyApp(ui = ui, server = server)
