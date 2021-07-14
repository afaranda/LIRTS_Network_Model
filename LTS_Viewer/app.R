library(shiny)
library(dplyr)
library(pheatmap)
setwd("~/Documents/LEC_Time_Series/LTS_Viewer")
source("load_global_wildtype.R")


gene_id_symbol<-gene_meta %>%
  filter(gene_id %in% fpkm$gene_id & ! duplicated(SYMBOL)) %>%
  select(gene_id, SYMBOL)
genes <- gene_id_symbol$gene_id
names(genes) <- gene_id_symbol$SYMBOL

ui<-fluidPage(
  titlePanel("LEC Time Series Viewer"),
  sidebarLayout(
    sidebarPanel(
      h4("Plot time-series for a gene"),
      p(
        paste(
          "Use this app to visualize changes in a gene's expression in",
          "lens epithelial cells over time after injury."
        )
      ),
      fluidRow(
        column(6,
               h5("Select the gene of interest"),
               selectInput("gene", "Gene", choices = names(genes))
        ),
        column(4,
               h5("Box plots or individual FPKM measurements"),
               radioButtons("ptype", "Plot Type", choices = c("Box","Scatter")))
      ),
      
      h4("Plot target clusters for a transcription factor"), 
      fluidRow(
        column(6,
               h5("Select the transcription factor of interest"),
               selectInput("factor", "Transcription Factor", choices = tf_tgt$Factor)
        ),
        column(3,
               h5("Log Transform"),
               radioButtons("uselog", "", choices = c("Yes","No"))
        ),
        column(3,
               h5("Scale and center"),
               radioButtons("stdize", "", choices = c("Yes","No")))
        
      ), width = 8
    ),
    mainPanel(
      tabsetPanel(
        type="tabs",
        tabPanel("Plot", plotOutput("timeSeries",height = "300px")),
        tabPanel("Table", tableOutput("tsTable")),
        tabPanel("TF Cluster", plotOutput("tfClust", height="600px")),
        tabPanel("Example Network", plotOutput("graphExample"))
      )
    )
  )
)

server <- function(input, output){
  
 
  dataInput <- reactive({
    g <- genes[input$gene]
    plt %>% filter(gene_id == g)
  })
  
  clustInput <- reactive({
    x<-tf_tgt%>% filter(Factor==input$factor) %>% as.data.frame()
    
    x<-x %>% inner_join(
      plt, by=c(tgt_gid="gene_id")
    ) %>% select(sample, Target, fpkm)

    x <- x %>%
      mutate(
        fpkm = switch(
          input$uselog,
          Yes=log(fpkm + min(x %>% filter(fpkm > 0) %>% pull(fpkm)),2),
          No=fpkm
        )
      ) %>% 
      group_by(Target) %>%
      mutate(
        fpkm = switch(
          input$stdize,
          Yes=(fpkm - mean(fpkm)),
          No = fpkm
        )
      ) %>%
      mutate(
        fpkm = switch(
          input$stdize,
          Yes=(fpkm / sd(fpkm)),
          No = fpkm
        )
      ) %>%
      pivot_wider(id_cols=Target, names_from=sample, values_from=fpkm) %>%
      as.data.frame()
    
    row.names(x)<-x$Target
    x<-as.data.frame(x[,setdiff(1:ncol(x), grep("Target", names(x)))])
  })
  
  finalInput <- reactive({
    switch(
      input$ptype,
      Box=ggplot(dataInput(), aes(x=hours_pcs, y=fpkm, group=hours_pcs)) +
        geom_boxplot() + ggtitle(input$gene) + theme(
          plot.title = element_text(size = 25),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15)
        ),
      
      Scatter=ggplot(dataInput(), aes(x=hours_pcs, y=fpkm, group=hours_pcs)) + 
        geom_point(aes(color=batch), size=5) + ggtitle(input$gene) + theme(
          plot.title = element_text(size = 25),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15)
        )
    )
    
  })
  
  # output$selectedGene <- renderText({
  #   paste("Selected:", input$gene, genes[input$gene])
  # })
  # 
  output$timeSeries <- renderPlot({
    finalInput()
  })
  
  output$tsTable <- renderTable({
    dataInput() %>% select(
      sample, batch, hours_pcs, fpkm
    ) %>% arrange(hours_pcs, batch)
  })
  
  output$tfClust<-renderPlot({
    pheatmap(clustInput(), annotation_col=annot[2])
  })
  
  output$graphExample<-renderPlot({
    library(igraph)
    g<-graph(
      edges=c(
        'Gata3', 'Nfkb1',
        'Gata3', 'Rela',
        'Gata3', 'Cxcl1',
        'Gata3', 'Cdh1',
        'Gata3', 'Gata3',
        'Gata3', 'Runx3',
        'Runx3', 'Runx1', 
        'Runx3', 'Thsb1',
        'Runx3', 'Timp1',
        'Runx3', 'Wnt4'
      )
    )
    plot(g)
  })
  # 
  
  # output$selectedRange <- renderText({
  #   paste("Ranging from", input$range[1], "to", input$range[2])
  # })
}

shinyApp(ui=ui, server=server)
