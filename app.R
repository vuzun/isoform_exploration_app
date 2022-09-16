library(shiny)
library(tidyverse)
library(Seurat)
library(schex) # if missing, install from bioconductor: BiocManager::install("schex")


ui <- fluidPage(
    
    #titlePanel("Single gene inspection"),
    navbarPage("Isoform exploration",
    # sidebar with inputs
    tabPanel("Single gene check",
    fluidRow(
        column(width=5, #tags$style(".well {background-color: white;}"),
            
            textInput("gene",
                      "Gene",
                      "Enter gene present in the data..."
            ),
            actionButton("getbutton","Get"),
            plotOutput("distPlot")
        ),
        
        
        # plots 
        column(width=7,

                tableOutput("gene_info"),
                plotOutput("agePlot")
    
            )
        )
    ),
    tabPanel("UMAP",
             fluidRow(
                 column(width = 6,
                        plotOutput("featurePlot")),
                 column(width = 6,
                        selectInput(inputId = "metaselect",label = "Select metadata variable:",
                                    choices = c("Cluster"="seurat_clusters",
                                                "Age"="group")),
                        plotOutput("clusterPlot"))
             )
    ),
    tabPanel("Most changed isoforms",
             sidebarLayout(
                sidebarPanel(
                    numericInput('threshold_counts', 'Minimum number of counts', 30,
                                 min = 1, max = 1000),
                    radioButtons("change_dir", "Change direction: ",
                                 c("Old/Young" = "fromY",
                                   "Young/Old" = "fromO")),
                    downloadButton("downloadData", "Download")
                ),
                mainPanel(tableOutput("changeTable")))
             )

    )
)

# server logic
server <- function(input, output, session) {
    
    withProgress(message =  "Preparing count data...", {
    
    counts <- readRDS("b_counts.rds")
    num_cells <- ncol(counts)
    incProgress(0.3)
    counts$total_counts <- rowSums(counts)
    incProgress(0.1)
    
    # 3 new columns - tx, gene, total counts per the corresponding gene
    counts <- counts %>% 
        mutate(tx_name=rownames(.)) %>% 
        mutate(gene_name=tx_name %>% sapply(function(n) substr(n,1,nchar(n)-4))) %>% 
        group_by(gene_name) %>% 
        mutate(total_per_gene=sum(total_counts)) %>% 
        ungroup()
    
    })
    
    withProgress(message = "Preparing 10X object...", min = 0.2, {
        sobj_10x <- readRDS("schex40_sobj_10x_RNA.rds")#
    })
    
    
    # this part first to break on new data probably
    # young(1,3,5), old (2,4,6) tx sums
    ids_samples <- colnames(counts)[1:num_cells] %>% 
        gsub("[[:alpha:]]|_","",.) %>% as.numeric
    inds_young <- which(ids_samples%%2==1)
    inds_old <- which(ids_samples%%2==0)
    counts$total_young <- counts[,inds_young] %>% rowSums()
    counts$total_old <- counts[,inds_old] %>% rowSums()

    # this probably needs some lower count threshold
    # e.g. look at Slc43a3
    counts <- counts %>% group_by(gene_name) %>% 
            mutate(per_gene_old=sum(total_old)) %>% 
            mutate(per_gene_young=sum(total_young)) %>% 
            mutate(tx_ageing_change= (total_old*per_gene_young)/(total_young*per_gene_old)) %>% ungroup()
        
        

    
    
    extract_iso <- function(gene="Mpl", lr_df=counts, ageing_exp=TRUE){
        gene <- paste0("^",gene,"$")
        fetch_columns <- c("tx_name","total_counts","total_per_gene")
        if(ageing_exp){
            fetch_columns <- c(fetch_columns,"total_young","total_old",
                               "tx_ageing_change")
        }
        
        lr_df <- lr_df[grep(gene, lr_df$gene_name, ignore.case = T),fetch_columns]
        colnames(lr_df) <- c("Name","Counts","Gene_level", "Young", "Old",
                             "Age_change")
        
        return(lr_df)
    }
    

    
    iso_data <- reactive(extract_iso(input$gene))
    
    piechart_gen <- eventReactive(input$getbutton, {
        x    <- iso_data()
        pie_step <- x$Counts %>% sum %>% `/`(10)
        x %>% ggplot(aes(x="",y=Counts,fill=Name)) +
            geom_bar(width=1,stat="identity")+ coord_polar("y", start=0) +
            labs(x="") + scale_fill_discrete(name = "Isoform name") +
            scale_y_continuous(breaks=seq(0,pie_step*10,pie_step))
    })
    
    agebar_gen <- eventReactive(input$getbutton,{
        iso_data() %>% ungroup() %>%  select(Name, Old, Young) %>% 
            gather(age, counts, -Name) %>% 
            ggplot(aes(x=as.factor(age), y=counts, fill=Name)) +
            geom_bar(position="fill",stat="identity") + labs(x="", y="Scaled counts")
    })
    
    featureplot_gen <- eventReactive(input$getbutton,{
        plot_hexbin_gene(sobj_10x, "counts", input$gene, "mean")
        #FeaturePlot(sobj_10x, input$gene, cols = c("grey90","darkred"))
    }) 
    metaplot_gen <- eventReactive(input$metaselect,{
        plot_hexbin_meta(sobj_10x, col = input$metaselect, "majority")
        #DimPlot(sobj_10x, reduction = "umap")  
    })
    
    changetable_gen <- reactive({
        
        direction <- 1
        if(input$change_dir=="fromO"){direction <- (-1)}
        
        counts  %>% filter(is.finite(tx_ageing_change)) %>%
            filter(total_counts>input$threshold_counts) %>% 
            select(tx_name, total_counts, tx_ageing_change) %>%
            mutate(tx_ageing_change=(tx_ageing_change**direction)) %>% 
            arrange(desc(tx_ageing_change)) %>% head(20)

        #Hnrnpa2b1
    })
    
    observe({
        metacol <- input$metaselect
        metacols <- colnames(sobj_10x@meta.data)
        updateSelectInput(session,"metaselect",
                          choices=metacols[!grepl("ID",metacols)],
                          selected = metacol)
    })
    
    output$distPlot <- renderPlot({
        piechart_gen()
    })
    
    output$agePlot <- renderPlot({
        agebar_gen()
    })
    
    output$featurePlot <- renderPlot({
        featureplot_gen()  
    })
    
    output$clusterPlot <- renderPlot({
        metaplot_gen()
    })
    
    output$gene_info <- renderTable({
        iso_data()
    })
    
    output$changeTable <- renderTable({
        changetable_gen()
    })
    
    output$downloadData <- downloadHandler(
        filename = function(){
            paste0("Transcripts_by_ageing_change_",
                   input$threshold_counts,
                   "counts_min.csv",collapse = "")
            },
        content = function(file){
            write.csv(changetable_gen(), file, row.names = F)
        })
}

shinyApp(ui = ui, server = server)

