### Library loading
	cran_packages <- c("shinyjs", "ggplot2", "ggrepel", "plotly", "randomForest", "DT", "bslib")
	bioc_packages <- c("DESeq2")
	install_if_missing_cran <- function(pkg) {
 		if (!requireNamespace(pkg, quietly = TRUE)) {
    			install.packages(pkg, dependencies = TRUE, quiet = TRUE)
		}
	}
	install_if_missing_bioc <- function(pkg) {
  		if (!requireNamespace(pkg, quietly = TRUE)) {
    			if (!requireNamespace("BiocManager", quietly = TRUE)) {
				install.packages("BiocManager", dependencies = TRUE, quiet = TRUE)
			}
			BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
		}
	}
	invisible(lapply(cran_packages, install_if_missing_cran))
	invisible(lapply(bioc_packages, install_if_missing_bioc))
	library(shiny)
	library(shinyjs)
	library(ggplot2)
	library(scico)
	library(ggrepel)
	library(plotly)
	library(DESeq2)
	library(randomForest)
	library(DT)
	library(bslib)
	result_env <- new.env()

options(shiny.maxRequestSize = 200*1024^2)  # 100 MB


# === CUSTOM THEME ===
theme_custom <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  base_font = font_google("Fira Sans"),
  heading_font = font_google("Fira Sans Condensed"),
  primary = "#5b9bd5",
  success = "#28a745",
  info = "#17a2b8",
  warning = "#ffc107",
  danger = "#dc3545"
)

# === UI ===
ui <- fluidPage(
  useShinyjs(),
  theme = theme_custom,
  tags$head(
    tags$style(HTML("body {background-color: #dbe5ff !important;}",".footer {display: flex; justify-content: space-between; align-items: center; margin-top: 40px; padding: 20px; font-weight: bold; font-size: 14px;} .about-box {margin-top: 10px; margin-bottom: 30px; padding: 20px; border: 1px solid #ddd; border-radius: 12px; background-color: #fdfdfd; box-shadow: 0 2px 6px rgba(0,0,0,0.05);} .dataTables_wrapper td {white-space: normal !important; word-break: break-word;} .card {margin-top: 20px;} .card-title {font-weight: bold; font-size: 18px;} .section-title {margin-top: 30px; font-size: 26px; font-weight: bold; border-bottom: 2px solid #eee; padding-bottom: 5px;}"))
  ),
  
  fluidRow(
    column(12, align = "center",
           tags$img(src = "images/logo.png", style = "height: 250px; margin-top: 0px; margin-bottom: 0px; max-width: 500px; width: 80%")
    )
  ),
  
  fluidRow(
    column(8, offset = 2,
           div(class = "about-box",
	   style = "margin-top: 0px;",
               HTML("<h4><i class='fa-solid fa-dna'></i> About TASC:</h4> TASC (T-ALL Subtype Classifier) is a user-friendly R Shiny application designed to classify T-cell acute lymphoblastic leukemia samples into molecular subtypes using RNA-seq count data. It leverages a Random Forest model trained on curated expression signatures to provide subtype predictions and feature visualizations. Visit our <a href='https://github.com/CBenetti/TASC/blob/main/wiki.md' target='_blank'>GitHub repository</a> to know more about TASC usage or to contact us.")
           )
    )
  ),
  
  fluidRow(
    column(6, offset = 3,
           card(
             card_header("Upload Input Data", class = "bg-primary text-white"),
             card_body(
             fluidRow(
                 column(6,
                        selectInput("inputSource", "Source of count matrix:",
                                    choices = c("Select..." = "", "Salmon" = "salmon", "featureCounts" = "featureCounts", "Custom format"="custom"),
                                    width = "100%")
                 ),
                 	column(6,
                              uiOutput("file_input_ui")
			      #fileInput("countsfile", "Upload count matrix (.tsv):", accept = c(".tsv", ".txt"), width = "100%",id = "upload_input")
		  )
                 ),
		fluidRow(
               	  column(6, offset = 3,
                      actionButton("test_button", "Run Test on Example Data", 
                                   class = "btn btn-outline-primary", 
                                   style = "width: 100%; margin-top: 1px;"))
               ),
               tags$div(id = "info_source_help", style = "font-size: 13px; background-color: #eef7fb; padding: 8px; border-radius: 5px; margin-top: 10px;",
                        HTML("<b>Info:</b> Please choose a valid count matrix source and upload a file with gene IDs (starting with <code>ENSG</code>) in the first column.")),
               textOutput("gene_id_error"),
               textOutput("input_format_error")
             )
           )
    )
  ),
  
  div(id = "waiting_message",
      style = "text-align: center; font-size: 20px; padding: 50px; color: #555;",
      "\U0001F4E5 Please upload data to display the results."
  ),
  
  hidden(div(id = "results_section",
             fluidRow(
               column(12,
                      hr(),
                      uiOutput("Imputation_alert"),
                      downloadButton("downloadGenes", "Download Imputed Gene List")
               )
             ),
		br(),
             fluidRow(
               column(6,
                      div(class = "section-title", "Predicted Classes"),
                      card(
                        card_body(DTOutput("classTable"))
                      )
               ),
               column(6,
                      div(class = "section-title", "PCA Plot"),
                      card(
                        card_body(
                          plotlyOutput("pcaPlot", height = "500px"),
                          br(),
                          downloadButton("downloadPlot", "Download PCA Plot")
                        )
                      )
               )
             ),
     
	fluidRow(
               column(12,
                      div(class = "section-title", "Metagene Feature Plots"),
                      uiOutput("metaPlots"),
                      downloadButton("downloadMetaPlots", "Download Metagene Plots (PDF)")
               )
             ),
	fluidRow(
  column(12,
    div(class = "section-title", "Explore Signature on Demand"),
    selectInput("manual_sig", "Choose a Signature Class:",
                choices = c("Select a class..." = ""), width = "50%"),
    conditionalPanel(
      condition = "input.manual_sig != ''",
      plotOutput("manualMetaPlot", height = "800px"),
      br(),
    downloadButton("downloadManualMetaPlot", "Download Metagene Plot")	
    )
  )
),

             fluidRow(
               column(12,
                      downloadButton("downloadProbabilities", "Download the complete table of probabilities")
               )
             )
  )),
  
  tags$div(class = "footer",
           tags$div(HTML("Authors: <b>Cinzia Benetti</b>, <b>Omar Almolla</b>, <b>Francesco Boccalatte</b>")),
           tags$div(tags$img(src = "images/irccs_logo.png", height = "150px"))
  )
)

# === SERVER ===
server <- function(input, output, session) {
	shinyjs::useShinyjs()
	output$file_input_ui <- renderUI({
	  fileInput("countsfile", "Upload count matrix (.tsv):", accept = c(".tsv", ".txt"), width = "100%")
	})
	load("data/unique_signatures.rda")
	test_mode <- reactiveVal(FALSE)
	counts_data <- reactiveVal(NULL)
	custom_data <- reactiveVal(NULL)

  observeEvent(result(), {
    shinyjs::hide("waiting_message")
    shinyjs::show("results_section")
    updateSelectInput(session, "manual_sig",
     choices = c("Select a class..." = "",gsub("B:SAMP:subtype_","",sort(unique(unlist(unique_signatures[, 3]))[-1])))
)

  })
  
  observeEvent(input$test_button, {
	
    shinyjs::hide("waiting_message")
    shinyjs::show("results_section")
    updateSelectInput(session, "manual_sig",
    choices = c("Select a class..." = "",gsub("B:SAMP:subtype_","",sort(unique(unlist(unique_signatures[, 3]))[-1]))))
    test_mode(TRUE)
    counts <- read.table("test/test.tsv", header = TRUE, sep = "\t", check.names = FALSE)
    gene_ids <- counts[, 1]
    numeric_counts <- counts[, -(1:6)]
    counts_data(list(counts = numeric_counts, gene_ids = gene_ids))
  })

  observeEvent(input$countsfile, {
    shinyjs::enable("countsfile")
    shinyjs::enable("inputSource")
    test_mode(FALSE)
    req(input$inputSource != "", input$countsfile)
    counts <- read.table(input$countsfile$datapath, header = TRUE, sep = "\t", check.names = FALSE)
    gene_ids <- counts[, 1]
    if (!all(grepl("^ENSG", gene_ids))) {
      output$gene_id_error <- renderText({ "❌ gene_id should be in Ensembl format (e.g., ENSG...)" })
      return(NULL)
    } else {
      output$gene_id_error <- renderText({ "" })
    }
    numeric_counts <- switch(input$inputSource,
                             "featureCounts" = counts[, -(1:6)],
                             "salmon" = counts[,-c(1:2)],
                             "custom" = counts[, -1])
    if(is.numeric(numeric_counts[,1])==FALSE){output$input_format_error <- renderText({ "The format of the input counts is not correct, please check documentation and try again." })
      return(NULL)}else{output$input_format_error <- renderText({ "" })}
    custom_data(list(counts = numeric_counts, gene_ids = gene_ids))
  })

   
  result <- reactive({
    if (test_mode()) {
      pc <- counts_data()
	test_mode <- reactiveVal(FALSE)
        counts_data <- reactiveVal(NULL)
    }else{
     pc <- custom_data()
	custom_data <- reactiveVal(NULL)
    }   
    req(pc)
    counts <- pc$counts
    gene_ids <- pc$gene_ids
    rm(pc)
    load("data/random_forest_model.rda")
    load("data/unique_signatures.rda")
withProgress(message = "Processing input file...", value = 0, {
    DGE_OBJ <- DESeqDataSetFromMatrix(round(counts), design = ~0,
                                      rowData = data.frame(gene_id = gene_ids),
                                      colData = data.frame(samples = colnames(counts)))
    vst_data <- assay(vst(DGE_OBJ, blind = TRUE))
    rownames(vst_data) <- sub("\\..*", "", gene_ids)
   colnames(vst_data) <- DGE_OBJ$samples
})    
    gen <- sub("\\..*", "", model$subset_genes)
    missing_genes <- gen[!gen %in% rownames(vst_data)]
    
    if (length(missing_genes) > 0) {
withProgress(message = "Performing gene imputation...", value = 0, {    
      idx <- which(!gen %in% rownames(vst_data))
      sd_imput <- apply(model$original_data_matrix[model$subset_genes[idx], ], 1, sd)
      mean_imput <- apply(model$original_data_matrix[model$subset_genes[idx], ], 1, mean)
      sub <- t(sapply(seq_along(mean_imput), function(i) rnorm(ncol(vst_data), mean_imput[i], sd_imput[i])))
      rownames(sub) <- missing_genes
      colnames(sub) <- colnames(vst_data)
      vst_sub <- rbind(vst_data, sub)[gen, ]
      rownames(vst_sub) <- model$subset_genes
})
    } else {
      vst_sub <- vst_data[gen, ]
      rownames(vst_sub) <- model$subset_genes
    }

  withProgress(message = "Performing predictions...", value = 0, {  
    class_tmp <- predict(model$rf, newdata = t(vst_sub))
    class_p <- predict(model$rf, newdata = t(vst_sub),type="prob")
    colnames(class_p)<-gsub("B:SAMP:subtype_","",colnames(class_p))
    class <- class_tmp
    class[which(class_p[,grep("LMO2",colnames(class_p))]>model$threshold & class=="B:SAMP:subtype_ETP-like")]<-"B:SAMP:subtype_LMO2 γδ-like"
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999933", "#CC6666",
                   "#339933", "#9966CC", "#66CCEE")
    names(cbPalette) <- gsub("B:SAMP:subtype_", "", levels(class))
    
    pca_coord <- as.data.frame(prcomp(t(vst_data))$x)
    pca_coord$class <- gsub("B:SAMP:subtype_", "", as.vector(class))
  })  
    p <- ggplot(pca_coord, aes(PC1, PC2, color = class, label = rownames(pca_coord))) +
      geom_point(size = 3) +
      geom_label_repel(color = "black", show.legend = FALSE) +
      scale_color_manual(values = cbPalette) +
      theme_minimal(base_size = 14) +
      labs(
        x = paste0("PC1 (", round(summary(prcomp(t(vst_data)))$importance[2, 1] * 100, 1), "%)"),
        y = paste0("PC2 (", round(summary(prcomp(t(vst_data)))$importance[2, 2] * 100, 1), "%)"),
        title = "PCA of subtype prediction results"
      )+theme(plot.margin = grid::unit(c(35, 5, 5, 5), "pt"))
    
    unique_signatures[,3] <- gsub("B:SAMP:subtype_", "", unique_signatures[,3])
    meta_plots <- list()
    for (sig in unique(pca_coord$class)) {
      genes <- unique_signatures[unique_signatures[, 3] == sig, 2]
      metagenes <- apply(vst_data, 2, function(x) mean(x[na.omit(match(genes, rownames(vst_data)))]))
      pca_temp <- pca_coord
      pca_temp$metagene <- metagenes
      class_corrected <- pca_coord$class
      names(class_corrected) <- rownames(pca_coord)
	if(length(which(genes%in%gene_ids==T))>(length(genes)/2)){
	 meta_plots[[sig]] <- ggplot(pca_temp, aes(PC1, PC2, color = metagene, label = rownames(pca_temp))) +
        geom_point(size = 3) +
        geom_text_repel(
          color = "black",
          size = 4,
          box.padding = unit(1, "lines"),
          point.padding = unit(1.5, "lines"),
          force = 10,                      # forza di repulsione molto elevata
          force_pull = 0.1,                # poco effetto di attrazione verso il punto
          max.overlaps = Inf,
          segment.color = "grey50",
          segment.curvature = -0.1,
          segment.ncp = 3,
          nudge_y = 15,                    # sposta le etichette verso l'alto
          show.legend = FALSE
        )  +
	scale_color_scico(palette="buda",midpoint=mean(pca_temp$metagene),name=paste(sig, "metagene"))+
        coord_fixed() +
	theme_minimal(base_size=14)+
        theme(plot.title = element_text(size = 20, face = "bold"),panel.border = element_rect(color = "black", fill = NA, size = 1)) +
        labs(title = paste("Feature plot of", sig, "signature"),
             x = paste0("PC1 (", round(summary(prcomp(t(vst_data)))$importance[2, 1] * 100, 1), "%)"),
             y = paste0("PC2 (", round(summary(prcomp(t(vst_data)))$importance[2, 2] * 100, 1), "%)"))+facet_wrap(vars(class))+
	theme(plot.margin = grid::unit(c(30, 10, 30, 10), "pt"))
	}else{meta_plots[[sig]]<-sig}
    }
	names(meta_plots)[grep("LMO2",names(meta_plots))]<-"LMO2"
	names(meta_plots)[grep("αβ",names(meta_plots))]<-"TAL1_ab"
	result_env$vst_data <- vst_data
	result_env$pca_coord <- pca_coord
	result_env$gene_ids <- rownames(vst_data)


    structure(list(
      class = class_corrected,
      p = p,
      prob = class_p,
      imp = missing_genes,
      plot_data = p,
      meta_plots = meta_plots
    ), model_genes = model$subset_genes)
  })
  
#  output$classTable <- renderDT({
#    req(result())
#    df <- data.frame(Sample = names(result()$class), Prediction = result()$class)
#    datatable(df, extensions = 'Buttons',
#              options = list(dom = 'Blfrtip', buttons = list('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 15),
#              rownames = FALSE)
#  })
 
output$classTable <- renderDT({
  req(result())
  df <- data.frame(Sample = names(result()$class), Prediction = result()$class)
  datatable(df,
    extensions = 'Buttons',
    options = list(
      dom = 'Blfrtip',
      buttons = list(
        list(
          extend = 'copy',
          text = 'Copy'
        ),
        list(
          extend = 'csv',
          text = 'CSV',
          filename = JS("function() {
            let d = new Date();
            let timestamp = d.toISOString().replace(/[:.]/g, '-');
            return 'TASC_class_predictions_' + timestamp;
          }")
        ),
        list(
          extend = 'excel',
          text = 'Excel',
          filename = JS("function() {
            let d = new Date();
            let timestamp = d.toISOString().replace(/[:.]/g, '-');
            return 'TASC_class_predictions_' + timestamp;
          }")
        ),
        list(
          extend = 'pdf',
          text = 'PDF',
          filename = JS("function() {
            let d = new Date();
            let timestamp = d.toISOString().replace(/[:.]/g, '-');
            return 'TASC_class_predictions_' + timestamp;
          }")
        ),
        list(
          extend = 'print',
          text = 'Print'
        )
      ),
      pageLength = 15
    ),
    rownames = FALSE
  )
})
 
  output$pcaPlot <- renderPlotly({ ggplotly(result()$p) })
  
  output$Imputation_alert <- renderUI({
    res <- result()
    readRDS("data/top_30.rds")->top30
    imp_genes <- res$imp
    model_genes <- attr(res, "model_genes")
    if (is.null(model_genes)) return(NULL)
    total_genes <- length(model_genes)
    imp_ratio <- length(imp_genes) / total_genes
    if (length(imp_genes) > 0) {
      tagList(
        tags$p(strong(paste("Imputation required for", length(imp_genes), "genes:"))),
        verbatimTextOutput("imputedGenes"),
        if (imp_ratio > 0.10) {
          tags$p("⚠️ More than 10% of genes were imputed using model-based mean & std dev",
                 style = "color: red; font-weight: bold; font-size: 16px;")
        },
	if (length(which(imp_genes%in%top30)) > 0) {
        tags$p("⚠️ Some of the imputed genes are among the top 30 important features!",
               style = "color: darkorange; font-weight: bold; font-size: 16px;")
        }
      )
    }
  })
  
  output$imputedGenes <- renderPrint({ result()$imp })
  
output$metaPlots <- renderUI({
    req(result()$meta_plots)          # do nothing until meta_plots is non-NULL
    pl <- result()$meta_plots         # now safe to grab the list

    # build one UI element per entry
    ui_list <- lapply(names(pl), function(name) {
      output_id <- paste0("meta_plot_", name)
      plot_obj  <- pl[[name]]

      if (inherits(plot_obj, "ggplot")) {
        column(
          width = 12, align = "center",
          plotOutput(output_id, height = "800px")
        )
      } else {
        column(
          width = 12, align = "center",
          div(
            class = "alert alert-warning",
            paste0(
              "Impossible to calculate metagene for “", name,
              "”: not enough signature genes found in your dataset."
            )
          )
        )
      }
    })

    do.call(tagList, ui_list)
  })


  # 2) observeEvent on the same meta_plots—register one renderPlot per ggplot
  observeEvent(result()$meta_plots, {
    pl <- result()$meta_plots

    # loop and assign each to output[[...]] in its own local() closure
    for (name in names(pl)) {
      local({
        nm   <- name
        p_obj <- pl[[nm]]
        out_id <- paste0("meta_plot_", nm)

        output[[out_id]] <- renderPlot({
          p_obj
        })
      })
    }
  })  

  output$downloadPred <- downloadHandler(
    filename = function() paste0("class_predictions_", Sys.Date(), ".tsv"),
    content = function(file) {
      write.table(data.frame(Sample = names(result()$class), Prediction = result()$class),
                  file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  output$downloadProbabilities <- downloadHandler(
    filename = function() paste0("full_probability_matrix", Sys.Date(), ".tsv"),
    content = function(file) {
      write.table(as.data.frame(result()$prob),
                  file, sep = "\t", row.names = TRUE, col.names=NA, quote = FALSE)
    }
  )
  
  output$downloadPlot <- downloadHandler(
    filename = function() paste0("PCA_plot_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file); print(result()$plot_data); dev.off()
    }
  )
  
  output$downloadGenes <- downloadHandler(
    filename = function() paste0("imputed_genes_", Sys.Date(), ".txt"),
    content = function(file) writeLines(result()$imp, file)
  )
  
  output$downloadMetaPlots <- downloadHandler(
    filename = function() paste0("metagene_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 6)
      lapply(result()$meta_plots, function(x){if(is.character(x)){NULL}else{print(x)}})
      dev.off()
    }
  )
  output$manualMetaUI <- renderUI({
  req(input$manual_sig)
  plotname <- paste0("manual_plot_", input$manual_sig)
  uiOutput(plotname)
  })

output$manualMetaPlot <- renderPlot({
  req(result(), input$manual_sig)

  sig <- paste0("B:SAMP:subtype_", input$manual_sig)
  vst_data <- result_env$vst_data
  pca_coord <- result_env$pca_coord
  gene_ids <- result_env$gene_ids

  genes <- unique_signatures[unique_signatures[, 3] == sig, 2]
  found_genes <- which(genes %in% gene_ids)

  if (length(found_genes) <= (length(genes) * 2 / 3)) {
    plot.new()
    text(0.5, 0.5, paste("⚠️ Not enough signature genes for", input$manual_sig), cex = 1.5)
    return()
  }

  metagenes <- apply(vst_data, 2, function(x) mean(x[na.omit(match(genes, rownames(vst_data)))]))
  pca_temp <- pca_coord
  pca_temp$metagene <- metagenes

  ggplot(pca_temp, aes(PC1, PC2, color = metagene, label = rownames(pca_temp))) +
    geom_point(size = 3) +
    geom_text_repel(
      color = "black", size = 4, box.padding = unit(1, "lines"),
      point.padding = unit(1.5, "lines"), force = 10, force_pull = 0.1,
      max.overlaps = Inf, segment.color = "grey50", segment.curvature = -0.1,
      segment.ncp = 3, nudge_y = 15, show.legend = FALSE
    ) +
    scale_color_scico(palette = "buda", midpoint = mean(pca_temp$metagene),
                      name = paste(input$manual_sig, "metagene")) +
    coord_fixed() +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.margin = grid::unit(c(30, 10, 30, 10), "pt")
    ) +
    labs(
      title = paste("Feature plot of", input$manual_sig, "signature"),
      x = paste0("PC1 (", round(summary(prcomp(t(vst_data)))$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(prcomp(t(vst_data)))$importance[2, 2] * 100, 1), "%)")
    ) +
    facet_wrap(vars(class))
})
output$downloadManualMetaPlot <- downloadHandler(
  filename = function() {
    paste0("metagene_plot_", input$manual_sig, "_", Sys.Date(), ".pdf")
  },
  content = function(file) {
    sig <- paste0("B:SAMP:subtype_", input$manual_sig)
    vst_data <- result_env$vst_data
    pca_coord <- result_env$pca_coord
    gene_ids <- result_env$gene_ids

    genes <- unique_signatures[unique_signatures[, 3] == sig, 2]
    found_genes <- which(genes %in% gene_ids)

    if (length(found_genes) <= (length(genes) * 2 / 3)) {
      pdf(file)
      plot.new()
      text(0.5, 0.5, paste("⚠️ Not enough signature genes for", input$manual_sig), cex = 1.5)
      dev.off()
      return()
    }

    metagenes <- apply(vst_data, 2, function(x) mean(x[na.omit(match(genes, rownames(vst_data)))]))
    pca_temp <- pca_coord
    pca_temp$metagene <- metagenes

    pdf(file, width = 10, height = 6)
    print(
      ggplot(pca_temp, aes(PC1, PC2, color = metagene, label = rownames(pca_temp))) +
        geom_point(size = 3) +
        geom_text_repel(
          color = "black", size = 4,
          box.padding = unit(1, "lines"),
          point.padding = unit(1.5, "lines"),
          force = 10, force_pull = 0.1,
          max.overlaps = Inf, segment.color = "grey50",
          segment.curvature = -0.1, segment.ncp = 3, nudge_y = 15,
          show.legend = FALSE
        ) +
        scale_color_scico(palette = "buda", midpoint = mean(pca_temp$metagene),
                          name = paste(input$manual_sig, "metagene")) +
        coord_fixed() +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(size = 20, face = "bold"),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.margin = grid::unit(c(30, 10, 30, 10), "pt")
        ) +
        labs(
          title = paste("Feature plot of", input$manual_sig, "signature"),
          x = paste0("PC1 (", round(summary(prcomp(t(vst_data)))$importance[2, 1] * 100, 1), "%)"),
          y = paste0("PC2 (", round(summary(prcomp(t(vst_data)))$importance[2, 2] * 100, 1), "%)")
        ) +
        facet_wrap(vars(class))
    )
    dev.off()
  }
)

}
shinyApp(ui, server)


