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
	library(ggrepel)
	library(plotly)
	library(DESeq2)
	library(randomForest)
	library(DT)
	library(bslib)

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
    tags$style(HTML(".footer {display: flex; justify-content: space-between; align-items: center; margin-top: 40px; padding: 20px; font-weight: bold; font-size: 14px;} .about-box {margin-top: 10px; margin-bottom: 30px; padding: 20px; border: 1px solid #ddd; border-radius: 12px; background-color: #fdfdfd; box-shadow: 0 2px 6px rgba(0,0,0,0.05);} .dataTables_wrapper td {white-space: normal !important; word-break: break-word;} .card {margin-top: 20px;} .card-title {font-weight: bold; font-size: 18px;} .section-title {margin-top: 30px; font-size: 26px; font-weight: bold; border-bottom: 2px solid #eee; padding-bottom: 5px;}"))
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
               HTML("<h4><i class='fa-solid fa-dna'></i> About TASC:</h4> TASC (T-ALL Subtype Classifier) is a user-friendly R Shiny application designed to classify T-cell acute lymphoblastic leukemia samples into molecular subtypes using RNA-seq count data. It leverages a Random Forest model trained on curated expression signatures to provide subtype predictions and feature visualizations.")
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
                                    choices = c("Select..." = "", "Salmon" = "salmon", "featureCounts" = "featureCounts"),
                                    width = "100%")
                 ),
                 column(6,
                        fileInput("countsfile", "Upload count matrix (.tsv):", accept = c(".tsv", ".txt"), width = "100%")
                 )
               ),
               tags$div(id = "info_source_help", style = "font-size: 13px; background-color: #eef7fb; padding: 8px; border-radius: 5px; margin-top: 10px;",
                        HTML("<b>Info:</b> Please choose a valid count matrix source and upload a file with gene IDs (starting with <code>ENSG</code>) in the first column.")),
               textOutput("gene_id_error")
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
                      hr(),
                      uiOutput("Imputation_alert"),
                      downloadButton("downloadGenes", "Download Imputed Gene List")
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
  observeEvent(result(), {
    shinyjs::hide("waiting_message")
    shinyjs::show("results_section")
  })
  
  processedCounts <- reactive({
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
                             counts[, -1])
    list(counts = numeric_counts, gene_ids = gene_ids)
  })
  
  result <- reactive({
    pc <- processedCounts()
    req(pc)
    counts <- pc$counts
    gene_ids <- pc$gene_ids
    load("data/random_forest_model.rda")
    load("data/unique_signatures.rda")
    
    DGE_OBJ <- DESeqDataSetFromMatrix(round(counts), design = ~0,
                                      rowData = data.frame(gene_id = gene_ids),
                                      colData = data.frame(samples = colnames(counts)))
    vst_data <- assay(vst(DGE_OBJ, blind = TRUE))
    rownames(vst_data) <- sub("\\..*", "", gene_ids)
    colnames(vst_data) <- DGE_OBJ$samples
    
    gen <- sub("\\..*", "", model$subset_genes)
    missing_genes <- gen[!gen %in% rownames(vst_data)]
    
    if (length(missing_genes) > 0) {
      idx <- which(!gen %in% rownames(vst_data))
      sd_imput <- apply(model$original_data_matrix[model$subset_genes[idx], ], 1, sd)
      mean_imput <- apply(model$original_data_matrix[model$subset_genes[idx], ], 1, mean)
      sub <- t(sapply(seq_along(mean_imput), function(i) rnorm(ncol(vst_data), mean_imput[i], sd_imput[i])))
      rownames(sub) <- missing_genes
      colnames(sub) <- colnames(vst_data)
      vst_sub <- rbind(vst_data, sub)[gen, ]
      rownames(vst_sub) <- model$subset_genes
    } else {
      vst_sub <- vst_data[gen, ]
      rownames(vst_sub) <- model$subset_genes
    }

    
    class_tmp <- predict(model$rf, newdata = t(vst_sub))
    class_p <- predict(model$rf, newdata = t(vst_sub),type="prob")
    class <- class_tmp
    class[which(class_p[,grep("LMO2",colnames(class_p))]>model$threshold)]<-"B:SAMP:subtype_LMO2 γδ-like"
    cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                   "#0072B2", "#D55E00", "#CC79A7", "#999933", "#CC6666",
                   "#339933", "#9966CC", "#66CCEE")
    names(cbPalette) <- gsub("B:SAMP:subtype_", "", levels(class))
    
    pca_coord <- as.data.frame(prcomp(t(vst_data))$x)
    pca_coord$class <- gsub("B:SAMP:subtype_", "", as.vector(class))
    
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
        scale_color_gradient2(low = "purple", mid = "black", high = "yellow",
                              midpoint = mean(pca_temp$metagene),
                              name = paste(sig, "metagene")) +
        coord_fixed() +
	theme_minimal(base_size=14)+
        theme(plot.title = element_text(size = 20, face = "bold"),panel.border = element_rect(color = "black", fill = NA, size = 1)) +
        labs(title = paste("Feature plot of", sig, "signature"),
             x = paste0("PC1 (", round(summary(prcomp(t(vst_data)))$importance[2, 1] * 100, 1), "%)"),
             y = paste0("PC2 (", round(summary(prcomp(t(vst_data)))$importance[2, 2] * 100, 1), "%)"))+facet_wrap(vars(class))+
	theme(plot.margin = grid::unit(c(30, 10, 30, 10), "pt"))

    }
    
    structure(list(
      class = class_corrected,
      p = p,
      imp = missing_genes,
      plot_data = p,
      meta_plots = meta_plots
    ), model_genes = model$subset_genes)
  })
  
  output$classTable <- renderDT({
    req(result())
    df <- data.frame(Sample = names(result()$class), Prediction = result()$class)
    datatable(df, extensions = 'Buttons',
              options = list(dom = 'Blfrtip', buttons = list('copy', 'csv', 'excel', 'pdf', 'print'), pageLength = 15),
              rownames = FALSE)
  })
  
  output$pcaPlot <- renderPlotly({ ggplotly(result()$p) })
  
  output$Imputation_alert <- renderUI({
    res <- result()
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
        }
      )
    }
  })
  
  output$imputedGenes <- renderPrint({ result()$imp })
  
  output$metaPlots <- renderUI({
    lapply(names(result()$meta_plots), function(name) {
      column(width = 12, align = "center",
             plotOutput(outputId = paste0("meta_plot_", name), height = "800px"))
    }) |> tagList()
  })
  
  observe({
    plots <- result()$meta_plots
    lapply(names(plots), function(name) {
      output[[paste0("meta_plot_", name)]] <- renderPlot({ plots[[name]] })
    })
  })
  
  output$downloadPred <- downloadHandler(
    filename = function() paste0("class_predictions_", Sys.Date(), ".tsv"),
    content = function(file) {
      write.table(data.frame(Sample = names(result()$class), Prediction = result()$class),
                  file, sep = "\t", row.names = FALSE, quote = FALSE)
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
      lapply(result()$meta_plots, print)
      dev.off()
    }
  )
}

shinyApp(ui, server)


