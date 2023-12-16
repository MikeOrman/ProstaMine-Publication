library(shiny)
library(DT)
library(shinybusy)
library(shinyalert)
library(curatedPCaData)
library(ComplexHeatmap)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(shinyAce)
library(rio)
library(S4Vectors)
library("edgeR")
library("dplyr")
library("VennDiagram")
library(ggplot2)
library(ggrepel)
library(DT)
library(metap)
library(ggpubr)
library(gridExtra)
library(survival)
library(survminer)
library(karyoploteR)

source("subsetting functions.R")
source("manhattan.R")
source("oncoprint.R")

# taylor.mae <- getPCa("taylor")
# taylor.CNA <- taylor.mae[["cna.gistic"]]
# barbieri.mae <- getPCa("barbieri")
# barbieri.CNA <- barbieri.mae[["cna.gistic"]]
# baca.mae <- getPCa("baca")
# baca.CNA <- baca.mae[["cna.gistic"]]
# hieronymus.mae <- getPCa("hieronymus")
# hieronymus.CNA <- hieronymus.mae[["cna.gistic"]]
# TCGA.mae <- getPCa("TCGA")
# TCGA.CNA <- TCGA.mae[["cna.gistic"]]
# abida.mae <- getPCa("abida")
# abida.CNA <- abida.mae[["cna.gistic"]]

Subtypes <- read.table("alteration subtypes.txt")
colnames(Subtypes) <- c("Subtype", "Alteration")

shinyServer <- function(input, output, session) {
  # Script for timeout feature
  options(shiny.maxRequestSize = 100 * 1024^2)
  #----Reactive Vals----
  vals <- reactiveValues(
    ws = NULL,
    gene_3=NULL,
    a1=NULL,
    a2=NULL,
    a3=NULL,
    a4=NULL,
    a5 = NULL,
    a6 = NULL,
    kp=NULL,
    conc_dge=NULL,
    genomic_plot=NULL,
    gleason_barplot=NULL)
  #----Tabs----
  # Select Background tab
  output$Subtype_table <- DT::renderDataTable({
    return(DT::datatable(Subtypes, selection = 'single', rownames = FALSE, options = list(scrollX = T)))
    shinyjs::show(id = "Subtypes")
    js$enableTabs();})
  output$Subtype <- renderText({
    paste(Subtypes[input$Subtype_table_rows_selected,1], 
          Subtypes[input$Subtype_table_rows_selected,2], sep = " ")
  })
  # Hits tab
  observe({
    stem <- paste(Subtypes[input$Subtype_table_rows_selected,1], 
                  Subtypes[input$Subtype_table_rows_selected,2], sep = "")
    filename <- paste(stem, ".RData", sep = "")
    load(filename)
    # Hits table
    index <- ((hits$`Primary Alteration Ratio Difference` >= input$cutoff1) &
                (hits$`Primary Enrichment FDR` <= input$cutoff2) &
                (hits$`Metastatic Alteration Ratio Difference` >= input$cutoff3) &
                (hits$`Metastatic Enrichment FDR` <= input$cutoff4) &
                (hits$`Primary DGE FDR` <= input$cutoff5) &
                (hits$`Metastatic DGE FDR` <= input$cutoff6) &
                (hits$`TCGA ST PFS pval` <= input$cutoff7))
    hits <- hits[index,]
    hits[hits$Alteration=="LOF",4] <- "Loss"
    hits[hits$Alteration=="GOF",4] <- "Gain"
    hits$`ProstaMine Score` <- signif(hits$`ProstaMine Score`, 2)
    hits[,5:20] <- signif(hits[,5:20], 2)
    colnames(hits)[c(5,7)] <- c("Primary Co-Alteration Frequency Difference", 
                                "Metastatic Co-Alteration Frequency Difference")
    output$hits_table <- DT::renderDataTable({
      return(DT::datatable(hits, selection = 'single', options = list(scrollX = T)))
      shinyjs::show(id = "hits_table")
      js$enableTabs();})
    # Rank and Karyoplot
    output$rank <- renderPlot({
      scores <- hits
      scores <- scores[,c(1,2)]
      scores$Gene <- factor(scores$Gene, levels = scores$Gene)
      scores$Rank <- 1:nrow(scores)
      gene_3=hits[[1]][input$hits_table_rows_selected]
      index <- match(gene_3, scores$Gene)
      vals$rank <- (ggplot() + geom_point(scores, mapping = aes(`Rank`, `ProstaMine Score`), size = 2)
                    + geom_point(scores[index,], mapping = aes(`Rank`, `ProstaMine Score`, color = "red"), size = 2)
                    + geom_text_repel(scores[index,], mapping = aes(`Rank`, `ProstaMine Score`, label = `Gene`),
                                      nudge_x = 5, nudge_y = 0.02, point.padding = 0.4, size = 6)
                    + xlab("ProstaMine Rank")
                    + ylim(0,1) + theme_classic()
                    + theme(axis.title.y.left  = element_text(size=26),
                            axis.text.y.left =  element_text(size=24),
                            axis.title.x.bottom = element_text(size=26),
                            axis.text.x.bottom = element_blank(),
                            axis.line = element_line(size = .1), 
                            axis.ticks = element_line(size = 0.2),
                            legend.position = "none"))
      return(vals$rank)
    })
    output$karyoplot <- renderPlot({
      a<-hits
      a <- a[,c(1,3)]
      a[,2] <- 0.2
      mp <- manhattan.object(a)
      mp <- sortSeqlevels(mp)
      gene_3=hits[[1]][input$hits_table_rows_selected]
      index <- match(gene_3, mp@ranges@NAMES)
      vals$kp <- plotKaryotype(plot.type = 2, chromosomes = levels(mp@seqnames), cex = 1.7)
      kpPoints(vals$kp, mp, col = "Black", cex = 1.3, ymin = 0, ymax = .8)
      kpPoints(vals$kp, mp[index], col = "Red", cex = 1.3, ymin = 0, ymax = .8)
      kpPlotMarkers(vals$kp, chr = as.character(mp@seqnames)[index], x = mp@ranges@start[index], 
                    labels = names(mp)[index], text.orientation = "horizontal", marker.parts = c(0,0,0))
      return(vals$kp)
    })
    # Oncoprint tab
    # output$plot1 <- renderPlot({
    #   gene_3=hits[[1]][input$hits_table_rows_selected]
    #   {vals$a1 <- oncoprint1(query.gene, gene_3, taylor.CNA)}
    #   return(vals$a1)
    #   js$enableTabs();})
    # output$plot2 <- renderPlot({
    #   gene_3=hits[[1]][input$hits_table_rows_selected]
    #   vals$a2<-oncoprint1(query.gene, gene_3, barbieri.CNA)
    #   return(vals$a2)
    #   js$enableTabs();})
    # output$plot3 <- renderPlot({
    #   gene_3=hits[[1]][input$hits_table_rows_selected]
    #   vals$a3 <-oncoprint1(query.gene, gene_3, baca.CNA)
    #   return(vals$a3)
    #   js$enableTabs();})
    # output$plot4 <- renderPlot({
    #   gene_3=hits[[1]][input$hits_table_rows_selected]
    #   vals$a4 <-oncoprint1(query.gene, gene_3, hieronymus.CNA)
    #   return(vals$a4)
    #   js$enableTabs();})
    # output$plot5 <- renderPlot({
    #   gene_3=hits[[1]][input$hits_table_rows_selected]
    #   vals$a5 <-oncoprint1(query.gene, gene_3, TCGA.CNA)
    #   return(vals$a5)
    #   js$enableTabs();})
    # output$plot6 <- renderPlot({
    #   gene_3=hits[[1]][input$hits_table_rows_selected]
    #   vals$a6 <-oncoprint1(query.gene, gene_3, abida.CNA)
    #   return(vals$a6)
    #   js$enableTabs();})
    #----Integrated Analysis tab----
    output$genomic_plot<-renderPlot({
      gene_3 = hits[[1]][input$hits_table_rows_selected]
      type = hits[[4]][input$hits_table_rows_selected]
      if (type == "Loss"){
        index <- match(gene_3, LOF.coalteration.enrichment$`Hugo Symbol`)
        primary.WT.altered <- LOF.coalteration.enrichment$`Set 2 altered count`[index]
        primary.WT.unaltered <- LOF.coalteration.enrichment$`Set 2 unaltered count`[index]
        primary.WT.rate <- primary.WT.altered / (primary.WT.altered + primary.WT.unaltered)
        primary.ST.altered <- LOF.coalteration.enrichment$`Set 1 altered count`[index]
        primary.ST.unaltered <- LOF.coalteration.enrichment$`Set 1 unaltered count`[index]
        primary.ST.rate <- primary.ST.altered / (primary.ST.altered + primary.ST.unaltered)
        met.ST.altered <- LOF.met.enrichment$`Set 1 altered count`[index]
        met.ST.unaltered <- LOF.met.enrichment$`Set 1 unaltered count`[index]
        met.ST.rate <- met.ST.altered / (met.ST.altered + met.ST.unaltered)
        met.WT.altered <- LOF.met.enrichment$`Set 2 altered count`[index]
        met.WT.unaltered <- LOF.met.enrichment$`Set 2 unaltered count`[index]
        met.WT.rate <- met.WT.altered / (met.WT.altered + met.WT.unaltered)
        x <- factor(c("Primary\nTumors", "Primary\nTumors", "Metastatic\nTumors", "Metastatic\nTumors"), 
                    levels = c("Primary\nTumors", "Metastatic\nTumors"))
        y <- c(primary.WT.rate, primary.ST.rate, met.WT.rate, met.ST.rate)
        st <- factor(c(query.gene.wt, query.gene.alteration, query.gene.wt, query.gene.alteration),
                     levels = c(query.gene.wt, query.gene.alteration))
        genomic.plot.data <- data.frame(x, y, st)
        index <-       index <- match(gene_3, hits$Gene)
        stats.met <- t(c("Alteration Frequency", "Primary\nTumors", "Metastatic\nTumors",
                         paste("FDR=", formatC(hits$`Metastatic Enrichment FDR`[index], format = "e", digits = 2), sep = "")))
        stats.met <- as.data.frame(stats.met)
        colnames(stats.met) <- c(".y.", "group1", "group2", "p")
        stats.primary <- t(c("Alteration Frequency", "Primary\nTumors", "Primary\nTumors",
                             paste("FDR=", formatC(hits$`Primary Enrichment FDR`[index], format = "e", digits = 2), sep = "")))
        stats.primary <- as.data.frame(stats.primary)
        colnames(stats.primary) <- c(".y.", "group1", "group2", "p")
        vals$genomic_plot <- (ggplot(genomic.plot.data, aes(x, y))
                              + geom_bar(aes(fill = st), stat='identity', position = "dodge", color = NA)
                              + scale_fill_manual(values = c("black", "blue"))
                              + ggtitle(paste(gene_3, "Loss Alterations", sep = " "))
                              + ylab("Alteration Frequency") + ylim(0,1)
                              + theme_classic()
                              + theme(legend.position = "bottom", legend.title = element_blank(),
                                      axis.title.x.bottom = element_blank())
                              + stat_pvalue_manual(stats.primary, label = "p", tip.length = 0.01, 
                                                   y.position = max(genomic.plot.data$y[genomic.plot.data$x == "Primary\nTumors"]) 
                                                   + .04, size = 6)
                              + stat_pvalue_manual(stats.met, label = "p", tip.length = 0.01,
                                                   y.position = max(genomic.plot.data$y[genomic.plot.data$x == "Metastatic\nTumors"])
                                                   + .04, size = 6)
                              + theme(axis.text.x = element_text(size=20),
                                      axis.text.y = element_text(size=20),
                                      axis.title.y = element_text(size = 20),
                                      axis.line = element_line(size = 1),
                                      legend.key.size = unit(1, "cm"),
                                      legend.text = element_text(size=20),
                                      title = element_text(size=20),
                                      legend.position = c(.45, .9)))
        return(vals$genomic_plot)
      }
      if (type == "Gain"){
        index <- match(gene_3, GOF.coalteration.enrichment$`Hugo Symbol`)
        primary.WT.altered <- GOF.coalteration.enrichment$`Set 2 altered count`[index]
        primary.WT.unaltered <- GOF.coalteration.enrichment$`Set 2 unaltered count`[index]
        primary.WT.rate <- primary.WT.altered / (primary.WT.altered + primary.WT.unaltered)
        primary.ST.altered <- GOF.coalteration.enrichment$`Set 1 altered count`[index]
        primary.ST.unaltered <- GOF.coalteration.enrichment$`Set 1 unaltered count`[index]
        primary.ST.rate <- primary.ST.altered / (primary.ST.altered + primary.ST.unaltered)
        met.ST.altered <- GOF.met.enrichment$`Set 1 altered count`[index]
        met.ST.unaltered <- GOF.met.enrichment$`Set 1 unaltered count`[index]
        met.ST.rate <- met.ST.altered / (met.ST.altered + met.ST.unaltered)
        met.WT.altered <- GOF.met.enrichment$`Set 2 altered count`[index]
        met.WT.unaltered <- GOF.met.enrichment$`Set 2 unaltered count`[index]
        met.WT.rate <- met.WT.altered / (met.WT.altered + met.WT.unaltered)
        x <- factor(c("Primary\nTumors", "Primary\nTumors", "Metastatic\nTumors", "Metastatic\nTumors"), 
                    levels = c("Primary\nTumors", "Metastatic\nTumors"))
        y <- c(primary.WT.rate, primary.ST.rate, met.WT.rate, met.ST.rate)
        st <- factor(c(query.gene.wt, query.gene.alteration, query.gene.wt, query.gene.alteration),
                     levels = c(query.gene.wt, query.gene.alteration))
        genomic.plot.data <- data.frame(x, y, st)
        index <-       index <- match(gene_3, hits$Gene)
        stats.met <- t(c("Alteration Frequency", "Primary\nTumors", "Metastatic\nTumors",
                         paste("FDR=", formatC(hits$`Metastatic Enrichment FDR`[index], format = "e", digits = 2), sep = "")))
        stats.met <- as.data.frame(stats.met)
        colnames(stats.met) <- c(".y.", "group1", "group2", "p")
        stats.primary <- t(c("Alteration Frequency", "Primary\nTumors", "Primary\nTumors",
                             paste("FDR=", formatC(hits$`Primary Enrichment FDR`[index], format = "e", digits = 2), sep = "")))
        stats.primary <- as.data.frame(stats.primary)
        colnames(stats.primary) <- c(".y.", "group1", "group2", "p")
        vals$genomic_plot <- (ggplot(genomic.plot.data, aes(x, y))
                              + geom_bar(aes(fill = st), stat='identity', position = "dodge", color = NA)
                              + scale_fill_manual(values = c("black", "blue"))
                              + ggtitle(paste(gene_3, "Gain Alterations", sep = " "))
                              + ylab("Alteration Frequency") + ylim(0,1)
                              + theme_classic()
                              + theme(legend.position = "bottom", legend.title = element_blank(),
                                      axis.title.x.bottom = element_blank())
                              + stat_pvalue_manual(stats.primary, label = "p", tip.length = 0.01, 
                                                   y.position = max(genomic.plot.data$y[genomic.plot.data$x == "Primary\nTumors"]) 
                                                   + .04, size = 6)
                              + stat_pvalue_manual(stats.met, label = "p", tip.length = 0.01,
                                                   y.position = max(genomic.plot.data$y[genomic.plot.data$x == "Metastatic\nTumors"])
                                                   + .04, size = 6)
                              + theme(axis.text.x = element_text(size=20),
                                      axis.text.y = element_text(size=20),
                                      axis.title.y = element_text(size = 20),
                                      axis.line = element_line(size = 1),
                                      legend.key.size = unit(1, "cm"),
                                      legend.text = element_text(size=20),
                                      title = element_text(size=20),
                                      legend.position = c(.45, .9)))
        return(vals$genomic_plot)
      }
    })
    output$gex_plot <- renderPlot({
      gene_3=hits[[1]][input$hits_table_rows_selected]
      Alteration <- hits[[4]][input$hits_table_rows_selected]
      DGE.boxplot <- data.frame()
      DGE.boxplot[1,1] <- gene_3
      if (Alteration == "Loss"){
        DGE.boxplot[1,2] <- "Loss"
        primary.altered <- alteration_subset(primary.ST, c(gene_3, "altered"))
        primary.unaltered <- alteration_subset(primary.ST, c(gene_3, "unaltered"))
        primary.altered.name <- paste(gene_3, "Loss", sep = "\n")
        primary.unaltered.name <- paste(gene_3, "WT", sep = "\n")
        met.altered <- alteration_subset(met.ST, c(gene_3, "altered"))
        met.unaltered <- alteration_subset(met.ST, c(gene_3, "unaltered"))
        met.altered.name <- paste(gene_3, "Loss", sep = "\n")
        met.unaltered.name <- paste(gene_3, "WT", sep = "\n")
      }
      if (Alteration == "Gain"){
        DGE.boxplot[1,2] <- "Gain"
        primary.altered <- alteration_subset(primary.ST, c(gene_3, "altered"))
        primary.unaltered <- alteration_subset(primary.ST, c(gene_3, "unaltered"))
        primary.altered.name <- paste(gene_3, "Gain", sep = "\n")
        primary.unaltered.name <- paste(gene_3, "WT", sep = "\n")
        met.altered <- alteration_subset(met.ST, c(gene_3, "altered"))
        met.unaltered <- alteration_subset(met.ST, c(gene_3, "unaltered"))
        met.altered.name <- paste(gene_3, "Gain", sep = "\n")
        met.unaltered.name <- paste(gene_3, "WT", sep = "\n")
      }
      #---------------------------BOXPLOTS OF CONCORDANT DGE--------------------------
      primary.unaltered.patients <- colnames(primary.unaltered)
      primary.altered.patients <- colnames(primary.altered)
      met.unaltered.patients <- colnames(met.unaltered)
      met.altered.patients <- colnames(met.altered)
      #-----Subset Primary mRNA data-----
      primary.altered.GEX.patients <- intersect(primary.altered.patients, colnames(primary.mRNA)[1:ncol(primary.mRNA)])
      primary.unaltered.GEX.patients <- intersect(primary.unaltered.patients, colnames(primary.mRNA)[1:ncol(primary.mRNA)])
      index <- c(TRUE)
      for (k in 1:ncol(primary.mRNA)){
        if (sum(colnames(primary.mRNA)[k] == primary.altered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      primary.altered.GEX <- primary.mRNA[, index]
      index <- c(TRUE)
      for (k in 1:ncol(primary.mRNA)){
        if (sum(colnames(primary.mRNA)[k] == primary.unaltered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      primary.unaltered.GEX <- primary.mRNA[, index]
      # Combine primary data into single data frame
      index <- match(gene_3, rownames(primary.altered.GEX))
      primary.altered.vec <- as.numeric(primary.altered.GEX[index,1:ncol(primary.altered.GEX)])
      primary.unaltered.vec <- as.numeric(primary.unaltered.GEX[index,1:ncol(primary.unaltered.GEX)])
      data.primary <- data.frame()
      data.primary[1:length(primary.unaltered.vec),1] <- primary.unaltered.name
      data.primary[1:length(primary.unaltered.vec),2] <- primary.unaltered.vec
      data.primary[((length(primary.unaltered.vec)+1):(length(primary.unaltered.vec)+length(primary.altered.vec))),1] <- primary.altered.name
      data.primary[((length(primary.unaltered.vec)+1):(length(primary.unaltered.vec)+length(primary.altered.vec))),2] <- primary.altered.vec
      data.primary[,3] <- paste("Primary Tumors\n", query.gene.alteration)
      colnames(data.primary) <- c("Alteration Status",
                                  paste(gene_3, "Normalized Gene Expression", sep = " "),
                                  "Patient Background")
      #-----Subset met mRNA data-----
      met.altered.GEX.patients <- intersect(met.altered.patients, colnames(met.mRNA)[1:ncol(met.mRNA)])
      met.unaltered.GEX.patients <- intersect(met.unaltered.patients, colnames(met.mRNA)[1:ncol(met.mRNA)])
      index <- c(TRUE)
      for (k in 1:ncol(met.mRNA)){
        if (sum(colnames(met.mRNA)[k] == met.altered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      met.altered.GEX <- met.mRNA[, index]
      index <- c(TRUE)
      for (k in 1:ncol(met.mRNA)){
        if (sum(colnames(met.mRNA)[k] == met.unaltered.GEX.patients) == 1){
          index[k] <- TRUE} else {index[k] <- FALSE}
      }
      met.unaltered.GEX <- met.mRNA[, index]
      # Combine met data into single data frame
      index <- match(gene_3, rownames(met.altered.GEX))
      met.altered.vec <- as.numeric(met.altered.GEX[index,1:ncol(met.altered.GEX)])
      met.unaltered.vec <- as.numeric(met.unaltered.GEX[index,1:ncol(met.unaltered.GEX)])
      data.met <- data.frame()
      data.met[1:length(met.unaltered.vec),1] <- met.unaltered.name
      data.met[1:length(met.unaltered.vec),2] <- met.unaltered.vec
      data.met[((length(met.unaltered.vec)+1):(length(met.unaltered.vec)+length(met.altered.vec))),1] <- met.altered.name
      data.met[((length(met.unaltered.vec)+1):(length(met.unaltered.vec)+length(met.altered.vec))),2] <- met.altered.vec
      data.met[,3] <- paste("Metastatic Tumors\n", query.gene.alteration)
      colnames(data.met) <- c("Alteration Status",
                              paste(gene_3, "Normalized Gene Expression", sep = " "),
                              "Patient Background")
      #-----Bind datasets. Define order of factor levels for plotting----
      data <- rbind(data.primary, data.met)
      if (Alteration == "Loss"){
        data$`Alteration Status` <- factor(data$`Alteration Status`, 
                                           levels = c(paste(gene_3, "WT", sep = "\n"), 
                                                      paste(gene_3, "Loss", sep = "\n")))
      }
      if (Alteration == "Gain"){
        data$`Alteration Status` <- factor(data$`Alteration Status`, 
                                           levels = c(paste(gene_3, "WT", sep = "\n"), 
                                                      paste(gene_3, "Gain", sep = "\n")))
      }
      data$`Patient Background` <- factor(data$`Patient Background`,
                                          levels = c(paste("Primary Tumors\n", query.gene.alteration), 
                                                     paste("Metastatic Tumors\n", query.gene.alteration)))
      #--------PLOT DGE-------------
      # Format statistics table for plotting
      facet <- c(paste("Primary Tumors\n", query.gene.alteration, sep = " "),
                 paste("Metastatic Tumors\n", query.gene.alteration, sep = " "))
      y <- paste(gene_3, "mRNA Expression z-score")
      if (Alteration == "Loss"){
        group1 <- paste(gene_3, "WT", sep = " ")
        group2 <- paste(gene_3, "Loss", sep = " ")
      }
      if (Alteration == "Gain"){
        group1 <- paste(gene_3, "WT", sep = " ")
        group2 <- paste(gene_3, "Gain", sep = " ")
      }
      FDR <- c(paste("FDR=", formatC(hits$`Primary DGE FDR`[hits$`Gene`==gene_3],
                                     format = "e", digits = 2), sep = ""),
               paste("FDR=", formatC(hits$`Metastatic DGE FDR`[hits$`Gene`==gene_3],
                                     format = "e", digits = 2), sep = ""))
      xmin <- c(1,1)
      xmax <- c(2,2)
      primary.y.position <- max(data[data$`Patient Background` == paste("Primary Tumors\n", query.gene.alteration, sep = " "),2]) + 0.2
      met.y.position <- max(data[data$`Patient Background` == paste("Metastatic Tumors\n", query.gene.alteration, sep = " "),2]) + 0.2
      y.positions <- c(primary.y.position, met.y.position)
      stats <- data.frame(facet, "Alteration Status", y, group1, group2, FDR, xmin, xmax, y.positions)
      colnames(stats) <- c("Patient Background", "Alteration Status", ".y.", "group1", "group2", "p", "xmin", "xmax", "y.position")
      #Make boxplot
      vals$gex_plot <- (ggplot(data, aes(`Alteration Status`, data[,2], color = `Alteration Status`))
                        + scale_color_manual(values=c("Black", "Blue"))
                        + geom_boxplot() + geom_jitter() 
                        + ylab(paste(gene_3, "mRNA Expression z-score"))
                        + facet_wrap(~`Patient Background`)
                        + theme_classic()
                        + theme(axis.title.x = element_blank(),
                                strip.text = element_text(size=20),
                                axis.title.y.left = element_text(size=20),
                                axis.text.x.bottom = element_text(size=20),
                                axis.text.y.left = element_text(size=20),
                                axis.line = element_line(size = 1),
                                strip.background = element_rect(size = 0.25),
                                legend.text = element_text(size=20),
                                legend.key.size = unit(1, "cm"),
                                legend.title = element_text(size=20))
                        + stat_pvalue_manual(stats, label = "p", size = 6))
      return(vals$gex_plot)
    })
    output$gleason_wt_plot <- renderPlot({
      gene_3=hits[[1]][input$hits_table_rows_selected]
      sample.names <- colnames(combined.WT)
      data <- data.frame()
      for (i in 1:length(sample.names)){
        index <- match(sample.names[i], combined.gleason$`Sample ID`)
        if (is.na(index) == FALSE & is.na(combined.gleason$`Grade Group`[index]) == FALSE
            & length(combined.gleason$`Grade Group`[index]) > 0){
          data[i,1] <- combined.gleason$`Grade Group`[index]
          if (data[i,1] == ">=8") {data[i,2] <- "High"}
          if (data[i,1] == "4+3" ) {data[i,2] <- "Intermediate"}
          if (data[i,1] == "3+4" ) {data[i,2] <- "Intermediate"}
          if (data[i,1] <= "<=6") {data[i,2] <- "Low"}
          gene.index <- match(gene_3, rownames(combined.mRNA))
          sample.index <- match(sample.names[i], colnames(combined.mRNA))
          data[i,3] <- combined.mRNA[gene.index, sample.index]
        }
      }
      data <- na.omit(data)
      colnames(data) <- c("Grade Group", "Risk Group", "Expression")
      data <- data[order(data[,3], decreasing = TRUE),]
      quantiles <- quantile(data$Expression, c(0.5))
      high.expressing <- data[data$Expression >= quantiles[1],]
      high.expressing[,4] <- "Upper 50%"
      colnames(high.expressing)[4] <- "mRNA Expression"
      low.expressing <- data[data$Expression <= quantiles[1],]
      low.expressing[,4] <- "Lower 50%"
      colnames(low.expressing)[4] <- "mRNA Expression"
      data <- rbind(low.expressing, high.expressing)
      data$`mRNA Expression` <- factor(data$`mRNA Expression`, levels = c("Lower 50%", "Upper 50%"))
      FDR.index <- match(gene_3, hits$`Gene`)
      FDR <- hits$`WT Gleason Enrichment FDR`[FDR.index]
      data$`Grade Group` <- factor(data$`Grade Group`, levels =c(">=8", "4+3", "3+4", "<=6"))
      stats <- t(c("Count", "Lower 50%", "Upper 50%",
                   paste("FDR=", formatC(FDR, format = "e", digits = 2), sep = "")))
      stats <- as.data.frame(stats)
      colnames(stats) <- c(".y.", "group1", "group2", "p")
      vals$gleason.wt <- print((ggplot(data, aes(`mRNA Expression`)) + geom_bar(aes(fill=`Grade Group`)) 
                                + labs(title = wt, x = paste(gene_3, "mRNA Expression"), y = "Count")
                                + ylim(0,(nrow(data)*.60)) + stat_pvalue_manual(stats, label = "p", y.position = nrow(data)*.55, size = 7)
                                + theme_classic()
                                + theme(axis.text.y.left = element_text(size=20), 
                                        axis.text.x.bottom = element_text(size=20),
                                        axis.line = element_line(size = 1),
                                        title = element_text(size=20), 
                                        legend.text = element_text(size=20))
                                + scale_fill_brewer(palette = "Set1", direction = 1)))
      return(vals$gleason.wt)
    })
    output$gleason_st_plot <- renderPlot({
      gene_3=hits[[1]][input$hits_table_rows_selected]
      sample.names <- colnames(combined.ST)
      data <- data.frame()
      for (i in 1:length(sample.names)){
        index <- match(sample.names[i], combined.gleason$`Sample ID`)
        if (is.na(index) == FALSE & is.na(combined.gleason$`Grade Group`[index]) == FALSE
            & length(combined.gleason$`Grade Group`[index]) > 0){
          data[i,1] <- combined.gleason$`Grade Group`[index]
          if (data[i,1] == ">=8") {data[i,2] <- "High"}
          if (data[i,1] == "4+3" ) {data[i,2] <- "Intermediate"}
          if (data[i,1] == "3+4" ) {data[i,2] <- "Intermediate"}
          if (data[i,1] <= "<=6") {data[i,2] <- "Low"}
          gene.index <- match(gene_3, rownames(combined.mRNA))
          sample.index <- match(sample.names[i], colnames(combined.mRNA))
          data[i,3] <- combined.mRNA[gene.index, sample.index]
        }
      }
      data <- na.omit(data)
      colnames(data) <- c("Grade Group", "Risk Group", "Expression")
      data <- data[order(data[,3], decreasing = TRUE),]
      quantiles <- quantile(data$Expression, c(0.5))
      high.expressing <- data[data$Expression >= quantiles[1],]
      high.expressing[,4] <- "Upper 50%"
      colnames(high.expressing)[4] <- "mRNA Expression"
      low.expressing <- data[data$Expression <= quantiles[1],]
      low.expressing[,4] <- "Lower 50%"
      colnames(low.expressing)[4] <- "mRNA Expression"
      data <- rbind(low.expressing, high.expressing)
      data$`mRNA Expression` <- factor(data$`mRNA Expression`, levels = c("Lower 50%", "Upper 50%"))
      FDR.index <- match(gene_3, hits$`Gene`)
      FDR <- hits$`ST Gleason Enrichment FDR`[FDR.index]
      data$`Grade Group` <- factor(data$`Grade Group`, levels =c(">=8", "4+3", "3+4", "<=6"))
      stats <- t(c("Count", "Lower 50%", "Upper 50%",
                   paste("FDR=", formatC(FDR, format = "e", digits = 2), sep = "")))
      stats <- as.data.frame(stats)
      colnames(stats) <- c(".y.", "group1", "group2", "p")
      vals$gleason.st <- print((ggplot(data, aes(`mRNA Expression`)) + geom_bar(aes(fill=`Grade Group`)) 
                                + labs(title = st, x = paste(gene_3, "mRNA Expression"), y = "Count")
                                + ylim(0,(nrow(data)*.60)) + stat_pvalue_manual(stats, label = "p", y.position = nrow(data)*.55, size = 7)
                                + theme_classic()
                                + theme(axis.text.y.left = element_text(size=20), 
                                        axis.text.x.bottom = element_text(size=20),
                                        axis.line = element_line(size = 1),
                                        title = element_text(size=20), 
                                        legend.text = element_text(size=20))
                                + scale_fill_brewer(palette = "Set1", direction = 1)))
      return(vals$gleason.st)
    })
    output$survival_wt_plot <- renderPlot({
      gene_3=hits[[1]][input$hits_table_rows_selected]
      data <- data.frame()
      sample.names.WT <- colnames(TCGA.WT)
      for (i in 1:length(sample.names.WT)){
        index <- match(sample.names.WT[i], TCGA.PFS$`Sample ID`)
        sample <- TCGA.PFS$`Patient ID`[index]
        index2 <- match(sample, TCGA.PFS$`Patient ID`)
        if (is.na(index) == FALSE & is.na(index2) == FALSE){
          data[i,1] <- TCGA.PFS$`Patient ID`[index2]
          data[i,2] <- (TCGA.PFS$Time)[index2]
          data[i,3] <- TCGA.PFS$Status[index2]
        }
      }
      for (i in 1:length(sample.names.WT)){
        index <- match(sample.names.WT[i], colnames(TCGA.mRNA))
        index2 <- match(gene_3, rownames(TCGA.mRNA))
        if (is.na(index) == FALSE & is.na(index2) == FALSE){
          data[i,4] <- TCGA.mRNA[index2,index]
        }
      }
      data <- na.omit(data)
      colnames(data) <- c("Patient ID", "DFS Months", "DFS Status", "Expression")
      for (i in 1:nrow(data)){
        if (data[i,3] == 0) {data[i,3] <- 0}
        if (data[i,3] == 1) {data[i,3] <- 1}
      }
      data <- data[order(data$Expression, decreasing = TRUE),]
      quantiles <- quantile(data$Expression, c(0.5))
      high.expressing <- data[data$Expression >= quantiles[1],]
      high.expressing[,5] <- paste("Upper 50%", paste(gene_3, "mRNA Expression", sep = " "), sep = " ")
      colnames(high.expressing)[5] <- "Group"
      low.expressing <- data[data$Expression <= quantiles[1],]
      low.expressing[,5] <- paste("Lower 50%", paste(gene_3, "mRNA Expression", sep = " "), sep = " ")
      colnames(low.expressing)[5] <- "Group"
      high.expressing <- high.expressing[,c(1, 5, 2, 3)]
      low.expressing <- low.expressing[,c(1, 5, 2, 3)]
      surv_data.WT <- rbind(high.expressing, low.expressing)
      surv_data.WT$`DFS Months` <- as.numeric(surv_data.WT$`DFS Months`)
      surv_data.WT$`DFS Status` <- as.numeric(surv_data.WT$`DFS Status`)
      surv_obj.WT <- Surv(time = surv_data.WT$`DFS Months`, event = surv_data.WT$`DFS Status`)
      fit.WT <- do.call(survfit, list(surv_obj.WT ~ Group, data = surv_data.WT))
      vals$survival.wt.plot <- ggsurvplot(fit.WT, data = surv_data.WT,
                                          pval = TRUE, log.rank.weights = "S1",
                                          risk.table = TRUE,
                                          palette = c("gold3", "blue2"),
                                          size = 1,
                                          risk.table.y.text = FALSE,
                                          ylab = "Ratio Tumors Progression-Free",
                                          xlab = "Months",
                                          legend.title = "",
                                          legend.labs = c(low.expressing[1,2], high.expressing[1,2]),
                                          xlim = c(0,80),
                                          title = wt,
                                          fontsize = 6,
                                          legend = c(0.45,0.15),
                                          pval.size = 7,
                                          pval.coord = c(0.45,0.4),
                                          ggtheme = theme_survminer(font.y = c(20, "plain", "black"),
                                                                    font.x =  c(20, "plain", "black"),
                                                                    font.tickslab =  c(20, "plain", "black"),
                                                                    font.legend =  c(20, "plain", "black"),
                                                                    font.main =  c(22, "plain", "black")))
      return(vals$survival.wt.plot)
    })
    output$survival_st_plot <- renderPlot({
      gene_3=hits[[1]][input$hits_table_rows_selected]
      data <- data.frame()
      sample.names <- colnames(TCGA.ST)
      for (i in 1:length(sample.names)){
        index <- match(sample.names[i], TCGA.PFS$`Sample ID`)
        sample <- TCGA.PFS$`Patient ID`[index]
        index2 <- match(sample, TCGA.PFS$`Patient ID`)
        if (is.na(index) == FALSE & is.na(index2) == FALSE){
          data[i,1] <- TCGA.PFS$`Patient ID`[index2]
          data[i,2] <- (TCGA.PFS$Time)[index2]
          data[i,3] <- TCGA.PFS$Status[index2]
        }
      }
      for (i in 1:length(sample.names)){
        index <- match(sample.names[i], colnames(TCGA.mRNA))
        index2 <- match(gene_3, rownames(TCGA.mRNA))
        if (is.na(index) == FALSE & is.na(index2) == FALSE){
          data[i,4] <- TCGA.mRNA[index2,index]
        }
      }
      data <- na.omit(data)
      colnames(data) <- c("Patient ID", "DFS Months", "DFS Status", "Expression")
      data <- data[order(data$Expression, decreasing = TRUE),]
      quantiles <- quantile(data$Expression, c(0.5))
      high.expressing <- data[data$Expression >= quantiles[1],]
      high.expressing[,5] <- paste("Upper 50%", paste(gene_3, "mRNA Expression", sep = " "), sep = " ")
      colnames(high.expressing)[5] <- "Group"
      low.expressing <- data[data$Expression <= quantiles[1],]
      low.expressing[,5] <- paste("Lower 50%", paste(gene_3, "mRNA Expression", sep = " "), sep = " ")
      colnames(low.expressing)[5] <- "Group"
      high.expressing <- high.expressing[,c(1, 5, 2, 3)]
      low.expressing <- low.expressing[,c(1, 5, 2, 3)]
      surv_data <- rbind(high.expressing, low.expressing)
      surv_data$`DFS Months` <- as.numeric(surv_data$`DFS Months`)
      surv_data$`DFS Status` <- as.numeric(surv_data$`DFS Status`)
      surv_object <- Surv(time = surv_data$`DFS Months`, event = surv_data$`DFS Status`)
      fit <- do.call(survfit, list(surv_object ~ Group, data = surv_data))
      vals$survival.st.plot <- ggsurvplot(fit, data = surv_data, 
                                          pval = TRUE, log.rank.weights = "S1",
                                          risk.table = TRUE,
                                          palette = c("gold3", "blue2"),
                                          size = 1,
                                          risk.table.y.text = FALSE,
                                          ylab = "Ratio Tumors Progression-Free",
                                          xlab = "Months",
                                          legend.title = "",
                                          legend.labs = c(low.expressing[1,2], high.expressing[1,2]),
                                          xlim = c(0,80),
                                          title = st,
                                          fontsize = 6,
                                          legend = c(0.45,0.15),
                                          pval.size = 7,
                                          pval.coord = c(0.45,0.4),
                                          ggtheme = theme_survminer(font.y = c(20, "plain", "black"),
                                                                    font.x =  c(20, "plain", "black"),
                                                                    font.tickslab =  c(20, "plain", "black"),
                                                                    font.legend =  c(20, "plain", "black"),
                                                                    font.main =  c(22, "plain", "black")))
      return(vals$survival.st.plot)
    })
  })
}