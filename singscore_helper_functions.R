
# define singscore wrapper function
singscore_function <- function(ranks_data, metadata, genesets, analysis_name, output_folder, input_type, tool_colour){

  annot_colouring <- data.frame(row.names = colnames(ranks_data)) %>%
      mutate(annotation = metadata$annotation[match(colnames(ranks_data), metadata[, c("filename", "annotation")]$filename)]) %>% 
      dplyr::arrange(annotation)
  
  ranks_data_ordered <- ranks_data %>% 
      dplyr::select(rownames(annot_colouring))
  
  invisible(lapply(c('ES_matrix', 'ES_plots'), function(x) dir.create(file.path(output_folder, x), recursive = TRUE)))
  
  # singscore requires the input to be formatted as output of the following function
  rankData <- rankGenes(ranks_data_ordered)
  singscore_res <- multiScore(rankData, upSet = genesets, knownDirection = T)
  singscore_df <- as.data.frame(singscore_res$Scores)
  singscore_score <- singscore_df[order(rownames(singscore_df)), , drop = FALSE]
  write.table(rownames_to_column(singscore_score, var = 'GOI_Set'),
              file.path(output_folder, paste0(analysis_name, "-fullNES.tsv")),
              sep = '\t', quote = FALSE, row.names = FALSE)

  plot_GSEApheatmap_wNAs(singscore_score,
      file.path(output_folder, 'ES_plots', paste0(analysis_name, '.png')),
      paste0(analysis_name, " - singscore_", input_type),
      tool_colour, wannotation = annot_colouring)
  # with legend
  plot_GSEApheatmap_wNAs(singscore_score,
      file.path(output_folder, 'ES_plots', paste0(analysis_name, '_annot.png')),
      paste0(analysis_name, " - singscore_", input_type),
      tool_colour, wannotation = annot_colouring, wlegend = TRUE)

  # run for single annotation groups for direct comparison of scores and plotting
  annotation_list <- unique(metadata$annotation)
  for (annotation_group in annotation_list){

      if (annotation_group == "" || is.na(annotation_group)){
          next
      }

      annotation_group_og <- annotation_group
      # simpler name saving (avoids file saving issue with spaces)
      if (grepl("\\(", annotation_group) & grepl("\\)", annotation_group)) {
          annotation_group <- sub(".*\\((.*)\\).*", "\\1", annotation_group)
      }
      annotation_group <- gsub(' ', '_', gsub('/', '-', gsub(' / ', '-', annotation_group)))
      print(annotation_group)

      metadata_of_annotation <- metadata %>% 
          filter(annotation == annotation_group_og)
      singscore_score_annot <- singscore_score[, colnames(singscore_score) %in% metadata_of_annotation$filename, drop = FALSE]
      
      write.table(rownames_to_column(singscore_score_annot, var = 'GOI_Set'),
                  file.path(output_folder, 'ES_matrix', paste0(annotation_group, ".tsv")),
                  sep = '\t', quote = FALSE, row.names = FALSE)
      
      plot_GSEApheatmap_wNAs(singscore_score_annot,
          file.path(output_folder, 'ES_plots', paste0(annotation_group, '.png')),
          paste0(annotation_group_og, " - singscore_", input_type),
          tool_colour)
  }
}


plot_GSEApheatmap_wNAs <- function(ES_matrix, png_name, plot_title, tool_colour, wannotation = NA, wlegend = FALSE){
    # function to plot all GSEA pheatmaps the same way, handling NAs in ES matrix outputs
    ES_matrix[ES_matrix == "---"] <- NA
    for (i in seq_along(ES_matrix)) {
        ES_matrix[[i]] <- suppressWarnings(as.numeric(ES_matrix[[i]]))
    }
    
    ES_matrix <- as.matrix(ES_matrix)
    
    scale_columns <- function(mat) {
        out <- mat  # keep original values unless scaling is appropriate
        
        for (j in seq_len(ncol(mat))) {
            col <- mat[, j]
            
            # if all values are NA → leave column untouched
            if (all(is.na(col))) next
            
            m <- mean(col, na.rm = TRUE)
            s <- sd(col, na.rm = TRUE)
            
            # If no variation: leave column as-is
            if (is.na(s) || s == 0) {
                next
            } else {
                out[, j] <- (col - m) / s
            }
        }
        out
    }
    
    scaled_ES <- scale_columns(ES_matrix)

    # --- Create guaranteed-unique breaks ---
    data_min <- min(scaled_ES, na.rm = TRUE)
    data_max <- max(scaled_ES, na.rm = TRUE)

    if (data_min == data_max) {
        data_min <- data_min - 0.1
        data_max <- data_max + 0.1
    }

    myColor <- colorRampPalette(c("white", "white", tool_colour))(101)
    breaks <- seq(data_min, data_max, length.out = length(myColor) + 1)

    png(png_name)
    pheatmap::pheatmap(scaled_ES,
                    show_rownames = T, show_colnames = F,
                    treeheight_row = 0, treeheight_col = 0,
                    cluster_cols = F, cluster_rows = F, scale = 'none',
                    color = myColor, breaks = breaks,
                    main = plot_title, fontsize_row = 5,
                    annotation_col = wannotation,
                    annotation_legend = wlegend)
    dev.off()
}

