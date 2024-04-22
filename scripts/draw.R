library(gplots) 
library(tidyr)
library(dplyr)
library(ComplexHeatmap)

# Draw Heatmap #################################################################
draw_heatmap <- function(input_file, column_order, column_level_1, column_level_2, column_title){
  data <- read.csv(input_file, check.names = FALSE)
  
  data_sorted <- data %>%
    arrange(Pathway, Molecule) 
  

  # Select and reorder columns
  data_matrix <- data_sorted %>%
    select(all_of(column_order)) %>%
    as.matrix()
  
  # Set row names to 'Pathway_Molecule'
  rownames(data_matrix) <- data_sorted$Molecule
  
  # Create a HeatmapAnnotation object for column annotations
  column_annotation <- HeatmapAnnotation(
    Treatment = column_level_1,
    Gender = column_level_2,
    annotation_name_side = "left"
  )
  
  row_annotation <- rowAnnotation(
    Pathway = data_sorted$Pathway, 
    annotation_name_side = "bottom",  
    show_annotation_name = TRUE 
  )
  
  # Normalize the data to Z-scores
  z_score_matrix <- t(scale(t(data_matrix)))
  
  
  # Create the heatmap with annotations
  Heatmap(
    z_score_matrix,
    name = "Expression",
    col = colorRampPalette(c("blue", "white", "red"))(75),
    cluster_rows = TRUE,  
    cluster_columns = TRUE,  
    # show_column_dend = FALSE,
    top_annotation = column_annotation,
    left_annotation = row_annotation,
    column_title = column_title,
    row_title = "Pathway - Molecule",
    row_title_side = "left",
    row_names_gp = gpar(fontsize = 4),
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Expression Level"),
    column_names_rot = 90,  
    row_names_side = "left",
    column_names_side = "bottom"
  )
}

# Draw Bar Chart ###############################################################
draw_bar <- function(input_file, pathways_list, title){
  df <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)
  df <- df %>%
    filter(`Ingenuity Canonical Pathways` %in% pathways_list) %>%
    arrange(desc(`-log(p-value)`))
  
  df$`z-score` <- as.numeric(as.character(df$`z-score`)) 
  
  max_z_score <-  max(abs(df$`z-score`), na.rm = TRUE)
  min_z_score = -max_z_score
  
  ggplot(df, aes(x = `-log(p-value)`,y = reorder(`Ingenuity Canonical Pathways`, `-log(p-value)`, decreasing = FALSE), fill = `z-score`)) + 
    geom_bar(stat = "identity") +  # Bar plot
    scale_fill_gradient2(low ='blue', high = "red", mid = "white",limits = c(min_z_score, max_z_score)) + 
    geom_text(aes(label = round(`-log(p-value)`, 2)), hjust = -0.2) +
    labs(
      title = title,
      x = "-log(p-value)",
      y = "Canonical Pathways",
      fill = "z-score"
    ) +
    theme_minimal() 
}
