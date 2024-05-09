library(gplots)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(readxl)

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(gridExtra)


# Draw Heatmap #################################################################
draw_heatmap <- function(input_file, column_order, column_level_1, column_level_2, column_title, output){
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
  p <- Heatmap(
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
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    heatmap_legend_param = list(title = "Expression Level"),
    column_names_rot = 90,  
    row_names_side = "left",
    column_names_side = "bottom"
  )

  png(output, width = 6000, height = 4800, res = 300, bg = "white")

  print(p)

  dev.off()
}
# Draw Standard Heatmap #################################################################
draw_heatmap_std <- function(df, column_order, column_level_1, column_level_2, column_title, flag_row_name, output){

  
  # data_sorted <- df
    # arrange(Pathway, Molecule) 
  

  # Select and reorder columns
  data_matrix <- df %>%
    dplyr::select(all_of(column_order)) %>%
    as.matrix()
  
  # Set row names to 'Pathway_Molecule'
  rownames(data_matrix) <- df$Name
  
  # Create a HeatmapAnnotation object for column annotations
  column_annotation <- HeatmapAnnotation(
    Gene = column_level_1,
    Treatment = column_level_2,
    annotation_name_side = "right"
  )

  
  # Normalize the data to Z-scores
  z_score_matrix <- t(scale(t(data_matrix)))
  
  
  # Create the heatmap with annotations
  p <- Heatmap(
    z_score_matrix,
    name = "Expression",
    col = colorRampPalette(c("blue", "white", "red"))(75),
    cluster_rows = T,  
    cluster_columns = FALSE,  
    # show_column_dend = FALSE,
    show_row_names = flag_row_name,
    top_annotation = column_annotation,
    # left_annotation = row_annotation,
    column_title = column_title,
    row_title = "Gene",
    row_title_side = "left",
    row_names_gp = gpar(fontsize = 14),
    column_names_gp = gpar(fontsize = 14),
    heatmap_legend_param = list(title = "Expression Level"),
    column_names_rot = 60,  
    row_names_side = "left",
    column_names_side = "bottom"
  )

  png(output, width = 6000, height = 4800, res = 300, bg = "white")

  print(p)

  dev.off()
}

# Draw Bar Chart ###############################################################
draw_bar <- function(input_file, pathways_list, title, output){
  df <- readxl::read_excel(input_file, skip = 1, col_names = TRUE)

  df <- df %>%
    filter(`Ingenuity Canonical Pathways` %in% pathways_list) %>%
    arrange(desc(`-log(p-value)`))
  
  df$`z-score` <- as.numeric(as.character(df$`z-score`)) 
  
  max_z_score <-  max(abs(df$`z-score`), na.rm = TRUE)
  min_z_score = -max_z_score
  
  p <- ggplot(df, aes(x = `-log(p-value)`,y = reorder(`Ingenuity Canonical Pathways`, `-log(p-value)`, decreasing = FALSE), fill = `z-score`)) + 
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
  print(p)

  ggsave(filename = output, plot = p, width = 15, height =12 , dpi = 300, bg = "white")
}


# Draw Go analysis map #########################################################
draw_go_analysis <- function(genes, top_n, title, output){
  
  GO_results <- enrichGO(gene = genes, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont='ALL')

  plot <- dotplot(GO_results, showCategory = top_n, split = "ONTOLOGY", title= title) + facet_grid(ONTOLOGY~., scale="free")
  
  png(output, width = 2000, height = 4000, res = 300, bg = "white")
  print(plot)
  dev.off()
}

# Draw box plot ##########################

draw_box_plot <- function(df, title, output){ 
  ggplot(df, aes(x=Name, y=TPM, fill=Treatment)) + 
           geom_boxplot() +
    labs(title = title, x = "Gene", y = "TPM") +
    theme_minimal()
  
  ggsave(output, width = 10, height = 8, bg = "white")

}