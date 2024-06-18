# setwd("~/Documents/work/bio_analysis/bioinformatics_analysis")
source("scripts/format_data.R")
source("scripts/draw.R")
source("scripts/install.R")
gene_folder = "data/LPS/genes"
data_dir = "data/LPS"
out_dir = "out/LPS"

# heat map #####################################################################
## compare treatment KO - 5_g vs. PBS====
# var
pattern = "(.*KO.*5_g.*|.*KO.*PBS.*)"
feature = "TPM"

df_gene <- GeneFeatureExtractor$new(gene_folder, pattern, feature)$get_feature_hierarchy()

column_order <- df_gene$file_name_clean
column_level_2 <- df_gene$gender
column_level_1 <- df_gene$treatment

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_KO - 5_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = pattern)
  df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n, 

                                             pivot = T,
                                             feature = feature)

  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_KO - 5_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - KO - 5_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}


## compare treatment KO - 100_g vs. PBS====
# var
pattern = "(.*KO.*100_g.*|.*KO.*PBS.*)"
feature = "TPM"

df_gene <- GeneFeatureExtractor$new(gene_folder, pattern, feature)$get_feature_hierarchy()

column_order <- df_gene$file_name_clean
column_level_2 <- df_gene$gender
column_level_1 <- df_gene$treatment

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_KO - 100_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = pattern)
  df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n, 

                                             pivot = T,
                                             feature = feature)
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_KO - 100_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - KO - 100_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}

## compare treatment WT - 5_g vs. PBS====
# var
pattern = "(.*WT.*5_g.*|.*WT.*PBS.*)"
feature = "TPM"

df_gene <- GeneFeatureExtractor$new(gene_folder, pattern, feature)$get_feature_hierarchy()

column_order <- df_gene$file_name_clean
column_level_2 <- df_gene$gender
column_level_1 <- df_gene$treatment

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_WT - 5_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = pattern)
  df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n, 

                                             pivot = T,
                                             feature = feature)
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_WT - 5_g vs. PBS - top %d.png", top_n*2)), 

                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - WT - 5_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),

                                flag_row_name = flag_row_name
  )
}


## compare treatment WT - 100_g vs. PBS====
# var
pattern = "(.*WT.*100_g.*|.*WT.*PBS.*)"
feature = "TPM"

df_gene <- GeneFeatureExtractor$new(gene_folder, pattern, feature)$get_feature_hierarchy()
column_order <- df_gene$file_name_clean
column_level_2 <- df_gene$gender
column_level_1 <- df_gene$treatment

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_WT - 100_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = pattern)
  df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n, 
                                             pivot = T,
                                             feature = feature)
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_WT - 100_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - WT - 100_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}


## compare treatment GT - 5_g vs. PBS====
# var
pattern = "(.*GT.*5_g.*|.*GT.*PBS.*)"
feature = "TPM"

df_gene <- GeneFeatureExtractor$new(gene_folder, pattern, feature)$get_feature_hierarchy()
column_order <- df_gene$file_name_clean
column_level_2 <- df_gene$gender
column_level_1 <- df_gene$treatment

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_GT - 5_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = pattern)
  df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n, 
                                             pivot = T,
                                             feature = feature)
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_GT - 5_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - GT - 5_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}


## compare treatment GT - 100_g vs. PBS====
# var
pattern = "(.*GT.*100_g.*|.*GT.*PBS.*)"
feature = "TPM"

df_gene <- GeneFeatureExtractor$new(gene_folder, pattern, feature)$get_feature_hierarchy()
column_order <- df_gene$file_name_clean
column_level_2 <- df_gene$gender
column_level_1 <- df_gene$treatment

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_GT - 100_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = pattern)
  df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n, 
                                             pivot = T,
                                             feature = feature)
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_GT - 100_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - GT - 100_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}


## compare treatment SPURT - 5_g vs. PBS====
# var
pattern = "(.*SPURT.*5_g.*|.*SPURT.*PBS.*)"
feature = "TPM"

df_gene <- GeneFeatureExtractor$new(gene_folder, pattern, feature)$get_feature_hierarchy()
column_order <- df_gene$file_name_clean
column_level_2 <- df_gene$gender
column_level_1 <- df_gene$treatment

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_SPURT - 5_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = pattern)
  df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n, 
                                             pivot = T,
                                             feature = feature)
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_SPURT - 5_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - SPURT - 5_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}


## compare treatment SPURT - 100_g vs. PBS====
# var
pattern = "(.*SPURT.*100_g.*|.*SPURT.*PBS.*)"
feature = "TPM"

df_gene <- GeneFeatureExtractor$new(gene_folder, pattern, feature)$get_feature_hierarchy()
column_order <- df_gene$file_name_clean
column_level_2 <- df_gene$gender
column_level_1 <- df_gene$treatment

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_SPURT - 100_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = pattern)
  df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n, 
                                             pivot = T,
                                             feature = feature)
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_SPURT - 100_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - SPURT - 100_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}



# go analysis ##################################################################
## WT - 5_g vs. PBS =============
# var
log_2_FC = 1
top_n = 10

# process data
data = DifferentialExpressionDataProcessor$new(file_path = file.path(data_dir,"diff_expr_WT - 5_g vs. PBS.csv"))
genes = data$filter_by_log2_fc(log_2_FC)

# plot
go_plot <- GoAnalysisPlot$new(output = file.path(out_dir, "go_enrichment_WT - 5_g vs. PBS.png"), 
                              width = 2000, 
                              height = 4000, 
                              res = 300, 
                              bg = "white")
go_plot$draw(genes, top_n, "Go Analysis: WT - 5_g vs. PBS")

## KO - 5_g vs. PBS =============
# var
log_2_FC = 1
top_n = 10

# process data
data = DifferentialExpressionDataProcessor$new(file_path = file.path(data_dir,"diff_expr_KO - 5_g vs. PBS.csv"))
genes = data$filter_by_log2_fc(log_2_FC)

# plot
go_plot <- GoAnalysisPlot$new(output = file.path(out_dir, "go_enrichment_KO - 5_g vs. PBS.png"), 
                              width = 2000, 
                              height = 4000, 
                              res = 300, 
                              bg = "white")
go_plot$draw(genes, top_n, "Go Analysis: KO - 5_g vs. PBS")
## WT - 100_g vs. PBS =============
# var
log_2_FC = 1
top_n = 10

# process data
data = DifferentialExpressionDataProcessor$new(file_path = file.path(data_dir,"diff_expr_WT - 100_g vs. PBS.csv"))
genes = data$filter_by_log2_fc(log_2_FC)

# plot
go_plot <- GoAnalysisPlot$new(output = file.path(out_dir, "go_enrichment_WT - 100_g vs. PBS.png"), 
                              width = 2000, 
                              height = 4000, 
                              res = 300, 
                              bg = "white")
go_plot$draw(genes, top_n, "Go Analysis: WT - 100_g vs. PBS")

## KO - 100_g vs. PBS =============
# var
log_2_FC = 1
top_n = 10

# process data
data = DifferentialExpressionDataProcessor$new(file_path = file.path(data_dir,"diff_expr_KO - 100_g vs. PBS.csv"))
genes = data$filter_by_log2_fc(log_2_FC)

# plot
go_plot <- GoAnalysisPlot$new(output = file.path(out_dir, "go_enrichment_KO - 100_g vs. PBS.png"), 
                              width = 2000, 
                              height = 4000, 
                              res = 300, 
                              bg = "white")
go_plot$draw(genes, top_n, "Go Analysis: KO - 100_g vs. PBS")

# box plot #####################################################################
## KO - 5_g vs. PBS ===============================
# Var
top_n = 6
feature = "TPM"
pattern = "(.*KO.*5_g.*|.*KO.*PBS.*)"

# process data
data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_KO - 5_g vs. PBS.csv"), 
                          gene_folder = gene_folder, 
                          pattern = pattern
                          )
df_diff_expr_unpivot <- data$extract_top_n_gene_feature(top_n = top_n, 

                                                   pivot = F,
                                                   feature = feature)

# plot
box_plot <- BoxPlot$new(file.path(out_dir, "box_plot_KO - 5_g vs. PBS.png"), 
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw(df = df_diff_expr_unpivot, 
                       title = "Box Plot KO - 5_g vs. PBS")

## WT - 5_g vs. PBS ===============================
# Var
top_n = 6
feature = "TPM"
pattern = "(.*WT.*5_g.*|.*WT.*PBS.*)"

# process data
data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_WT - 5_g vs. PBS.csv"), 
                          gene_folder = gene_folder, 
                          pattern = pattern
)
df_diff_expr_unpivot <- data$extract_top_n_gene_feature(top_n = top_n, 

                                                   pivot = F,
                                                   feature = feature)

# plot
box_plot <- BoxPlot$new(file.path(out_dir, "box_plot_WT - 5_g vs. PBS.png"), 
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw(df = df_diff_expr_unpivot, 
                       title = "Box Plot WT - 5_g vs. PBS")

## KO - 100_g vs. PBS ===============================
# Var
top_n = 6
feature = "TPM"
pattern = "(.*KO.*100_g.*|.*KO.*PBS.*)"

# process data
data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_KO - 100_g vs. PBS.csv"), 
                          gene_folder = gene_folder, 
                          pattern = pattern
)
df_diff_expr_unpivot <- data$extract_top_n_gene_feature(top_n = top_n, 

                                                   pivot = F,
                                                   feature = feature)

# plot
box_plot <- BoxPlot$new(file.path(out_dir, "box_plot_KO - 100_g vs. PBS.png"), 
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw(df = df_diff_expr_unpivot, 
                       title = "Box Plot KO - 100_g vs. PBS")

## WT - 100_g vs. PBS ===============================
# Var
top_n = 6
feature = "TPM"
pattern = "(.*WT.*100_g.*|.*WT.*PBS.*)"

# process data
data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_WT - 100_g vs. PBS.csv"), 
                          gene_folder = gene_folder, 
                          pattern = pattern
)
df_diff_expr_unpivot <- data$extract_top_n_gene_feature(top_n = top_n, 

                                                   pivot = F,
                                                   feature = feature)

# plot
box_plot <- BoxPlot$new(file.path(out_dir, "box_plot_WT - 100_g vs. PBS.png"), 
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw(df = df_diff_expr_unpivot, 

                       title = "Box Plot WT - 100_g vs. PBS")