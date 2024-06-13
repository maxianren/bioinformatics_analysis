# setwd("~/Documents/work/bio_analysis/bioinformatics_analysis")
source("scripts/format_data.R")
source("scripts/draw.R")
gene_folder = "data/LPS/genes"
data_dir = "data/LPS"
out_dir = "out/LPS"

# heat map #####################################################################
## compare treatment KO - 5_g vs. PBS====
# var
column_level_1 = c( "PBS","PBS","LPS_5_ug","LPS_5_ug","LPS_5_ug","LPS_5_ug")
column_level_2 = c( "M","F","M","M","F","F")
column_order = c("_11_pKO-M-PBS_TPM", 
                 "_15_pKO-F-PBS_TPM", 
                 "_12_pKO-M-5_g_TPM", 
                 "_13_pKO-M-5_g_TPM", 
                 "_16_pKO-F-5_g_TPM", 
                 "_17_pKO-F-5_g_TPM")
# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {
    flag_row_name = T
  } else {
    flag_row_name = F
  }
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_KO - 5_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = "(.*KO.*5_g.*|.*KO.*PBS.*)", 
                            top_n = top_n, 
                            pivot = T)
  df_diff_expr <- data$extract_diff_expr_tpm()

  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_KO - 5_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw_heatmap_std(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - KO - 5_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}
## compare treatment WT - 5_g vs. PBS====
column_level_1 = c("PBS","PBS","LPS_5_ug","LPS_5_ug","LPS_5_ug")
column_level_2 = c("M","F","M","F","F")
column_order = c("_1_pWT-M-PBS_TPM", 
                 "_6_pWT-F-PBS_TPM", 
                 "_3_pWT-M-5_g_TPM", 
                 "_7_pWT-F-5_g_TPM", 
                 "_8_pWT-F-5_g_TPM")

for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {
    flag_row_name = T
  } else {
    flag_row_name = F
  }
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_WT - 5_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = "(.*WT.*5_g.*|.*WT.*PBS.*)", 
                            top_n = top_n, 
                            pivot = T)
  df_diff_expr <- data$extract_diff_expr_tpm()
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_WT - 5_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw_heatmap_std(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - WT - 5_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}

## compare treatment WT - 5_100 vs. PBS====
column_level_1 = c("PBS","PBS","LPS_100_ug","LPS_100_ug","LPS_100_ug")
column_level_2 = c("M","F","M","M","F")
column_order = c("_1_pWT-M-PBS_TPM", 
                 "_6_pWT-F-PBS_TPM", 
                 "_4_pWT-M-100_g_TPM", 
                 "_5_pWT-M-100_g_TPM",
                 "_10_pWT-F-100_g_TPM")

for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {
    flag_row_name = T
  } else {
    flag_row_name = F
  }
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_WT - 100_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = "(.*WT.*100_g.*|.*WT.*PBS.*)", 
                            top_n = top_n, 
                            pivot = T)
  df_diff_expr <- data$extract_diff_expr_tpm()
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_WT - 100_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw_heatmap_std(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - WT - 100_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
                                flag_row_name = flag_row_name
  )
}

## compare treatment KO - 5_100 vs. PBS====
column_level_1 = c( "PBS","PBS","LPS_100_ug","LPS_100_ug")
column_level_2 = c( "M","F","M","F")
column_order = c("_11_pKO-M-PBS_TPM", 
                 "_15_pKO-F-PBS_TPM", 
                 "_14_pKO-M-100_g_TPM", 
                 "_40_pKO-F-100_g_TPM")

for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {
    flag_row_name = T
  } else {
    flag_row_name = F
  }
  # process data
  data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_KO - 100_g vs. PBS.csv"), 
                            gene_folder = gene_folder, 
                            pattern = "(.*KO.*100_g.*|.*KO.*PBS.*)", 
                            top_n = top_n, 
                            pivot = T)
  df_diff_expr <- data$extract_diff_expr_tpm()
  
  # plot
  heatmap_plot <- HeatmapPlot$new(output = file.path(out_dir, sprintf("heatmap_KO - 100_g vs. PBS - top %d.png", top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw_heatmap_std(df = df_diff_expr, 
                                column_order =column_order, 
                                column_level_1 = column_level_1, 
                                column_level_2 = column_level_2, 
                                column_title = sprintf("Gene Expression - KO - 100_g vs. PBS, Sorted by Log2 fold change Top %d", top_n*2),
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
go_plot$draw_go_analysis(genes, top_n, "Go Analysis: WT - 5_g vs. PBS")

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
go_plot$draw_go_analysis(genes, top_n, "Go Analysis: KO - 5_g vs. PBS")
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
go_plot$draw_go_analysis(genes, top_n, "Go Analysis: WT - 100_g vs. PBS")

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
go_plot$draw_go_analysis(genes, top_n, "Go Analysis: KO - 100_g vs. PBS")

# box plot #####################################################################
## KO - 5_g vs. PBS ===============================
# Var
top_n = 6

# process data
data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_KO - 5_g vs. PBS.csv"), 
                          gene_folder = gene_folder, 
                          pattern = "(.*KO.*5_g.*|.*KO.*PBS.*)", 
                          top_n = top_n, 
                          pivot = F)
df_diff_expr_unpivot <- data$extract_diff_expr_tpm()

# plot
box_plot <- BoxPlot$new(file.path(out_dir, "box_plot_KO - 5_g vs. PBS.png"), 
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw_box_plot(df = df_diff_expr_unpivot, 
                       title = "Box Plot KO - 5_g vs. PBS")
## WT - 5_g vs. PBS ===============================
# Var
top_n = 6

# process data
data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_WT - 5_g vs. PBS.csv"), 
                          gene_folder = gene_folder, 
                          pattern = "(.*WT.*5_g.*|.*WT.*PBS.*)", 
                          top_n = top_n, 
                          pivot = F)
df_diff_expr_unpivot <- data$extract_diff_expr_tpm()

# plot
box_plot <- BoxPlot$new(file.path(out_dir, "box_plot_WT - 5_g vs. PBS.png"), 
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw_box_plot(df = df_diff_expr_unpivot, 
                       title = "Box Plot WT - 5_g vs. PBS")

## KO - 100_g vs. PBS ===============================
# Var
top_n = 6

# process data
data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_KO - 100_g vs. PBS.csv"), 
                          gene_folder = gene_folder, 
                          pattern = "(.*KO.*100_g.*|.*KO.*PBS.*)", 
                          top_n = top_n, 
                          pivot = F)
df_diff_expr_unpivot <- data$extract_diff_expr_tpm()

# plot
box_plot <- BoxPlot$new(file.path(out_dir, "box_plot_KO - 100_g vs. PBS.png"), 
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw_box_plot(df = df_diff_expr_unpivot, 
                       title = "Box Plot KO - 100_g vs. PBS")

## WT - 100_g vs. PBS ===============================
# Var
top_n = 6

# process data
data <- DataProcessor$new(diff_expr_file = file.path(data_dir,"diff_expr_WT - 100_g vs. PBS.csv"), 
                          gene_folder = gene_folder, 
                          pattern = "(.*WT.*100_g.*|.*WT.*PBS.*)", 
                          top_n = top_n, 
                          pivot = F)
df_diff_expr_unpivot <- data$extract_diff_expr_tpm()

# plot
box_plot <- BoxPlot$new(file.path(out_dir, "box_plot_WT - 100_g vs. PBS.png"), 
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw_box_plot(df = df_diff_expr_unpivot, 
                       title = "Box Plot WT - 100_g vs. PBS")
