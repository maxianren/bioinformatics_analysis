source("scripts/install.R")
source("scripts/draw.R")
gene_folder = "data/LPS/genes"
data_dir = "data/LPS"
out_dir = "out/LPS"

# heat map #####################################################################
## compare treatment KO - 100_g vs. PBS====
# var
diff_expr_file = "diff_expr_KO - 100_g vs. PBS.csv" # Change as you need
pattern = "(.*KO.*100_g.*|.*KO.*PBS.*)" # Change as you need

feature = "TPM"
control_group = str_remove_all(diff_expr_file, "diff_expr_|\\.csv")

# top 2Ns
for (top_n in c(30, 100, 250)) {
  if (top_n < 50) {flag_row_name = T} else {flag_row_name = F}
 
  # plot
  heatmap_plot <- HeatmapPlot$new(data_dir = data_dir,
                                  gene_folder = gene_folder,
                                  output = file.path(out_dir, sprintf("heatmap_%s - top %d.png", control_group, top_n*2)), 
                                  width = 6000, 
                                  height = 4800, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(diff_expr_file = diff_expr_file,
                    column_title = sprintf("Gene Expression - %s, Sorted by Log2 fold change Top %d", control_group, top_n*2),
                    flag_row_name = flag_row_name,
                    pattern = pattern,
                    feature = feature)
}

# go analysis ##################################################################
## WT - 5_g vs. PBS =============
# var
diff_expr_file = "diff_expr_WT - 5_g vs. PBS.csv" # Change as you need
log_2_FC = 1 # Change as you need
top_n = 10 # Change as you need

control_group = str_remove_all(diff_expr_file, "diff_expr_|\\.csv")

# plot
go_plot <- GoAnalysisPlot$new(data_dir = data_dir,
                              gene_folder = gene_folder,
                              output = file.path(out_dir, sprintf("go_enrichment_%s.png", control_group)), 
                              width = 2000, 
                              height = 4000, 
                              res = 300, 
                              bg = "white")
go_plot$draw(diff_expr_file = diff_expr_file, 
             top_n = top_n, 
             title = sprintf("Go Analysis: %s", control_group))


# box plot #####################################################################

## KO - 5_g vs. PBS ===============================
# Var
diff_expr_file = "diff_expr_WT - 5_g vs. PBS.csv" # Change as you need
pattern = "(.*WT.*5_g.*|.*WT.*PBS.*)" # Change as you need
top_n = 6 # Change as you need

feature = "TPM"
control_group = str_remove_all(diff_expr_file, "diff_expr_|\\.csv")

# plot
box_plot <- BoxPlot$new(data_dir = data_dir,
                        gene_folder = gene_folder,
                        output = file.path(out_dir, sprintf("box_plot_KO - 5_g vs. PBS.png", control_group)),
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw(diff_expr_file = diff_expr_file, 
              title = sprintf("Box Plot %s", control_group),
              pattern = pattern,
              feature = feature)

