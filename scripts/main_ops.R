source("scripts/install.R")
source("scripts/draw.R")
gene_folder = "data/LPS/genes"
data_dir = "data/LPS"
out_dir = "out/LPS"

# heat map #####################################################################
## compare treatment KO - 100_g vs. PBS====
# var
diff_expr_file = "diff_expr_KO,WT,SPURT,TG - 100_g vs. PBS.csv" # Change as you need
pattern = "(.*KO.*100_g.*|.*KO.*PBS.*|.*WT.*100_g.*|.*WT.*PBS.*|.*TG.*100_g.*|.*TG.*PBS.*|.*SPURT.*100_g.*|.*SPURT.*PBS.*)" # Change as you need

control_group = str_remove_all(diff_expr_file, "diff_expr_|\\.csv")

# top 2Ns
for (top_n in c(30, 0)) {
  if (top_n < 50 & top_n > 0) {flag_row_name = T} else {flag_row_name = F}
 
  # plot
  heatmap_plot <- HeatmapPlot$new(data_dir = data_dir,
                                  gene_folder = gene_folder,
                                  output = file.path(out_dir, sprintf("heatmap_%s - top %d.png", control_group, top_n*2)), 
                                  width = 6000, 
                                  height = 5000, 
                                  res = 300, 
                                  bg = "white")
  heatmap_plot$draw(diff_expr_file = diff_expr_file,
                    column_title = sprintf("Gene Expression - %s, Sorted by Log2 fold change Top %d", control_group, top_n*2),
                    flag_row_name = flag_row_name,
                    top_n = top_n,
                    pattern = pattern)
}

# go analysis ##################################################################
## WT - 5_g vs. PBS =============
# var
diff_expr_file = "diff_expr_KO,WT,SPURT - 100_g vs. PBS.csv" # Change as you need
top_n = 10 # Change as you need


control_group = str_remove_all(diff_expr_file, "diff_expr_|\\.csv")

# plot
go_plot <- GoAnalysisPlot$new(data_dir = data_dir,
                              gene_folder = gene_folder,
                              output = file.path(out_dir, sprintf("go_enrichment_%s - top %d.png", control_group, top_n)), 
                              width = 2000, 
                              height = 4000, 
                              res = 300, 
                              bg = "white")
go_plot$drawDotPlot(diff_expr_file = diff_expr_file, 
             top_n = top_n, 
             title = sprintf("Go Analysis: %s - top %d", control_group, top_n))
# KEGG analysis ##################################################################
## WT - 5_g vs. PBS =============
# var

diff_expr_file = "diff_expr_KO,WT,SPURT - 100_g vs. PBS.csv" # Change as you need
top_n = 15 # Change as you need


control_group = str_remove_all(diff_expr_file, "diff_expr_|\\.csv")

# plot
kegg_plot <- KEGGPlot$new(data_dir = data_dir,
                              gene_folder = gene_folder,
                              output = file.path(out_dir, sprintf("kegg_enrichment_dot_plot_%s - top %d.png", control_group, top_n)), 
                              width = 2000, 
                              height = 4000, 
                              res = 300, 
                              bg = "white")
kegg_plot$drawDotPlot(diff_expr_file = diff_expr_file, 
             top_n = top_n, 
             title = sprintf("KEGG Analysis: %s - top %d", control_group, top_n))

# KEGG analysis - CNetPlot ##################################################################
## WT - 5_g vs. PBS =============
# var

diff_expr_file = "diff_expr_KO,WT,SPURT - 100_g vs. PBS.csv" # Change as you need
top_n = 5 # Change as you need


control_group = str_remove_all(diff_expr_file, "diff_expr_|\\.csv")

# plot
kegg_plot <- KEGGPlot$new(data_dir = data_dir,
                          gene_folder = gene_folder,
                          output = file.path(out_dir, sprintf("kegg_enrichment_cnet_plot_%s - top %d.png", control_group, top_n)), 
                          width = 6000, 
                          height = 4000, 
                          res = 300, 
                          bg = "white")
kegg_plot$drawCNetPlot(diff_expr_file = diff_expr_file, 
                       top_n = top_n, 
                       title = sprintf("KEGG Analysis: %s - top %d", control_group, top_n))
# box plot #####################################################################
## KO - 5_g vs. PBS ===============================
# Var
diff_expr_file = "diff_expr_WT - 5_g vs. PBS.csv" # Change as you need
pattern = "(.*WT.*5_g.*|.*WT.*PBS.*)" # Change as you need
top_n = 6 # Change as you need

control_group = str_remove_all(diff_expr_file, "diff_expr_|\\.csv")

# plot
box_plot <- BoxPlot$new(data_dir = data_dir,
                        gene_folder = gene_folder,
                        output = file.path(out_dir, sprintf("box_plot_%s.png", control_group)),
                        width = 3000, 
                        height = 2400, 
                        res = 300, 
                        bg = "white")
box_plot$draw(diff_expr_file = diff_expr_file, 
              title = sprintf("Box Plot %s", control_group),
              pattern = pattern)

