source("scripts/format_data.R")
library(stringr)
library(R6)
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

PlottingBase <- R6Class("PlottingBase",
                        public = list(
                          data_dir = NULL,
                          gene_folder = NULL,
                          output = NULL,
                          width = NULL,
                          height = NULL,
                          res = NULL,
                          bg = NULL,

                          initialize = function(data_dir, gene_folder, output, width, height, res, bg) {
                            self$data_dir <- data_dir
                            self$gene_folder <- gene_folder
                            self$output <- output
                            self$width <- width
                            self$height <- height
                            self$res <- res
                            self$bg <- bg
                          },

                          save_plot = function(plot) {
                            png(self$output, width = self$width, height = self$height, res = self$res, bg = self$bg)
                            print(plot)
                            dev.off()
                          }
                        )
)

HeatmapPlot <- R6Class("HeatmapPlot",
                       inherit = PlottingBase,

                       public = list(
                         draw = function(diff_expr_file, column_title, flag_row_name, pattern, feature) {
                           # var
                           df_gene <- GeneFeatureExtractor$new(self$gene_folder, pattern, feature)$get_feature_hierarchy()

                           column_order <- df_gene$file_name_clean
                           column_level_1 <- df_gene$treatment
                           column_level_2 <- df_gene$genotype
                           column_level_3 <- df_gene$gender

                           column_annotation <- HeatmapAnnotation(
                             Treatment = column_level_1,
                             Genotype = column_level_2,
                             Gender = column_level_3,
                             annotation_name_side = "right"
                           )
                           # get dataframe
                           data <- DataProcessor$new(diff_expr_file = file.path(self$data_dir,diff_expr_file),
                                                     gene_folder = gene_folder,
                                                     pattern = pattern)
                           df_diff_expr <- data$extract_top_n_gene_feature(top_n = top_n,
                                                                           pivot = T,
                                                                           feature = feature)
                           # plot
                           data_matrix <- df_diff_expr %>%
                             dplyr::select(all_of(column_order)) %>%
                             as.matrix()

                           rownames(data_matrix) <- df_diff_expr$Name

                           z_score_matrix <- t(scale(t(data_matrix)))

                           p <- Heatmap(
                             z_score_matrix,
                             name = "Expression",
                             col = colorRampPalette(c("blue", "white", "red"))(75),
                             cluster_rows = TRUE,
                             cluster_columns = FALSE,
                             show_row_names = flag_row_name,
                             top_annotation = column_annotation,
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
                           # save the image
                           self$save_plot(p)
                         }
                       )
)


GoAnalysisPlot <- R6Class("GoAnalysisPlot",
                          inherit = PlottingBase,
                          
                          public = list(
                            draw = function(diff_expr_file, top_n, title) {
                              # get dataframe
                              data = DifferentialExpressionDataProcessor$new(file_path = file.path(data_dir,diff_expr_file))
                              genes = data$filter_by_log2_fc(log_2_FC)
                              
                              # plot
                              GO_results <- enrichGO(gene = genes, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = 'ALL')
                              plot <- dotplot(GO_results, showCategory = top_n, split = "ONTOLOGY", title = title) + facet_grid(ONTOLOGY ~ ., scale = "free")

                              # save the image
                              self$save_plot(plot)
                            }
                          )
)

BoxPlot <- R6Class("BoxPlot",
                   inherit = PlottingBase,
                   
                   public = list(
                     draw = function(diff_expr_file, title, pattern, feature) {
                       # process data
                       data <- DataProcessor$new(diff_expr_file = file.path(self$data_dir,diff_expr_file), 
                                                 gene_folder = self$gene_folder, 
                                                 pattern = pattern)
                       df_diff_expr_unpivot <- data$extract_top_n_gene_feature(top_n = top_n, 
                                                                               pivot = F,
                                                                               feature = feature)
                       # plot
                       plot <- ggplot(df_diff_expr_unpivot, aes(x = Name, y = TPM, fill = Treatment)) +
                         geom_boxplot() +
                         labs(title = title, x = "Gene", y = "TPM") +
                         theme_minimal()

                       # save the image
                       self$save_plot(plot)
                     }
                   )
)