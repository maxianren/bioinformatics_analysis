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
                          output = NULL,
                          width = NULL,
                          height = NULL,
                          res = NULL,
                          bg = NULL,
                          
                          initialize = function(output, width, height, res, bg) {
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
                         draw_heatmap_std = function(df, column_order, column_level_1, column_level_2, column_title, flag_row_name) {
                           data_matrix <- df %>%
                             dplyr::select(all_of(column_order)) %>%
                             as.matrix()
                           
                           rownames(data_matrix) <- df$Name
                           
                           column_annotation <- HeatmapAnnotation(
                             Treatment = column_level_1,
                             Sex = column_level_2,
                             annotation_name_side = "right"
                           )
                           
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
                           
                           self$save_plot(p)
                         }
                       )
)

GoAnalysisPlot <- R6Class("GoAnalysisPlot",
                          inherit = PlottingBase,
                          
                          public = list(
                            draw_go_analysis = function(genes, top_n, title) {
                              GO_results <- enrichGO(gene = genes, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = 'ALL')
                              plot <- dotplot(GO_results, showCategory = top_n, split = "ONTOLOGY", title = title) + facet_grid(ONTOLOGY ~ ., scale = "free")
                              
                              self$save_plot(plot)
                            }
                          )
)

BoxPlot <- R6Class("BoxPlot",
                   inherit = PlottingBase,
                   
                   public = list(
                     draw_box_plot = function(df, title) {
                       plot <- ggplot(df, aes(x = Name, y = TPM, fill = Treatment)) +
                         geom_boxplot() +
                         labs(title = title, x = "Gene", y = "TPM") +
                         theme_minimal()
                       
                       self$save_plot(plot)
                     }
                   )
)